from typing import Tuple
import pandas as pd
import numpy as np
import os
import io
import time
import zipfile
import itertools
import pathlib
from utils.utils import dbg
from utils import utils, unique_id_generator
import urllib.request, urllib.error, urllib.parse

DBFILE_NAME = "cellphonedb"
MULTIDATA_TABLE_BOOLEAN_COLS = ['receptor','other','secreted_highlight',\
                                'transmembrane','secreted','peripheral','integrin','is_complex']
INPUT_FILE_NAMES = ['complex_input','gene_input','interaction_input','protein_input']

def get_interactions_genes_complex(user_dir_root, db_version) -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    db_path = os.path.join(user_dir_root, "releases", db_version)
    # Extract csv files from db_files_path/cellphonedb.zip into dbTableDFs
    dbTableDFs = extract_dataframes_from_db(db_path)
    # Convert dbTableDFs into interactions, genes, complex_composition, complex_expanded data frames
    mtTable = dbTableDFs['multidata_table']
    dbg(mtTable.dtypes)
    # Convert all MULTIDATA_TABLE_BOOLEAN_COLS from Integer (0/1) to Boolean
    for col in MULTIDATA_TABLE_BOOLEAN_COLS:
        mtTable[col] = mtTable[col].astype(bool)
    # Read genes 'table' - c.f. old CellphoneDB: GeneRepository.get_all_expanded()
    # Drop 'protein_name' column from dbTableDFs['gene_table'] as dbTableDFs['protein_table'] already has it
    dbTableDFs['gene_table'] = dbTableDFs['gene_table'].drop('protein_name', axis=1)
    genes = pd.merge(dbTableDFs['gene_table'], dbTableDFs['protein_table'], left_on='protein_id', right_on='id_protein')
    genes = pd.merge(genes, mtTable, left_on='protein_multidata_id', right_on='id_multidata')
    dbg("genes columns: ", genes.columns)
    # Read complex compositions
    complex_composition = dbTableDFs['complex_composition_table']
    # Read 'interactions' - c.f. old CellphoneDB: InteractionRepository.get_all_expanded()
    suffixes = ('_1', '_2')
    includeGene = False
    if includeGene:
        multidata_simple = genes.copy()
    else:
        multidata_simple = pd.merge(dbTableDFs['protein_table'], mtTable, left_on='protein_multidata_id', right_on='id_multidata')
    multidata_complex = pd.merge(mtTable, dbTableDFs['complex_table'], left_on='id_multidata', right_on='complex_multidata_id')
    if multidata_complex.empty:
        multidata_expanded = multidata_simple
    else:
        multidata_expanded = pd.concat([multidata_simple, multidata_complex], ignore_index=True, sort=True)
    # C.f. old CellphoneDB: InteractionRepository.get_all_expanded()
    multidata_expanded = multidata_expanded.astype({'id_multidata': 'int64'})
    dbg(multidata_expanded.columns)
    interactions = pd.merge(dbTableDFs['interaction_table'], multidata_expanded, left_on=['multidata_1_id'], right_on=['id_multidata'])
    interactions = pd.merge(interactions, multidata_expanded, left_on=['multidata_2_id'], right_on=['id_multidata'], suffixes=suffixes)
    dbg("interactions columns: ", interactions.columns)
    # Generate complex_expanded - c.f. ComplexRepository.get_all_expanded()
    complex_expanded = pd.merge(dbTableDFs['complex_table'], mtTable, left_on='complex_multidata_id', right_on='id_multidata')
    # dbg(complexes_expanded.columns)
    # index interactions and complex data frames
    # C.f. old CellphoneDB: method_launcher.get_interactions_genes_complex()
    interactions.set_index('id_interaction', drop=True, inplace=True)
    complex_composition.set_index('id_complex_composition', inplace=True, drop=True)

    return interactions, genes, complex_composition, complex_expanded

def extract_dataframes_from_db(db_files_path):
    dfs = {}
    start = time.time()
    for tuple in unzip(os.path.join(db_files_path,'{}.{}'.format(DBFILE_NAME, "zip"))):
        file_name = tuple[0]
        file_handle = tuple[1]
        dbg("Retrieving from zip file: " + file_name)
        dfs[file_name.replace('.csv','')] = pd.read_csv(file_handle)
    duration = time.time() - start
    dbg("Loaded DB into memory in " + str(round(duration,2)) + "s")
    return dfs

def unzip(zip_file_path):
    """
    Download a ZIP file and extract its contents in memory
    yields (filename, file-like object) pairs
    """
    with open(zip_file_path,'br') as file:
        # Note - eval unescapes the double-quotes
        content = file.read()
        dbg(type(content))
        with zipfile.ZipFile(io.BytesIO(content)) as thezip:
            for zipinfo in thezip.infolist():
                with thezip.open(zipinfo) as thefile:
                    yield zipinfo.filename, thefile

def get_db_path(user_dir_root, db_version):
    return os.path.join(user_dir_root, "releases", db_version)

def create_db(user_dir_root, db_version, \
              gene_input=None, protein_input=None, complex_input=None, interaction_input=None):
    db_path = os.path.join(user_dir_root, "releases", db_version)
    pathlib.Path(db_path).mkdir(parents=True, exist_ok=True)
    dataDFs = getDFs(gene_input=gene_input, protein_input=protein_input, complex_input=complex_input,
                     interaction_input=interaction_input)

    # Collect protein data
    protein_db_df = dataDFs['protein_input'][['protein_name', 'tags', 'tags_reason', 'tags_description', 'uniprot']]
    num_proteins = protein_db_df.shape[0]
    multidata_id_list_so_far = list(range(num_proteins))
    protein_db_df.insert(0, 'id_protein', multidata_id_list_so_far, False)
    protein_db_df.insert(len(protein_db_df.columns), 'protein_multidata_id', protein_db_df['id_protein'].tolist(), True)
    # dbg(protein_db_df.info)

    # Collect gene data
    gene_db_df = dataDFs['gene_input'][['ensembl', 'gene_name', 'hgnc_symbol', 'uniprot']]
    num_genes = gene_db_df.shape[0]
    gene_db_df.insert(0, 'id_gene', list(range(num_genes)), False)
    # Assign values from protein_db_df['protein_multidata_id'] into gene_db_df['protein_id']
    # via join between 'uniprot' and 'protein_name'
    gene_db_df = pd.merge(gene_db_df, protein_db_df[['protein_name', 'protein_multidata_id', 'uniprot']], on='uniprot')
    gene_db_df = gene_db_df.drop('uniprot', axis=1)
    protein_db_df = protein_db_df.drop('uniprot', axis=1)
    gene_db_df.rename(columns = {'protein_multidata_id':'protein_id'}, inplace = True)
    # print(gene_db_df.info)

    # Collect complex data
    complex_db_df = dataDFs['complex_input'] \
        [['uniprot_1','uniprot_2','uniprot_3','uniprot_4','pdb_structure','pdb_id','stoichiometry','comments_complex']]
    # Note that uniprot_* cols will be dropped after complex_composition_df has been constructed
    num_complexes = complex_db_df.shape[0]
    complex_db_df.insert(0, 'id_complex', list(range(num_complexes)), False)
    next_md_id = multidata_id_list_so_far[-1] + 1
    complex_multidata_ids = list(range(next_md_id, next_md_id + num_complexes))
    complex_db_df.insert(1, 'complex_multidata_id', complex_multidata_ids, False)
    multidata_id_list_so_far.extend(complex_multidata_ids)
    # dbg(multidata_id_list_so_far[-1])
    # dbg(complex_db_df.info)

    # Collect multidata
    # Insert proteins into multidata
    multidata_db_df = \
         dataDFs['protein_input'][['uniprot','receptor','receptor_desc','other','other_desc','secreted_highlight','secreted_desc','transmembrane','secreted','peripheral','integrin']].copy()
    multidata_db_df.rename(columns={'uniprot':'name'}, inplace=True)
    multidata_ids = pd.merge(dataDFs['protein_input'][['protein_name']], \
                             protein_db_df[['protein_name', 'protein_multidata_id']], on='protein_name')['protein_multidata_id'].tolist()
    multidata_db_df.insert(0, 'id_multidata', multidata_ids, False)
    multidata_db_df.insert(len(multidata_db_df.columns), 'is_complex', list(itertools.repeat(False, multidata_db_df.shape[0])), True)
    dbg(multidata_db_df.shape, multidata_db_df.index, multidata_db_df.columns)
    # Insert complexes into multidata
    complex_aux_df = \
         dataDFs['complex_input'][['complex_name','receptor','receptor_desc','other','other_desc','secreted_highlight','secreted_desc','transmembrane','secreted','peripheral','integrin']].copy()
    complex_aux_df.rename(columns={'complex_name': 'name'}, inplace=True)
    complex_aux_df.insert(0, 'id_multidata', complex_multidata_ids, False)
    complex_aux_df.insert(len(complex_aux_df.columns), 'is_complex', list(itertools.repeat(True, complex_aux_df.shape[0])), True)
    dbg(complex_aux_df.shape, complex_aux_df.index, complex_aux_df.columns)
    # Append complex_aux_df to multidata_db_df
    multidata_db_df = pd.concat([multidata_db_df, complex_aux_df], ignore_index=True, verify_integrity=True)
    dbg(multidata_db_df.shape, multidata_db_df.index, multidata_db_df.columns)

    # First collect total_protein counts for each complex in complex_db_df
    uniprot_cols = ['uniprot_1','uniprot_2','uniprot_3','uniprot_4']
    total_protein_cnt_list = np.apply_along_axis(lambda s: sum(type(x) == str for x in s), 1, complex_db_df[uniprot_cols].values).tolist()
    complex_db_df.insert(len(complex_db_df.columns), 'total_protein', total_protein_cnt_list, True)
    dbg(complex_db_df.info)
    # Next collect all complex_composition data into cc_list
    cc_list = []
    for r in complex_db_df[uniprot_cols + ['complex_multidata_id','total_protein']].values.tolist():
        for acc in filter(lambda x: type(x) == str, r):
            protein_multidata_id = \
                multidata_db_df.loc[(multidata_db_df['is_complex'] == False) & (multidata_db_df['name'] == acc), ['id_multidata']] \
                    .iat[0,0]
            complex_multidata_id = r[4]
            total_protein = r[5]
            cc_list.append([complex_multidata_id, protein_multidata_id, total_protein])

    complex_composition_df = pd.DataFrame(cc_list, columns=['complex_multidata_id', 'protein_multidata_id', 'total_protein'])
    complex_composition_df.insert(0, 'id_complex_composition', list(range(len(cc_list))), False)
    dbg(complex_composition_df.shape, complex_composition_df.index, complex_composition_df.columns, complex_composition_df.info)
    # Next drop the auxiliary columns from complex_db_df: uniprot_cols and 'total_protein'
    for col in uniprot_cols + ['total_protein']:
        complex_db_df = complex_db_df.drop(col, axis=1)

    # Collect interaction data
    interactions_aux_df = pd.merge(dataDFs['interaction_input'], multidata_db_df, \
                                   left_on=['partner_a'],right_on=['name'])
    interactions_aux_df = pd.merge(interactions_aux_df, multidata_db_df, \
                                               left_on=['partner_b'], right_on=['name'], \
                                               suffixes=['_x', '_y'])
    dbg(interactions_aux_df.shape)
    # Remove interactions non-CPDB interactors
    interactions_aux_df = interactions_aux_df[
        interactions_aux_df.apply(
            lambda interaction: utils.is_cellphonedb_interactor(interaction, ('_x', '_y')), axis=1)]
    interactions_aux_df['id_cp_interaction'] = interactions_aux_df.apply(
        lambda interaction: unique_id_generator.interaction(interaction, ('_x', '_y')), axis=1)
    dbg(interactions_aux_df.info)
    dbg(interactions_aux_df.columns)
    interactions_df = interactions_aux_df[['id_cp_interaction','id_multidata_x','id_multidata_y', \
                                          'source','annotation_strategy']].copy()
    interactions_df.rename(columns={'id_multidata_x': 'multidata_1_id', 'id_multidata_y': 'multidata_2_id'}, inplace=True)
    interactions_df.insert(0, 'id_interaction', list(range(interactions_df.shape[0])), False)
    dbg(interactions_df.shape, interactions_df.index, interactions_df.columns)

    # Save all DFs as csv files inside a DB zip file
    zip_buffer = io.BytesIO()
    with zipfile.ZipFile(zip_buffer, "a",
                         zipfile.ZIP_DEFLATED, False) as zip_file:
        zip_file.writestr('protein_table.csv', protein_db_df.to_csv(index=False, sep=',').encode('utf-8'))
        zip_file.writestr('gene_table.csv', gene_db_df.to_csv(index=False, sep=',').encode('utf-8'))
        zip_file.writestr('complex_table.csv', complex_db_df.to_csv(index=False, sep=',').encode('utf-8'))
        zip_file.writestr('complex_composition_table.csv', complex_composition_df.to_csv(index=False, sep=',').encode('utf-8'))
        zip_file.writestr('multidata_table.csv', multidata_db_df.to_csv(index=False, sep=',').encode('utf-8'))
        zip_file.writestr('interaction_table.csv', interactions_df.to_csv(index=False, sep=',').encode('utf-8'))
    with open(os.path.join(db_path,'cellphonedb.zip'), 'wb') as f:
        f.write(zip_buffer.getvalue())
    print("Created CellphoneDB {} in {} successfully" \
          .format(db_version, os.path.join(db_path,'cellphonedb.zip')))

def download_input_files(user_dir_root, db_version):
    db_path = get_db_path(user_dir_root, db_version)
    db_data_path = os.path.join(db_path, "data")
    pathlib.Path(db_path).mkdir(parents=True, exist_ok=True)
    for fname in INPUT_FILE_NAMES:
        fname = fname+".csv"
        print("Downloading: " + fname)
        r = urllib.request.urlopen('https://raw.githubusercontent.com/ventolab/cellphonedb-data/master/data/' + fname)
        with open(os.path.join(db_data_path, fname), 'wb') as f:
            f.write(r.read())

def getDFs(gene_input=None, protein_input=None, complex_input=None, interaction_input=None):
    dfs = {}
    dfs['gene_input'] = utils.read_data_table_from_file(gene_input)
    dfs['protein_input'] = utils.read_data_table_from_file(protein_input)
    dfs['complex_input'] = utils.read_data_table_from_file(complex_input)
    dfs['interaction_input'] = utils.read_data_table_from_file(interaction_input)
    return dfs