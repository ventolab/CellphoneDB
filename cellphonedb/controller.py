import sys
import pandas as pd
import numpy as np
import anndata
import os
from typing import Tuple
import io
import zipfile
import time
import itertools
from utils import utils, unique_id_generator, generate_input_files
import urllib.request, urllib.error, urllib.parse
import pathlib
import json

from src.core.methods import cpdb_analysis_method, cpdb_statistical_analysis_method, cpdb_degs_analysis_method
from src.core.preprocessors import method_preprocessors
from src.exceptions.ParseCountsException import ParseCountsException

MULTIDATA_TABLE_BOOLEAN_COLS = ['receptor','other','secreted_highlight','transmembrane','secreted','peripheral','integrin','is_complex']
DBFILE_NAME = "cellphonedb"
INPUT_FILE_NAMES = ['complex_input','gene_input','interaction_input','protein_input']
KEY2USER_TEST_FILE = {'counts' : 'test_counts.txt', 'meta': 'test_meta.txt', \
                         'microenvs' : 'test_microenviroments.txt', 'degs' : 'test_degs.txt'}
DEBUG=False

RELEASED_VERSION="v5.0.0"
CPDB_ROOT = os.path.join(os.path.expanduser('~'),".cpdb")

def get_db_path(user_dir_root, db_version):
    return os.path.join(user_dir_root, "releases", db_version, "data")

def download_input_files(db_files_path):
    for fname in INPUT_FILE_NAMES:
        fname = fname+".csv"
        print("Downloading: " + fname)
        r = urllib.request.urlopen('https://raw.githubusercontent.com/ventolab/cellphonedb-data/master/data/' + fname)
        with open(os.path.join(db_files_path,fname), 'wb') as f:
            f.write(r.read())

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

def getDFs(path, files_list, file_extension):
    dfs = {}
    for fname in files_list:
        dfs[fname] = utils.read_data_table_from_file(os.path.join(path,'{}.{}'.format(fname, file_extension)))
    return dfs

def extract_dataframes_from_db(db_files_path):
    dfs = {}
    start = time.time()
    for tuple in unzip(os.path.join(db_files_path,'{}.{}'.format(DBFILE_NAME, "zip"))):
        file_name = tuple[0]
        file_handle = tuple[1]
        dbg("Retrieving from zip file: " + file_name)
        dfs[file_name.replace('.csv','')] = pd.read_csv(file_handle)
    duration = time.time() - start
    dbg("Retrieved from DB zip file CSV files as DataFrames in: " + str(round(duration,2)) + "s")
    return dfs

def dbg(*argv):
    if DEBUG:
        for arg in argv:
            print(arg)

def get_interactions_genes_complex(user_dir_root, db_version) -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    db_files_path = get_db_path(user_dir_root, db_version)
    # Extract csv files from db_files_path/cellphonedb.zip into dbTableDFs
    dbTableDFs = extract_dataframes_from_db(db_files_path)
    # Convert dbTableDFs into interactions, genes, complex_composition, complex_expanded data frames
    mtTable = dbTableDFs['multidata_table']
    dbg(mtTable.dtypes)
    # Convert all MULTIDATA_TABLE_BOOLEAN_COLS from Integer (0/1) to Boolean
    for col in MULTIDATA_TABLE_BOOLEAN_COLS:
        mtTable[col] = mtTable[col].astype(bool)
    # Read genes 'table' - c.f. old CellphoneDB: GeneRepository.get_all_expanded()
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
        multidata_expanded = multidata_simple.append(multidata_complex, ignore_index=True, sort=True)
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

def get_user_files(user_dir_root, \
        counts_fn=KEY2USER_TEST_FILE['counts'], meta_fn=KEY2USER_TEST_FILE['meta'], microenvs_fn=None, degs_fn=None) \
        -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    loaded_user_files=[]
    user_files_path = os.path.join(user_dir_root,"user_files")
    # Read user files
    counts = utils.read_data_table_from_file(os.path.join(user_files_path, counts_fn),
                                             index_column_first=True)
    loaded_user_files.append(KEY2USER_TEST_FILE['counts'])
    raw_meta = utils.read_data_table_from_file(os.path.join(user_files_path, meta_fn),
                                               index_column_first=False)
    meta = method_preprocessors.meta_preprocessor(raw_meta)
    loaded_user_files.append(KEY2USER_TEST_FILE['meta'])

    if microenvs_fn:
        microenvs = utils.read_data_table_from_file(os.path.join(user_files_path, microenvs_fn))
        loaded_user_files.append(KEY2USER_TEST_FILE['microenvs'])
    else:
        microenvs = pd.DataFrame()

    if degs_fn:
        degs = utils.read_data_table_from_file(os.path.join(user_files_path, degs_fn))
        loaded_user_files.append(KEY2USER_TEST_FILE['degs'])
    else:
        degs = pd.DataFrame()

    print("The following user files were loaded successfully:")
    for fn in loaded_user_files:
        print(fn)

    return counts, raw_meta, meta, microenvs, degs

def get_user_file(user_dir_root, h5ad_fn='test.h5ad') \
        -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    adata = utils.read_h5ad(os.path.join(user_dir_root,'user_files',h5ad_fn),)
    counts = adata.to_df().T
    counts.columns = adata.obs['sample']

    meta_str = adata.uns['meta']
    if meta_str:
        raw_meta = utils._read_data(io.StringIO(meta_str), separator='\t', index_column_first=False, dtype=None,
                                    na_values=None, compression=None)
        meta = method_preprocessors.meta_preprocessor(raw_meta)
        dbg(meta.info)

    microenvs_str = adata.uns['microenvs']
    if microenvs_str:
        microenvs = utils._read_data(io.StringIO(microenvs_str), separator='\t', index_column_first=False, dtype=None,
                                     na_values=None, compression=None)
        dbg(microenvs.info)
    else:
        microenvs = pd.DataFrame()

    degs_str = adata.uns['degs']
    if degs_str:
        degs = utils._read_data(io.StringIO(degs_str), separator='\t', index_column_first=False, dtype=None, na_values=None,
                                compression=None)
        dbg(degs.info)
    else:
        degs = pd.DataFrame()

    print("User file {} was loaded successfully".format(h5ad_fn))
    return counts, raw_meta, meta, microenvs, degs

def _counts_validations(counts: pd.DataFrame, meta: pd.DataFrame) -> pd.DataFrame:
    if not len(counts.columns):
        raise ParseCountsException('Counts values are not decimal values', 'Incorrect file format')
    try:
        if np.any(counts.dtypes.values != np.dtype('float32')):
            counts = counts.astype(np.float32)  # type: pd.DataFrame
    except:
        raise ParseCountsException

    meta.index = meta.index.astype(str)

    if np.any(~meta.index.isin(counts.columns)):
        raise ParseCountsException("Some cells in meta did not exist in counts",
                                   "Maybe incorrect file format")

    if np.any(~counts.columns.isin(meta.index)):
        core_logger.debug("Dropping counts cells that are not present in meta")
        counts = counts.loc[:, counts.columns.isin(meta.index)].copy()
    return counts

def testrun_analyses(user_dir_root, db_version):
    interactions, genes, complex_composition, complex_expanded = get_interactions_genes_complex(user_dir_root, db_version)
    counts, raw_meta, meta, microenvs, degs = get_user_files(user_dir_root, \
                                                                        counts_fn=KEY2USER_TEST_FILE['counts'],
                                                                        meta_fn=KEY2USER_TEST_FILE['meta'], \
                                                                        microenvs_fn=KEY2USER_TEST_FILE['microenvs'],
                                                                        degs_fn=None)
    # ************ Call analysis method
    means, significant_means, deconvoluted = cpdb_analysis_method.call(
        meta,
        counts,
        'ensembl',
        interactions,
        genes,
        complex_expanded,
        complex_composition,
        microenvs=microenvs,
        # Does not store results files in cellphonedb/out/debug_intermediate.pkl
        debug=False,
        output_path=os.path.join('.','out')
        );
    # dbg("means:", means.index, means.columns, means.info)
    # dbg("significant_means:", significant_means.index, significant_means.columns, significant_means.info)
    # dbg("deconvoluted:", deconvoluted.index, deconvoluted.columns, deconvoluted.info)

    # ************ Call statistical analysis method
    meta = method_preprocessors.meta_preprocessor(raw_meta)
    counts = _counts_validations(counts, meta)
    subsampler = None
    if subsampler is not None:
        counts = subsampler.subsample(counts)
    deconvoluted, means, pvalues, significant_means = \
        cpdb_statistical_analysis_method.call(meta,
                                              counts,
                                              'ensembl',
                                              interactions,
                                              genes,
                                              complex_expanded,
                                              complex_composition,
                                              microenvs=microenvs,
                                              iterations = 10,
                                              threshold = 0.05,
                                              threads = 4,
                                              debug_seed = -1,
                                              result_precision = 5,
                                              pvalue = 0.05,
                                              separator = '|',
                                              debug = False,
                                              output_path = '')

    # dbg("means:", means.index, means.columns, means.info)
    # dbg("significant_means:", significant_means.index, significant_means.columns, significant_means.info)
    # dbg("deconvoluted:", deconvoluted.index, deconvoluted.columns, deconvoluted.info)
    # dbg("pvalues:", pvalues.index, pvalues.columns, pvalues.info)

    # ************ Call degs analysis method
    # NB. Same data prep as for cpdb_statistical_analysis_method
    deconvoluted, means, relevant_interactions, significant_means = \
        cpdb_degs_analysis_method.call(meta,
                                    counts,
                                    degs,
                                    'ensembl',
                                    interactions,
                                    genes,
                                    complex_expanded,
                                    complex_composition,
                                    microenvs=microenvs,
                                    iterations = 10,
                                    threshold = 0.05,
                                    threads = 4,
                                    debug_seed = -1,
                                    result_precision = 5,
                                    separator = '|',
                                    debug = False,
                                    output_path = '')

    # dbg("means:", means.index, means.columns, means.info)
    # dbg("significant_means:", significant_means.index, significant_means.columns, significant_means.info)
    # dbg("deconvoluted:", deconvoluted.index, deconvoluted.columns, deconvoluted.info)
    # dbg("relevant_interactions:", relevant_interactions.index, relevant_interactions.columns, relevant_interactions.info)


def download_source_files(user_dir_root, db_version):
    sources_path = os.path.join(get_db_path(user_dir_root, db_version), "sources")
    print("Downloading cellphonedb-data/data/sources files into {}:".format(sources_path))
    pathlib.Path(sources_path).mkdir(parents=True, exist_ok=True)
    r = urllib.request.urlopen("https://api.github.com/repos/ventolab/cellphonedb-data/git/trees/master?recursive=1")
    files_data = json.load(r)['tree']
    for rec in files_data:
        if rec['path'].startswith('data/sources/'):
            fname = rec['path'].split('/')[-1]
            url = 'https://raw.githubusercontent.com/ventolab/cellphonedb-data/master/{}'.format(rec['path'])
            print("Downloading: " + fname)
            r = urllib.request.urlopen(url)
            with open(os.path.join(sources_path, fname), 'wb') as f:
                f.write(r.read())

def create_db(user_dir_root, db_version, use_local_files):
    db_files_path = get_db_path(user_dir_root, db_version)
    pathlib.Path(db_files_path).mkdir(parents=True, exist_ok=True)
    # Get input files
    if use_local_files == False:
        download_input_files(db_files_path)
    dataDFs = getDFs(db_files_path, INPUT_FILE_NAMES, "csv")

    # Collect protein data
    protein_db_df = dataDFs['protein_input'][['protein_name','tags','tags_reason','tags_description']]
    num_proteins = protein_db_df.shape[0]
    multidata_id_list_so_far = list(range(num_proteins))
    protein_db_df.insert(0, 'id_protein', multidata_id_list_so_far, False)
    protein_db_df.insert(len(protein_db_df.columns), 'protein_multidata_id', protein_db_df['id_protein'].tolist(), True)
    # dbg(protein_db_df.info)

    # Collect gene data
    gene_db_df = dataDFs['gene_input'][['ensembl', 'gene_name', 'hgnc_symbol']]
    num_genes = gene_db_df.shape[0]
    gene_db_df.insert(0, 'id_gene', list(range(num_genes)), False)
    # Assign values from protein_db_df['protein_multidata_id'] into gene_db_df['protein_id']
    # via join between 'gene_name'_HUMAN and 'protein_name'
    gene_db_df.insert(len(gene_db_df.columns), 'protein_name', gene_db_df['gene_name'].apply(lambda x: x + '_HUMAN').tolist(), True)
    gene_db_df = pd.merge(gene_db_df, protein_db_df[['protein_name', 'protein_multidata_id']], on='protein_name')
    gene_db_df = gene_db_df.drop('protein_name', axis=1)
    gene_db_df.rename(columns = {'protein_multidata_id':'protein_id'}, inplace = True)
    # dbg(gene_db_df.info)

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
    multidata_db_df = multidata_db_df.append(complex_aux_df, ignore_index=True, verify_integrity=True)
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
    with open(os.path.join(db_files_path,'cellphonedb.zip'), 'wb') as f:
        f.write(zip_buffer.getvalue())
    print("Created CellphoneDB {} in {} successfully" \
          .format(db_version, os.path.join(db_files_path,'cellphonedb.zip')))

def convert_to_h5ad(user_dir_root):
    counts, raw_meta, meta, microenvs, degs = get_user_files(user_dir_root, \
                                                                        counts_fn=KEY2USER_TEST_FILE['counts'],
                                                                        meta_fn=KEY2USER_TEST_FILE['meta'],
                                                                        microenvs_fn=KEY2USER_TEST_FILE['microenvs'],
                                                                        degs_fn='test_degs.txt')
    obs = pd.DataFrame()
    obs['sample'] = counts.columns
    # dataframe for annotating the variables
    var = pd.DataFrame(index=counts.index)
    adata = anndata.AnnData(counts.T.values, obs=obs, var=var, dtype='float64')
    for key in [x for x in KEY2USER_TEST_FILE.keys() if x != 'counts']:
        with open(os.path.join(user_dir_root,'user_files',KEY2USER_TEST_FILE[key]),
                  'r') as content_file:
            content = content_file.read()
            adata.uns[key] = content
    outputPath = os.path.join(user_dir_root,"user_files","test.h5ad")
    adata.write(outputPath)

if __name__ == '__main__':
    arg = sys.argv[1]
    if arg == 'a':
        testrun_analyses(CPDB_ROOT, RELEASED_VERSION)
    elif arg == 'db':
        create_db(CPDB_ROOT, RELEASED_VERSION)
    elif arg == 'c':
        convert_to_h5ad(CPDB_ROOT)
    elif arg == 'g':
        download_source_files(CPDB_ROOT, RELEASED_VERSION)
        data_dir = get_db_path(CPDB_ROOT, RELEASED_VERSION)
        generated_path = os.path.join(data_dir, "generated")
        print("Generating gene_generated.csv file into {}".format(generated_path))
        pathlib.Path(generated_path).mkdir(parents=True, exist_ok=True)
        generate_input_files.generate_genes(data_dir,
                       user_gene=None,
                       fetch_uniprot=False,
                       fetch_ensembl=False,
                       result_path=generated_path,
                       project_name=None,
                       )
        print("Generating proteins_generated.csv file into {}".format(generated_path))
        generate_input_files.generate_proteins(data_dir,
                          user_protein=None,
                          fetch_uniprot=False,
                          result_path=generated_path,
                          log_file="log.txt",
                          project_name=None)
        print("Generating complex_generated.csv file into {}".format(generated_path))
        generate_input_files.generate_complex(data_dir,
                     user_complex=None,
                     result_path=generated_path,
                     log_file='log.txt',
                     project_name=None)
        print("Generating interactions_input.csv file into {}".format(generated_path))
        generate_input_files.generate_interactions(data_dir,
                              proteins=os.path.join(generated_path, 'protein_generated.csv'),
                              genes=os.path.join(generated_path, 'gene_generated.csv'),
                              complex=os.path.join(generated_path, 'complex_generated.csv'),
                              user_interactions=None,
                              user_interactions_only=False,
                              result_path=generated_path,
                              fetch_imex=False,
                              fetch_iuphar=False,
                              project_name=None,
                              release=False)
        print("Generating gene, protein and complex input files file into {}".format(generated_path))
        generate_input_files.filter_all(input_path=generated_path,
               result_path=generated_path,
               project_name=None)
    else:
        print("Arguments can be a (perform analysis) or db (create database)")

