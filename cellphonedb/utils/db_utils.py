from typing import Tuple
import pandas as pd
import numpy as np
import os
import io
import re
import time
import zipfile
import itertools
import pathlib
from cellphonedb.utils.file_utils import dbg
from cellphonedb.utils import file_utils, unique_id_generator
import urllib.request
import urllib.error
import urllib.parse
from zipfile import ZipFile
from cellphonedb.src.core.exceptions.DatabaseCreationException import DatabaseCreationException

MULTIDATA_TABLE_BOOLEAN_COLS = ['receptor', 'other', 'secreted_highlight',
                                'transmembrane', 'secreted', 'peripheral', 'integrin', 'is_complex']

PROTEIN_INFO_FIELDS_FOR_WEB = ['transmembrane', 'secreted', 'secreted_desc', 'receptor', 'integrin', 'other_desc']
COMPLEX_INFO_FIELDS_FOR_WEB = ['transmembrane', 'peripheral', 'secreted', 'secreted_desc', 'receptor', 'integrin',
                               'other_desc']
COMPLEX_CROSSREFERENCE_FIELDS_FOR_WEB = ['reactome_reaction', 'reactome_complex', 'complexPortal_complex',
                                         'rhea_reaction']

INPUT_FILE_NAMES = ['complex_input', 'gene_input', 'interaction_input', 'protein_input', 'transcription_factor_input']
# This is used to indicate CellPhoneDB released data (as opposed to user-added data
# when they create their own CellPhoneDB file)
CORE_CELLPHONEDB_DATA = "CellPhoneDBcore"


def get_protein_and_complex_data_for_web(cpdb_file_path) -> Tuple[dict, dict, dict, dict]:
    # Extract csv files from db_files_path/cellphonedb.zip into dbTableDFs
    dbTableDFs = extract_dataframes_from_db(cpdb_file_path)
    mtTable = dbTableDFs['multidata_table'].copy()
    cpxTable = dbTableDFs['complex_table']
    proteinTable = dbTableDFs['protein_table']

    for col in set(PROTEIN_INFO_FIELDS_FOR_WEB + COMPLEX_INFO_FIELDS_FOR_WEB):
        mtTable = mtTable.astype({col: 'str'})
        mtTable.loc[mtTable[col] == "True", col] = col.capitalize()
        mtTable.loc[mtTable[col] == "False", col] = np.nan
        if col in ['other_desc']:
            # Sanitize values for displaying to the user
            mtTable[col] = mtTable[col].str.replace("_", " ").str.capitalize()

    mtp = mtTable[~mtTable['is_complex']]
    aux = pd.merge(mtp, proteinTable, left_on='id_multidata', right_on='protein_multidata_id')
    proteinAcc2Name = dict(zip(aux['name'], aux['protein_name']))

    mtc = mtTable[mtTable['is_complex']]

    aux = dict(zip(mtp['name'], mtp[PROTEIN_INFO_FIELDS_FOR_WEB].values))
    protein2Info = {k: [x for x in aux[k] if str(x) != 'nan'] for k in aux}
    aux = dict(zip(mtc['name'], mtc[COMPLEX_INFO_FIELDS_FOR_WEB].values))
    complex2Info = {k: [x for x in aux[k] if str(x) != 'nan'] for k in aux}

    aux = pd.merge(mtc, cpxTable, left_on='id_multidata', right_on='complex_multidata_id')
    complex_cols = list(aux.columns.values.tolist())
    resource2Complex2Acc = {}
    for col in COMPLEX_CROSSREFERENCE_FIELDS_FOR_WEB:
        if col in complex_cols:
            # The above test is in case the user created their own CellphoneDB database and had chosen to remove
            # COMPLEX_CROSSREFERENCE_FIELDS_FOR_WEB fields from complex_input.csv
            aux1 = aux.loc[pd.notna(aux[col])]
            resource2Complex2Acc[col.replace("_", " ").capitalize()] = dict(zip(aux1['name'], aux1[col]))

    return protein2Info, complex2Info, resource2Complex2Acc, proteinAcc2Name


def get_interactions_genes_complex(cpdb_file_path) -> \
        Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame, dict, dict]:
    """
    Returns a tuple of four DataFrames containing data from <cpdb_dir>/cellphonedb.zip.

    Parameters
    ----------
    cpdb_file_path: str
        CellphoneDB database file path

    Returns
    -------
    Tuple
        - interactions: pd.DataFrame
        - genes: pd.DataFrame
        - complex_composition: pd.DataFrame
        - complex_expanded: pd.DataFrame
        - gene_synonym2gene_name: dict
        - receptor2tfs: dict
    """
    # Extract csv files from db_files_path/cellphonedb.zip into dbTableDFs
    dbTableDFs = extract_dataframes_from_db(cpdb_file_path)
    # Convert dbTableDFs into interactions, genes, complex_composition, complex_expanded data frames
    # and gene_synonym2gene_name dict
    gene_synonym2gene_name = {}
    # Cater for DB version-dependent input files
    if 'gene_synonym_to_gene_name' in dbTableDFs:
        gs2gn = dbTableDFs['gene_synonym_to_gene_name']
        gene_synonym2gene_name = dict(zip(gs2gn['Gene Synonym'], gs2gn['Gene Name']))

    mtTable = dbTableDFs['multidata_table']
    dbg(mtTable.dtypes)
    # Convert all MULTIDATA_TABLE_BOOLEAN_COLS from Integer (0/1) to Boolean
    for col in MULTIDATA_TABLE_BOOLEAN_COLS:
        mtTable[col] = mtTable[col].astype(bool)
    # Read genes 'table' - c.f. old CellphoneDB: GeneRepository.get_all_expanded()
    # Drop 'protein_name' column from dbTableDFs['gene_table'] as dbTableDFs['protein_table'] already has it
    if 'protein_name' in dbTableDFs['gene_table'].columns:
        dbTableDFs['gene_table'] = dbTableDFs['gene_table'].drop('protein_name', axis=1)
    # First filter out entries where gene_name = X with no hgnc_symbol, e.g. in the case of IGF2, filter out ENSG00000284779
    dbTableDFs['gene_table'] = dbTableDFs['gene_table'][~dbTableDFs['gene_table']['hgnc_symbol'].isnull()]
    # Now merge gene_table with protein_table and multidata_table
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
        multidata_simple = pd.merge(dbTableDFs['protein_table'], mtTable,
                                    left_on='protein_multidata_id', right_on='id_multidata')
    multidata_complex = pd.merge(mtTable, dbTableDFs['complex_table'],
                                 left_on='id_multidata', right_on='complex_multidata_id')
    if multidata_complex.empty:
        multidata_expanded = multidata_simple
    else:
        multidata_expanded = pd.concat([multidata_simple, multidata_complex], ignore_index=True, sort=True)
    # C.f. old CellphoneDB: InteractionRepository.get_all_expanded()
    multidata_expanded = multidata_expanded.astype({'id_multidata': 'int64'})
    dbg(multidata_expanded.columns)
    interactions = pd.merge(dbTableDFs['interaction_table'], multidata_expanded, how='left',
                            left_on=['multidata_1_id'], right_on=['id_multidata'])
    interactions = pd.merge(interactions, multidata_expanded, how='left',
                            left_on=['multidata_2_id'], right_on=['id_multidata'], suffixes=suffixes)
    dbg("interactions columns: ", interactions.columns)
    # Generate complex_expanded - c.f. ComplexRepository.get_all_expanded()
    complex_expanded = pd.merge(dbTableDFs['complex_table'], mtTable,
                                left_on='complex_multidata_id', right_on='id_multidata')
    # dbg(complexes_expanded.columns)
    # index interactions and complex data frames
    # C.f. old CellphoneDB: method_launcher.get_interactions_genes_complex()
    interactions.set_index('id_interaction', drop=True, inplace=True)
    complex_composition.set_index('id_complex_composition', inplace=True, drop=True)

    # Cater for DB version-dependent input files
    receptor2tfs = {}
    if 'receptor_to_transcription_factor' in dbTableDFs:
        receptor_to_tf_df = dbTableDFs['receptor_to_transcription_factor'][['Receptor', 'TF']]
        for receptor, tf in receptor_to_tf_df.values:
            receptor2tfs.update({receptor: receptor2tfs.get(receptor, []) + [tf]})

    return interactions, genes, complex_composition, complex_expanded, gene_synonym2gene_name, receptor2tfs


def extract_dataframes_from_db(cpdb_file_path):
    dfs = {}
    start = time.time()
    for tuple in unzip(cpdb_file_path):
        file_name = tuple[0]
        file_handle = tuple[1]
        dbg("Retrieving from zip file: " + file_name)
        dfs[file_name.replace('.csv', '')] = pd.read_csv(file_handle)
    duration = time.time() - start
    dbg("Loaded DB into memory in " + str(round(duration, 2)) + "s")
    return dfs


def unzip(zip_file_path):
    """
    Download a ZIP file and extract its contents in memory
    yields (filename, file-like object) pairs
    """
    with open(zip_file_path, 'br') as file:
        # Note - eval unescapes the double-quotes
        content = file.read()
        dbg(type(content))
        with zipfile.ZipFile(io.BytesIO(content)) as thezip:
            for zipinfo in thezip.infolist():
                with thezip.open(zipinfo) as thefile:
                    yield zipinfo.filename, thefile


def get_db_path(user_dir_root, db_version):
    """
    Retrieves the path to the local database file corresponding to db_version

    Parameters
    ----------
    user_dir_root: str
        The directory in which user stores CellphoneDB files
    db_version: str
        CellphoneDB version (and the name of the subdirectory containing the
        curated input files from https://github.com/ventolab/cellphonedb-data)

    Returns
    -------
    str
        The path to the local database file corresponding to db_version
    """
    return os.path.join(user_dir_root, "releases", db_version)


# Cater for DB version-dependent column names
def get_column_names_for_db_version(complex_db_df, interactions_df, protein_df) -> tuple:
    protein_column_names = ['uniprot_1', 'uniprot_2', 'uniprot_3', 'uniprot_4']
    interaction_column_names1 = []
    interaction_column_names2 = []
    complex_columns = []
    version = []

    if 'directionality' in interactions_df.columns:
        interaction_column_names1 = ['directionality', 'classification']
        interaction_column_names2 = ['is_ppi', 'curator']
    if 'uniprot_5' in complex_db_df.columns:
        protein_column_names += ['uniprot_5']
        complex_columns = COMPLEX_CROSSREFERENCE_FIELDS_FOR_WEB
    if 'version' in protein_df.columns:
        version = ['version']
    return (protein_column_names, interaction_column_names1, interaction_column_names2, version, complex_columns)


def collect_protein_data(data_dfs: dict) -> pd.DataFrame:
    # Collect protein data
    protein_db_df = data_dfs['protein_input'][['protein_name', 'tags', 'tags_reason', 'tags_description', 'uniprot']]
    num_proteins = protein_db_df.shape[0]
    multidata_id_list_so_far = list(range(num_proteins))
    protein_db_df.insert(0, 'id_protein', multidata_id_list_so_far, False)
    protein_db_df.insert(len(protein_db_df.columns), 'protein_multidata_id', protein_db_df['id_protein'].tolist(), True)
    # dbg(protein_db_df.info)
    return protein_db_df, multidata_id_list_so_far


def collect_gene_data(data_dfs: dict, protein_db_df: pd.DataFrame) -> pd.DataFrame:
    gene_db_df = data_dfs['gene_input'][['ensembl', 'gene_name', 'hgnc_symbol', 'uniprot']]
    num_genes = gene_db_df.shape[0]
    gene_db_df.insert(0, 'id_gene', list(range(num_genes)), False)
    # Assign values from protein_db_df['protein_multidata_id'] into gene_db_df['protein_id']
    # via join between 'uniprot' and 'protein_name'
    gene_db_df = pd.merge(gene_db_df, protein_db_df[['protein_name', 'protein_multidata_id', 'uniprot']], on='uniprot')
    gene_db_df = gene_db_df.drop('uniprot', axis=1)
    protein_db_df = protein_db_df.drop('uniprot', axis=1)
    gene_db_df.rename(columns={'protein_multidata_id': 'protein_id'}, inplace=True)
    # print(gene_db_df.info)
    return gene_db_df


def collect_receptor_to_tf_mapping(data_dfs: dict) -> pd.DataFrame:
    receptor_to_tf_df = None
    # Cater for DB version-dependent input files
    if data_dfs['transcription_factor_input'] is not None:
        receptor_to_tf_df = data_dfs['transcription_factor_input'][['receptor_id', 'TF_symbol']].copy()
        receptor_to_tf_df.columns = ['Receptor', 'TF']
        # Strip any leading or trailing spaces
        receptor_to_tf_df.replace(r'\s*(.*?)\s*', r'\1', regex=True, inplace=True)
    return receptor_to_tf_df


def collect_gene_synonym_to_gene_name_mapping(data_dfs: dict, gene_db_df: pd.DataFrame) -> pd.DataFrame:
    gene_synonym_to_gene_name_db_df = None
    # Cater for DB version-dependent input files
    if data_dfs['gene_synonyms_input'] is not None:
        gene_synonym_to_gene_name = {}
        for gene_names in data_dfs['gene_synonyms_input']\
                .filter(regex=("Gene Names.*")).dropna().agg(' '.join, axis=1).tolist():
            gene_names_arr = re.split(';\\s*|\\s+', gene_names)
            for gene_name in gene_db_df[gene_db_df['gene_name'].isin(gene_names_arr)]['gene_name'].tolist():
                for gene_synonym in gene_names_arr:
                    if gene_synonym != gene_name:
                        gene_synonym_to_gene_name[gene_synonym] = gene_name
        gene_synonym_to_gene_name_db_df = pd.DataFrame(gene_synonym_to_gene_name.items(),
                                                       columns=['Gene Synonym', 'Gene Name'])
    return gene_synonym_to_gene_name_db_df


def create_db(target_dir) -> None:
    """
    Creates CellphoneDB databases file (cellphonedb.zip) in <target_dir> directory.
    The assumption is that <target_dir> contains the four *_input.csv files required to
    create the database file.
    This simple zip file contains a number of CSV files that CellphoneDB package reads into memory
    and uses to drive its analysis of user data.

    Parameters
    ----------
    target_dir: str
        The directory in which to place the database file (and which contains the four *_input.csv files required to
    create the database file

    Returns
    -------

    """
    gene_input = os.path.join(target_dir, "gene_input.csv")
    protein_input = os.path.join(target_dir, "protein_input.csv")
    complex_input = os.path.join(target_dir, "complex_input.csv")
    interaction_input = os.path.join(target_dir, "interaction_input.csv")
    transcription_factor_input = os.path.join(target_dir, "transcription_factor_input.csv")
    gene_synonyms_input = os.path.join(target_dir, "sources/uniprot_synonyms.tsv")

    pathlib.Path(target_dir).mkdir(parents=True, exist_ok=True)
    data_dfs = get_dfs(gene_input=gene_input, protein_input=protein_input, complex_input=complex_input,
                       interaction_input=interaction_input, transcription_factor_input=transcription_factor_input,
                       gene_synonyms_input=gene_synonyms_input)

    (protein_column_names, interaction_column_names1, interaction_column_names2, version, complex_columns) = \
        get_column_names_for_db_version(
            data_dfs['complex_input'], data_dfs['interaction_input'], data_dfs['protein_input'])

    # Perform sanity tests on *_input files and report any issues to the user as warnings
    run_sanity_tests(data_dfs, protein_column_names, version)

    # Collect protein data
    protein_db_df, multidata_id_list_so_far = collect_protein_data(data_dfs)

    # Collect gene data
    gene_db_df = collect_gene_data(data_dfs, protein_db_df)

    # Collect mapping: (receptor) gene name -> TF gene name (in transcription_factor_input.tsv)
    receptor_to_tf_df = collect_receptor_to_tf_mapping(data_dfs)

    # Collect mapping: gene synonym (not in gene_input.csv) -> gene name (in gene_input.csv)
    gene_synonym_to_gene_name_db_df = collect_gene_synonym_to_gene_name_mapping(data_dfs, gene_db_df)

    # Collect complex data
    cols = protein_column_names + ['pdb_structure', 'pdb_id', 'stoichiometry', 'comments_complex'] + complex_columns
    complex_db_df = data_dfs['complex_input'][cols]

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
        data_dfs['protein_input'][['uniprot', 'receptor', 'receptor_desc', 'other', 'other_desc', 'secreted_highlight',
                                  'secreted_desc', 'transmembrane', 'secreted', 'peripheral', 'integrin']].copy()
    multidata_db_df.rename(columns={'uniprot': 'name'}, inplace=True)
    multidata_ids = pd.merge(
        data_dfs['protein_input'][['protein_name']],
        protein_db_df[['protein_name', 'protein_multidata_id']], on='protein_name')['protein_multidata_id'].tolist()
    multidata_db_df.insert(0, 'id_multidata', multidata_ids, False)
    multidata_db_df.insert(len(multidata_db_df.columns), 'is_complex',
                           list(itertools.repeat(False, multidata_db_df.shape[0])), True)
    dbg(multidata_db_df.shape, multidata_db_df.index, multidata_db_df.columns)
    # Insert complexes into multidata
    cols = ['complex_name', 'receptor', 'receptor_desc', 'other', 'other_desc', 'secreted_highlight', 'secreted_desc',
            'transmembrane', 'secreted', 'peripheral', 'integrin']
    complex_aux_df = data_dfs['complex_input'][cols].copy()
    complex_aux_df.rename(columns={'complex_name': 'name'}, inplace=True)
    complex_aux_df.insert(0, 'id_multidata', complex_multidata_ids, False)
    complex_aux_df.insert(len(complex_aux_df.columns), 'is_complex',
                          list(itertools.repeat(True, complex_aux_df.shape[0])), True)
    dbg(complex_aux_df.shape, complex_aux_df.index, complex_aux_df.columns)
    # Append complex_aux_df to multidata_db_df
    multidata_db_df = pd.concat([multidata_db_df, complex_aux_df], ignore_index=True, verify_integrity=True)
    dbg(multidata_db_df.shape, multidata_db_df.index, multidata_db_df.columns)

    # First collect total_protein counts for each complex in complex_db_df
    total_protein_cnt_list = np.apply_along_axis(
        lambda s: sum(isinstance(x, str) for x in s), 1, complex_db_df[protein_column_names].values).tolist()
    complex_db_df.insert(len(complex_db_df.columns), 'total_protein', total_protein_cnt_list, True)
    dbg(complex_db_df.info)
    # Next collect all complex_composition data into cc_list
    cc_list = []
    pos = len(protein_column_names)
    for r in complex_db_df[protein_column_names + ['complex_multidata_id', 'total_protein']].values.tolist():
        for acc in filter(lambda x: isinstance(x, str), r):
            protein_multidata_id = \
                multidata_db_df.loc[(~multidata_db_df['is_complex']) &
                                    (multidata_db_df['name'] == acc), ['id_multidata']].iat[0, 0]
            complex_multidata_id = r[pos]
            total_protein = r[pos+1]
            cc_list.append([complex_multidata_id, protein_multidata_id, total_protein])

    complex_composition_df = pd.DataFrame(cc_list, columns=['complex_multidata_id', 'protein_multidata_id', 'total_protein'])
    complex_composition_df.insert(0, 'id_complex_composition', list(range(len(cc_list))), False)
    dbg(complex_composition_df.shape, complex_composition_df.index, complex_composition_df.columns,
        complex_composition_df.info)
    # Next drop the auxiliary columns from complex_db_df: protein_column_names and 'total_protein'
    for col in protein_column_names + ['total_protein']:
        complex_db_df = complex_db_df.drop(col, axis=1)

    # Collect interaction data
    interactions_aux_df = pd.merge(data_dfs['interaction_input'], multidata_db_df,
                                   left_on=['partner_a'], right_on=['name'])
    interactions_aux_df = pd.merge(interactions_aux_df, multidata_db_df,
                                   left_on=['partner_b'], right_on=['name'], suffixes=['_x', '_y'])
    dbg(interactions_aux_df.shape)
    # Remove interactions non-CPDB interactors
    interactions_aux_df = interactions_aux_df[
        interactions_aux_df.apply(
            lambda interaction: file_utils.is_cellphonedb_interactor(interaction, ('_x', '_y')), axis=1)]
    interactions_aux_df['id_cp_interaction'] = interactions_aux_df.apply(
        lambda interaction: unique_id_generator.interaction(interaction, ('_x', '_y')), axis=1)
    dbg(interactions_aux_df.info)
    dbg(interactions_aux_df.columns)
    interactions_df = interactions_aux_df[['id_cp_interaction', 'id_multidata_x', 'id_multidata_y',
                                          'source', 'annotation_strategy'] +
                                          interaction_column_names2 + interaction_column_names1].copy()
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
        if gene_synonym_to_gene_name_db_df is not None:
            # Cater for DB version-dependent input files
            zip_file.writestr('gene_synonym_to_gene_name.csv',
                              gene_synonym_to_gene_name_db_df.to_csv(index=False, sep=',').encode('utf-8'))
        if receptor_to_tf_df is not None:
            # Cater for DB version-dependent input files
            zip_file.writestr('receptor_to_transcription_factor.csv',
                              receptor_to_tf_df.to_csv(index=False, sep=',').encode('utf-8'))

    file_suffix = file_utils.get_timestamp_suffix()
    file_path = os.path.join(target_dir, 'cellphonedb_{}.zip'.format(file_suffix))
    with open(file_path, 'wb') as f:
        f.write(zip_buffer.getvalue())
    print("Created {} successfully".format(file_path))


def download_database(target_dir, cpdb_version):
    download_released_files(target_dir, cpdb_version, "cellphonedb.zip|_input|sources\\/uniprot_synonyms")


def download_released_files(target_dir, cpdb_version, regex):
    r = urllib.request.urlopen('https://github.com/ventolab/cellphonedb-data/archive/refs/tags/{}.zip'.format(cpdb_version))
    zipContent = ZipFile(io.BytesIO(r.read()))
    for fpath in zipContent.namelist():
        if re.search(regex, fpath):
            fname = fpath.split("/")[-1]
            if fname:
                if re.search("sources", fpath):
                    target_dir = os.path.join(target_dir, "sources")
                pathlib.Path(target_dir).mkdir(parents=True, exist_ok=True)
                with open(os.path.join(target_dir, fname), 'wb') as f:
                    f.write(zipContent.read(fpath))
                    print("Downloaded {} into {}".format(fname, target_dir))


def get_dfs(gene_input=None, protein_input=None, complex_input=None, interaction_input=None,
            transcription_factor_input=None, gene_synonyms_input=None):
    dfs = {}
    dfs['gene_input'] = file_utils.read_data_table_from_file(gene_input)
    dfs['protein_input'] = file_utils.read_data_table_from_file(protein_input)
    dfs['complex_input'] = file_utils.read_data_table_from_file(complex_input)
    dfs['interaction_input'] = file_utils.read_data_table_from_file(interaction_input)
    dfs['transcription_factor_input'] = file_utils.read_data_table_from_file(transcription_factor_input, optional=True)
    dfs['gene_synonyms_input'] = file_utils.read_data_table_from_file(gene_synonyms_input, optional=True)
    return dfs


def sanity_test_uniprot_accessions_map_to_multiple_gene_names(gene_db_df: pd.DataFrame):
    gene_names_uniprot_df = gene_db_df[['gene_name', 'uniprot']].copy()
    gene_names_uniprot_df.drop_duplicates(inplace=True)
    dups = gene_names_uniprot_df[gene_names_uniprot_df['uniprot'].duplicated()]
    if not dups.empty:
        # data_errors_found = True
        print("WARNING: The following UniProt ids map to multiple gene names (it is expected that " +
              "they should map to only one):")
        print(", ".join(dups['uniprot'].tolist()))


def sanity_test_report_complex_duplicates(complex_db_df: pd.DataFrame):
    test_complex_db_df = complex_db_df.set_index('complex_name')
    if test_complex_db_df.index.has_duplicates:
        print("WARNING: complex_input.csv has the following duplicates:")
        print("\n".join(complex_db_df[test_complex_db_df.index.duplicated(keep='first')]['complex_name'].tolist()) + "\n")


def sanity_test_report_complexes_with_same_participants(
        complex_db_df: pd.DataFrame, version: str, protein_column_names: list):
    participants_set_to_complex_names = {}
    participants_set_to_data_sources = {}
    cols = ['complex_name'] + version + protein_column_names
    start_idx = cols.index('uniprot_1')
    for row in complex_db_df[cols].itertuples(index=False):
        participants_set = frozenset([i for i in row[start_idx:] if str(i) != 'nan'])
        complex_name = row[0]
        data_source = row[1]
        m = re.search(r"^{}".format(CORE_CELLPHONEDB_DATA), data_source)
        if m:
            # Store in data_source only CORE_CELLPHONEDB_DATA (i.e. exclude any version information)
            # Here we just need to know that this complex was added by the CellphoneDB team
            data_source = m[0]
        if participants_set not in participants_set_to_complex_names:
            participants_set_to_complex_names[participants_set] = [complex_name]
            participants_set_to_data_sources[participants_set] = set([data_source])
        else:
            participants_set_to_complex_names[participants_set].append(complex_name)
            participants_set_to_data_sources[participants_set].add(data_source)
    complex_dups = ""
    for participants_set in participants_set_to_complex_names:
        data_sources = list(participants_set_to_data_sources[participants_set])
        if (len(data_sources) > 1 or not data_sources or data_sources[0] != CORE_CELLPHONEDB_DATA):
            complex_names = participants_set_to_complex_names[participants_set]
            if len(complex_names) > 1:
                complex_dups += ", ".join(complex_names) + " : " + ", ".join(participants_set) + "\n"
    if len(complex_dups) > 0:
        # data_errors_found = True
        print("WARNING: The following multiple complexes (left) appear to have the same composition (right):")
        print(complex_dups)


def sanity_test_report_interactions_with_same_participants(interaction_db_df: pd.DataFrame):
    partner_sets = [set([i for i in row]) for row in
                    interaction_db_df[['partner_a', 'partner_b']].itertuples(index=False)]
    # Find duplicate sets of partners
    seen = set()
    duplicate_partner_sets = [x for x in partner_sets if x in seen or seen.add(frozenset(x))]
    if duplicate_partner_sets:
        # data_errors_found = True
        print("WARNING: The following sets of interaction partners appear in multiple rows of interaction_input.csv file:")
        for dup in set([frozenset(x) for x in duplicate_partner_sets]):
            print(','.join(dup))
        print()


def sanity_test_report_uniprot_accession_duplicates(protein_db_df: pd.DataFrame):
    test_protein_db_df = protein_db_df.set_index('uniprot')
    if test_protein_db_df.index.has_duplicates:
        print("WARNING: protein_input.csv has the following UniProt accession duplicates:")
        print("\n".join(protein_db_df[test_protein_db_df.index.duplicated(keep='first')]['uniprot'].tolist()) + "\n")


def sanity_test_report_orphan_complexes(complex_db_df: pd.DataFrame, interaction_db_df: pd.DataFrame) -> set:
    all_complexes_set = set(complex_db_df['complex_name'].tolist())
    interaction_participants_set = set(interaction_db_df['partner_a'].tolist() + interaction_db_df['partner_b'].tolist())
    orphan_complexes = all_complexes_set - interaction_participants_set
    if orphan_complexes:
        print("WARNING: The following complexes are not found in interaction_input.txt:")
        print("\n".join(orphan_complexes))
    print()
    return orphan_complexes


def sanity_test_report_orphan_proteins(
        protein_db_df: pd.DataFrame,
        complex_db_df: pd.DataFrame,
        interaction_db_df: pd.DataFrame,
        protein_column_names: list,
        orphan_complexes: set):
    all_proteins_set = set(protein_db_df['uniprot'].tolist())
    interaction_participants_set = set(
        interaction_db_df['partner_a'].tolist() + interaction_db_df['partner_b'].tolist())
    proteins_in_complexes_participating_in_interactions = []
    for colName in protein_column_names:
        proteins_in_complexes_participating_in_interactions += \
            complex_db_df[~complex_db_df['complex_name'].isin(orphan_complexes)][colName].tolist()
    orphan_proteins = all_proteins_set - interaction_participants_set - \
        set(proteins_in_complexes_participating_in_interactions)
    if orphan_proteins:
        print("WARNING: The following proteins are not found in interaction_input.txt (either directly " +
              "or via complexes they are part of):")
        print("\n".join(orphan_proteins))


def sanity_test_report_unknown_interactors(
        protein_db_df: pd.DataFrame,
        complex_db_df: pd.DataFrame,
        interaction_db_df: pd.DataFrame
):
    unknown_interactors = set()
    for col in ['partner_a', 'partner_b']:
        aux_df = pd.merge(interaction_db_df, protein_db_df, left_on=col, right_on='uniprot', how='outer')
        unknown_interactor_proteins = set(aux_df[pd.isnull(aux_df['uniprot'])][col].tolist())
        aux_df = pd.merge(interaction_db_df, complex_db_df, left_on=col, right_on='complex_name', how='outer')
        unknown_interactor_complexes = set(aux_df[pd.isnull(aux_df['complex_name'])][col].tolist())
        unknown_interactors = unknown_interactors.union(
            unknown_interactor_proteins.intersection(unknown_interactor_complexes))
    if unknown_interactors:
        print("WARNING: The following interactors in interaction_input.txt could not be found in either " +
              "protein_input.csv or complex_indput.csv:")
        print("\n".join(sorted(unknown_interactors)) + "\n")


def sanity_test_report_unknown_proteins(
        protein_db_df: pd.DataFrame,
        complex_db_df: pd.DataFrame,
        protein_column_names: list):
    unknown_proteins = set()
    for col in protein_column_names:
        aux_df = pd.merge(complex_db_df, protein_db_df, left_on=col, right_on='uniprot', how='outer')
        unknown_complex_proteins = set(aux_df[pd.isnull(aux_df['uniprot']) & ~pd.isnull(aux_df[col])][col].tolist())
        unknown_proteins = unknown_proteins.union(unknown_complex_proteins)
    if unknown_proteins:
        print("WARNING: The following proteins in complex_input.txt could not be found in protein_input.csv:")
        print("\n".join(sorted(unknown_proteins)) + "\n")


def sanity_test_report_proteins_not_in_genes_file(
        protein_db_df: pd.DataFrame,
        gene_db_df: pd.DataFrame):
    proteins = set(protein_db_df['uniprot'].tolist())
    unknown_proteins = proteins.difference(set(gene_db_df['uniprot'].tolist()))
    if unknown_proteins:
        print("WARNING: The following proteins in protein_input.txt could not be found in gene_input.csv:")
        print("\n".join(sorted(unknown_proteins)) + "\n")
        print()


def sanity_test_report_tfs_not_in_gene_or_complex_files(
        gene_db_df: pd.DataFrame,
        complex_db_df: pd.DataFrame,
        tf_input_df: pd.DataFrame):
    if tf_input_df is not None:
        # Cater for DB version-dependent input files
        for (bioentity, df) in {"gene": gene_db_df, "complex": complex_db_df}.items():
            if bioentity == "gene":
                complex_filter = ~tf_input_df['receptor_id'].str.match('.*_.*')
            else:
                complex_filter = tf_input_df['receptor_id'].str.match('.*_.*')
            bioentities_in_tf_input = set([i.strip() for i in tf_input_df[complex_filter]['receptor_id'].values.tolist()])
            bioentities_in_input = set(df['{}_name'.format(bioentity)].values.tolist())
            # Below: bioentities in bioentities_in_tf_input but not in bioentities_in_input
            bioentities_not_in_input = bioentities_in_tf_input.difference(bioentities_in_input)
            if bioentities_not_in_input:
                print("WARNING: The following receptors in transcription_factor_input could not be found in " +
                      "{}_input.csv:".format(bioentity))
                print("\n".join(set(bioentities_not_in_input)))
                print()


def run_sanity_tests(data_dfs, protein_column_names, version):
    data_errors_found = False
    protein_db_df = data_dfs['protein_input']
    complex_db_df = data_dfs['complex_input']
    gene_db_df = data_dfs['gene_input']
    interaction_db_df = data_dfs['interaction_input']
    tf_input_df = None
    # Cater for DB version-dependent input files
    if 'transcription_factor_input' in data_dfs:
        tf_input_df = data_dfs['transcription_factor_input']

    # 1. Report any uniprot accessions that map to multiple gene_names
    sanity_test_uniprot_accessions_map_to_multiple_gene_names(gene_db_df)

    # 2. Warn about complex name duplicates in complex_db_df
    sanity_test_report_complex_duplicates(complex_db_df)

    # 3. Report complexes with (possibly) different names, but with the same uniprot
    # accession participants (though not necessarily in the same order - hence the use of set below)
    # NB. Use set below as we don't care about the order of participants when looking for duplicates
    # NB. Report duplicate complexes _only if_ at least one duplicate's complex_db_df['version']
    # does not start with CORE_CELLPHONEDB_DATA)
    sanity_test_report_complexes_with_same_participants(complex_db_df, version, protein_column_names)

    # 4. Report interactions with (possibly) a different name, but with the same participants
    # (though not necessarily in the same order - hence the use of set below)
    sanity_test_report_interactions_with_same_participants(interaction_db_df)

    # 5. Warn about uniprot accession duplicates in protein_db_df
    sanity_test_report_uniprot_accession_duplicates(protein_db_df)

    # 6. Warn the user if some complexes don't participate in any interactions
    orphan_complexes = sanity_test_report_orphan_complexes(complex_db_df, interaction_db_df)

    # 7. Warn the user if some proteins don't participate in any interactions directly,
    # or are part of some complex in orphan_complexes
    sanity_test_report_orphan_proteins(protein_db_df, complex_db_df, interaction_db_df, protein_column_names, orphan_complexes)

    # 8. Warn the user if some interactions contain interactors that are neither
    # in complex_input.csv or protein_input.csv
    sanity_test_report_unknown_interactors(protein_db_df, complex_db_df, interaction_db_df)

    # 9. Warn if some complexes contain proteins not in protein_input.csv
    sanity_test_report_unknown_proteins(protein_db_df, complex_db_df, protein_column_names)

    # 10. Warn if some proteins in protein_input.csv are not in gene_input.csv
    sanity_test_report_proteins_not_in_genes_file(protein_db_df, gene_db_df)

    # 11. Warn if some receptor ids in tf_input_df are in neither gene_input.csv or complex_input.csv
    sanity_test_report_tfs_not_in_gene_or_complex_files(gene_db_df, complex_db_df, tf_input_df)

    if data_errors_found:
        raise DatabaseCreationException
