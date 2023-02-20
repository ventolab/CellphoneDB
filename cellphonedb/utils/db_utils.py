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
import urllib.request, urllib.error, urllib.parse
from zipfile import ZipFile
from cellphonedb.src.core.exceptions.DatabaseCreationException import DatabaseCreationException

MULTIDATA_TABLE_BOOLEAN_COLS = ['receptor','other','secreted_highlight',\
                                'transmembrane','secreted','peripheral','integrin','is_complex']
INPUT_FILE_NAMES = ['complex_input','gene_input','interaction_input','protein_input']
PROTEIN_COLUMN_NAMES = ['uniprot_1','uniprot_2','uniprot_3','uniprot_4','uniprot_5']
# This is used to indicate CellPhoneDB released data (as opposed to user-added data when they create their own CellPhoneDB file)
CORE_CELLPHONEDB_DATA = "CellPhoneDBcore"

def get_interactions_genes_complex(cpdb_file_path) -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame, dict]:
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
    """
    # Extract csv files from db_files_path/cellphonedb.zip into dbTableDFs
    dbTableDFs = extract_dataframes_from_db(cpdb_file_path)
    # Convert dbTableDFs into interactions, genes, complex_composition, complex_expanded data frames
    # and gene_synonym2gene_name dict
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
        multidata_simple = pd.merge(dbTableDFs['protein_table'], mtTable, left_on='protein_multidata_id', right_on='id_multidata')
    multidata_complex = pd.merge(mtTable, dbTableDFs['complex_table'], left_on='id_multidata', right_on='complex_multidata_id')
    if multidata_complex.empty:
        multidata_expanded = multidata_simple
    else:
        multidata_expanded = pd.concat([multidata_simple, multidata_complex], ignore_index=True, sort=True)
    # C.f. old CellphoneDB: InteractionRepository.get_all_expanded()
    multidata_expanded = multidata_expanded.astype({'id_multidata': 'int64'})
    dbg(multidata_expanded.columns)
    interactions = pd.merge(dbTableDFs['interaction_table'], multidata_expanded, how='left', left_on=['multidata_1_id'], right_on=['id_multidata'])
    interactions = pd.merge(interactions, multidata_expanded, how='left', left_on=['multidata_2_id'], right_on=['id_multidata'], suffixes=suffixes)
    dbg("interactions columns: ", interactions.columns)
    # Generate complex_expanded - c.f. ComplexRepository.get_all_expanded()
    complex_expanded = pd.merge(dbTableDFs['complex_table'], mtTable, left_on='complex_multidata_id', right_on='id_multidata')
    # dbg(complexes_expanded.columns)
    # index interactions and complex data frames
    # C.f. old CellphoneDB: method_launcher.get_interactions_genes_complex()
    interactions.set_index('id_interaction', drop=True, inplace=True)
    complex_composition.set_index('id_complex_composition', inplace=True, drop=True)

    return interactions, genes, complex_composition, complex_expanded, gene_synonym2gene_name

def extract_dataframes_from_db(cpdb_file_path):
    dfs = {}
    start = time.time()
    for tuple in unzip(cpdb_file_path):
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
    gene_synonyms_input = os.path.join(target_dir, "sources/uniprot_synonyms.tsv")

    pathlib.Path(target_dir).mkdir(parents=True, exist_ok=True)
    dataDFs = getDFs(gene_input=gene_input, protein_input=protein_input, complex_input=complex_input,
                     interaction_input=interaction_input, gene_synonyms_input=gene_synonyms_input)

    # Perform sanity tests on *_input files and report any issues to the user as warnings
    run_sanity_tests(dataDFs)

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

    # Collect mapping: gene synonym (not in gene_input.csv) -> gene name (in gene_input.csv)
    gene_synonym_to_gene_name = {}
    for gene_names in dataDFs['gene_synonyms_input'].filter(regex=("Gene Names.*")).dropna().agg(' '.join, axis=1).tolist():
        gene_names_arr = re.split(';\s*|\s+', gene_names)
        for gene_name in gene_db_df[gene_db_df['gene_name'].isin(gene_names_arr)]['gene_name'].tolist():
            for gene_synonym in gene_names_arr:
                if gene_synonym != gene_name:
                    gene_synonym_to_gene_name[gene_synonym] = gene_name
    gene_synonym_to_gene_name_db_df = pd.DataFrame(gene_synonym_to_gene_name.items(), columns=['Gene Synonym', 'Gene Name'])

    # Collect complex data
    complex_db_df = dataDFs['complex_input'] \
        [PROTEIN_COLUMN_NAMES + ['pdb_structure','pdb_id','stoichiometry','comments_complex']]
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
    total_protein_cnt_list = np.apply_along_axis(lambda s: sum(type(x) == str for x in s), 1, complex_db_df[PROTEIN_COLUMN_NAMES].values).tolist()
    complex_db_df.insert(len(complex_db_df.columns), 'total_protein', total_protein_cnt_list, True)
    dbg(complex_db_df.info)
    # Next collect all complex_composition data into cc_list
    cc_list = []
    pos = len(PROTEIN_COLUMN_NAMES)
    for r in complex_db_df[PROTEIN_COLUMN_NAMES + ['complex_multidata_id','total_protein']].values.tolist():
        for acc in filter(lambda x: type(x) == str, r):
            protein_multidata_id = \
                multidata_db_df.loc[(multidata_db_df['is_complex'] == False) & (multidata_db_df['name'] == acc), ['id_multidata']] \
                    .iat[0,0]
            complex_multidata_id = r[pos]
            total_protein = r[pos+1]
            cc_list.append([complex_multidata_id, protein_multidata_id, total_protein])

    complex_composition_df = pd.DataFrame(cc_list, columns=['complex_multidata_id', 'protein_multidata_id', 'total_protein'])
    complex_composition_df.insert(0, 'id_complex_composition', list(range(len(cc_list))), False)
    dbg(complex_composition_df.shape, complex_composition_df.index, complex_composition_df.columns, complex_composition_df.info)
    # Next drop the auxiliary columns from complex_db_df: PROTEIN_COLUMN_NAMES and 'total_protein'
    for col in PROTEIN_COLUMN_NAMES + ['total_protein']:
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
            lambda interaction: file_utils.is_cellphonedb_interactor(interaction, ('_x', '_y')), axis=1)]
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
        zip_file.writestr('gene_synonym_to_gene_name.csv', gene_synonym_to_gene_name_db_df.to_csv(index=False, sep=',').encode('utf-8'))

    file_suffix = file_utils.get_timestamp_suffix()
    file_path = os.path.join(target_dir,'cellphonedb_{}.zip'.format(file_suffix))
    with open(file_path, 'wb') as f:
        f.write(zip_buffer.getvalue())
    print("Created {} successfully".format(file_path))

def download_database(target_dir, cpdb_version):
    download_released_files(target_dir, cpdb_version, "cellphonedb.zip|_input|sources\/uniprot_synonyms")

def download_released_files(target_dir, cpdb_version, regex):
    r = urllib.request.urlopen('https://github.com/prete/cellphonedb-data/archive/refs/tags/{}.zip'.format(cpdb_version))
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

def getDFs(gene_input=None, protein_input=None, complex_input=None, interaction_input=None, gene_synonyms_input=None):
    dfs = {}
    dfs['gene_input'] = file_utils.read_data_table_from_file(gene_input)
    dfs['protein_input'] = file_utils.read_data_table_from_file(protein_input)
    dfs['complex_input'] = file_utils.read_data_table_from_file(complex_input)
    dfs['interaction_input'] = file_utils.read_data_table_from_file(interaction_input)
    dfs['gene_synonyms_input'] = file_utils.read_data_table_from_file(gene_synonyms_input)
    return dfs

def run_sanity_tests(dataDFs):
    data_errors_found = False
    protein_db_df = dataDFs['protein_input']
    complex_db_df = dataDFs['complex_input']
    gene_db_df = dataDFs['gene_input']
    interaction_db_df = dataDFs['interaction_input']

    # 1. Report any uniprot accessions that map to multiple gene_names
    gene_names_uniprot_df = gene_db_df[['gene_name', 'uniprot']].copy()
    gene_names_uniprot_df.drop_duplicates(inplace=True)
    dups = gene_names_uniprot_df[gene_names_uniprot_df['uniprot'].duplicated() == True]
    if not dups.empty:
        # data_errors_found = True
        print("WARNING: The following UniProt ids map to multiple gene names (it is expected that they should map to only one):")
        print(", ".join(dups['uniprot'].tolist()))

    # 2. Warn about complex name duplicates in complex_db_df
    test_complex_db_df = complex_db_df.set_index('complex_name')
    if test_complex_db_df.index.has_duplicates:
        print("WARNING: complex_input.csv has the following duplicates:")
        print("\n".join(complex_db_df[test_complex_db_df.index.duplicated(keep='first')]['complex_name'].tolist()) + "\n")

    # 3. Report complexes with (possibly) different names, but with the same uniprot
    # accession participants (though not necessarily in the same order - hence the use of set below)
    # NB. Use set below as we don't care about the order of participants when looking for duplicates
    # NB. Report duplicate complexes _only if_ at least one duplicate's complex_db_df['version']
    # does not start with CORE_CELLPHONEDB_DATA)
    participants_set_to_complex_names = {}
    participants_set_to_data_sources = {}
    for row in complex_db_df[['complex_name','version'] + PROTEIN_COLUMN_NAMES].itertuples(index=False):
        participants_set = frozenset([i for i in row[2:] if str(i) != 'nan'])
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

    # 4. Report interactions with (possibly) a different name, but with the same participants
    # (though not necessarily in the same order - hence the use of set below)
    partner_sets = [set([i for i in row]) for row in \
            interaction_db_df[['partner_a','partner_b']].itertuples(index=False)]
    # Find duplicate sets of partners
    seen = set()
    duplicate_partner_sets = [x for x in partner_sets if x in seen or seen.add(frozenset(x))]
    if duplicate_partner_sets:
        # data_errors_found = True
        print("WARNING: The following sets of interaction partners appear in multiple rows of interaction_input.csv file:")
        for dup in set([frozenset(x) for x in duplicate_partner_sets]):
            print(','.join(dup))
        print()

    # 5. Warn about uniprot accession duplicates in protein_db_df
    test_protein_db_df = protein_db_df.set_index('uniprot')
    if test_protein_db_df.index.has_duplicates:
        print("WARNING: protein_input.csv has the following UniProt accession duplicates:")
        print("\n".join(protein_db_df[test_protein_db_df.index.duplicated(keep='first')]['uniprot'].tolist()) + "\n")

    # 6. Warn the user if some complexes don't participate in any interactions
    all_complexes_set = set(complex_db_df['complex_name'].tolist())
    interaction_participants_set = set (interaction_db_df['partner_a'].tolist() + interaction_db_df['partner_b'].tolist())
    orphan_complexes = all_complexes_set - interaction_participants_set
    if orphan_complexes:
        print("WARNING: The following complexes are not found in interaction_input.txt:")
        print ("\n".join(orphan_complexes))
    print()

    # 7. Warn the user if some proteins don't participate in any interactions directly, or are part some complex in orphan_complexes
    all_proteins_set = set(protein_db_df['uniprot'].tolist())
    proteins_in_complexes_participating_in_interactions = []
    for colName in PROTEIN_COLUMN_NAMES:
        proteins_in_complexes_participating_in_interactions += complex_db_df[~complex_db_df['complex_name'].isin(orphan_complexes)][colName].tolist()
    orphan_proteins = all_proteins_set - interaction_participants_set - set(proteins_in_complexes_participating_in_interactions)
    if orphan_proteins:
        print("WARNING: The following proteins are not found in interaction_input.txt (either directly or via complexes they are part of):")
        print ("\n".join(orphan_proteins))

    # 8. Warn the user if some interactions contain interactors that are neither in complex_input.csv or protein_input.csv
    unknown_interactors = set()
    for col in ['partner_a', 'partner_b']:
        aux_df = pd.merge(interaction_db_df, protein_db_df, left_on=col, right_on='uniprot',how='outer')
        unknown_interactor_proteins = set(aux_df[pd.isnull(aux_df['uniprot'])][col].tolist())
        aux_df = pd.merge(interaction_db_df, complex_db_df, left_on=col, right_on='complex_name', how='outer')
        unknown_interactor_complexes = set(aux_df[pd.isnull(aux_df['complex_name'])][col].tolist())
        unknown_interactors = unknown_interactors.union(unknown_interactor_proteins.intersection(unknown_interactor_complexes))
    if unknown_interactors:
        print("WARNING: The following interactors in interaction_input.txt could not be found in either protein_input.csv or complex_indput.csv:")
        print("\n".join(sorted(unknown_interactors)) + "\n")

    # 9. Warn if some complexes contain proteins not in protein_input.csv
    unknown_proteins = set()
    for col in PROTEIN_COLUMN_NAMES:
        aux_df = pd.merge(complex_db_df, protein_db_df, left_on=col, right_on='uniprot', how='outer')
        unknown_complex_proteins = set(aux_df[pd.isnull(aux_df['uniprot']) & ~pd.isnull(aux_df[col])][col].tolist())
        unknown_proteins = unknown_proteins.union(unknown_complex_proteins)
    if unknown_proteins:
        print("WARNING: The following proteins in complex_input.txt could not be found in protein_input.csv:")
        print("\n".join(sorted(unknown_proteins)) + "\n")

    # 10. Warn if some proteins in protein_input.csv are not in gene_input.csv
    proteins = set(protein_db_df['uniprot'].tolist())
    unknown_proteins = proteins.difference(set(gene_db_df['uniprot'].tolist()))
    if unknown_proteins:
        print("WARNING: The following proteins in protein_input.txt could not be found in gene_input.csv:")
        print("\n".join(sorted(unknown_proteins)) + "\n")

    if data_errors_found:
        raise DatabaseCreationException