import sys
import pandas as pd
import numpy as np
import anndata
import os
from typing import Tuple
import io
from utils import utils, generate_input_files, database_version_manager, search_utils, db_utils
from utils.utils import dbg
import urllib.request, urllib.error, urllib.parse

from src.core.methods import cpdb_analysis_method, cpdb_statistical_analysis_method, cpdb_degs_analysis_method
from src.core.preprocessors import method_preprocessors
from src.exceptions.ParseCountsException import ParseCountsException

KEY2USER_TEST_FILE = {'counts' : 'test_counts.txt', 'meta': 'test_meta.txt', \
                         'microenvs' : 'test_microenviroments.txt', 'degs' : 'test_degs.txt'}

RELEASED_VERSION="v5.0.0"
CPDB_ROOT = os.path.join(os.path.expanduser('~'),".cpdb")

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
    adata = utils.read_h5ad(os.path.join(user_dir_root,'user_files',h5ad_fn))

    counts = adata.to_df().T
    counts.columns = adata.obs['sample']

    meta_str = adata.uns['meta']
    if meta_str:
        raw_meta = utils._read_data(io.StringIO(meta_str), separator='\t', index_column_first=False, dtype=None,
                                    na_values=None, compression=None)
        meta = method_preprocessors.meta_preprocessor(raw_meta)
        # adata and the merge below are needed for plotting results via ktplotspy
        adata.obs = pd.merge(adata.obs, meta, left_on="sample", right_on="cell")
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
    return adata, counts, raw_meta, meta, microenvs, degs

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
    interactions, genes, complex_composition, complex_expanded = \
        db_utils.get_interactions_genes_complex(user_dir_root, db_version)
    counts, raw_meta, meta, microenvs, degs = get_user_files(user_dir_root, \
                                                                        counts_fn=KEY2USER_TEST_FILE['counts'],
                                                                        meta_fn=KEY2USER_TEST_FILE['meta'], \
                                                                        microenvs_fn=KEY2USER_TEST_FILE['microenvs'],
                                                                        degs_fn=KEY2USER_TEST_FILE['degs'])
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
                                              iterations = 1000,
                                              threshold = 0.1,
                                              threads = 4,
                                              debug_seed = -1,
                                              result_precision = 3,
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
                                    iterations = 1000,
                                    threshold = 0.1,
                                    threads = 4,
                                    debug_seed = -1,
                                    result_precision = 3,
                                    separator = '|',
                                    debug = False,
                                    output_path = '')

    # dbg("means:", means.index, means.columns, means.info)
    # dbg("significant_means:", significant_means.index, significant_means.columns, significant_means.info)
    # dbg("deconvoluted:", deconvoluted.index, deconvoluted.columns, deconvoluted.info)
    # dbg("relevant_interactions:", relevant_interactions.index, relevant_interactions.columns, relevant_interactions.info)


def download_source_files(user_dir_root, db_version):
    sources_path = os.path.join(db_utils.get_db_data_path(user_dir_root, db_version), "sources")
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
        db_utils.download_input_files(CPDB_ROOT, RELEASED_VERSION)
        db_dir = db_utils.get_db_path(CPDB_ROOT, RELEASED_VERSION)
        data_dir = os.path.join(db_dir, "data")
        gene_input_path = os.path.join(data_dir, "gene_input.csv")
        protein_input_path = os.path.join(data_dir, "protein_input.csv")
        complex_input_path = os.path.join(data_dir, "complex_input.csv")
        interaction_input_path = os.path.join(data_dir, "interaction_input.csv")
        db_utils.create_db(CPDB_ROOT, RELEASED_VERSION, \
                           gene_input=gene_input_path, protein_input=protein_input_path, complex_input=complex_input_path,
                           interaction_input=interaction_input_path)
    elif arg == "dbd":
        # database_version_manager.download_database("latest")
        database_version_manager.download_database("v4.0.0")
    elif arg == 'c':
        convert_to_h5ad(CPDB_ROOT)
    elif arg == 's':
        search_utils.search('ENSG00000134780,integrin_a10b1_complex', CPDB_ROOT, RELEASED_VERSION)
    elif arg == 'g':
        generate_input_files.generate_all(CPDB_ROOT, RELEASED_VERSION, \
                                          user_complex=None, user_interactions=None, user_interactions_only=False)

    else:
        print("Arguments can be a (perform analysis) or db (create database)")

