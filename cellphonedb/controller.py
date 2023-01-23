import sys
import pandas as pd
import anndata
import os
from typing import Tuple
import io
import urllib.request, urllib.error, urllib.parse
from cellphonedb.utils import utils, generate_input_files, database_version_manager, search_utils, db_utils, db_releases_utils
from cellphonedb.utils.utils import dbg
from cellphonedb.src.core.methods import cpdb_analysis_method, cpdb_statistical_analysis_method, cpdb_degs_analysis_method
from cellphonedb.src.core.preprocessors import method_preprocessors, counts_preprocessors
from cellphonedb.src.core.utils import subsampler
import time

KEY2USER_TEST_FILE = {'counts' : 'test_counts.txt', 'meta': 'test_meta.txt', \
                         'microenvs' : 'test_microenviroments.txt', 'degs' : 'test_degs.txt'}

RELEASED_VERSION="v5.0.0"
CPDB_ROOT = os.path.join(os.path.expanduser('~'),".cpdb")

def get_user_files(user_dir_root, \
        counts_fn=KEY2USER_TEST_FILE['counts'], meta_fn=KEY2USER_TEST_FILE['meta'], microenvs_fn=None, degs_fn=None) \
        -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """

    Parameters
    ----------
    user_dir_root: str
        The directory in which user stores CellphoneDB files
    counts_fn
        Path to the user's counts file, exemplified by \
        https://github.com/ventolab/CellphoneDB/blob/bare-essentials/example_data/test_counts.txt
    meta_fn
        Path to the user's meta file, exemplified by \
        https://github.com/ventolab/CellphoneDB/blob/bare-essentials/example_data/test_meta.txt
    microenvs_fn
        Path to the user's microenvironments file, exemplified by \
        https://github.com/ventolab/CellphoneDB/blob/bare-essentials/example_data/test_microenviroments.txt
    degs_fn
        Path to the user's differentially expresses genes (DEGs) file, exemplified by \
        https://github.com/ventolab/CellphoneDB/blob/bare-essentials/example_data/test_degs.txt

    Returns
    -------
    Tuple
        - counts: pd.DataFrame
        - raw_meta: pd.DataFrame
        - meta: pd.DataFrame
        - microenvs: pd.DataFrame
        - degs: pd.DataFrame

    """
    loaded_user_files=[]
    user_files_path = os.path.join(user_dir_root,"user_files")
    # Read user files
    counts = utils.read_data_table_from_file(os.path.join(user_files_path, counts_fn),
                                             index_column_first=True)
    loaded_user_files.append(counts_fn)
    raw_meta = utils.read_data_table_from_file(os.path.join(user_files_path, meta_fn),
                                               index_column_first=False)
    meta = method_preprocessors.meta_preprocessor(raw_meta)
    loaded_user_files.append(meta_fn)

    if microenvs_fn:
        microenvs = utils.read_data_table_from_file(os.path.join(user_files_path, microenvs_fn))
        loaded_user_files.append(microenvs_fn)
    else:
        microenvs = pd.DataFrame()

    if degs_fn:
        degs = utils.read_data_table_from_file(os.path.join(user_files_path, degs_fn))
        loaded_user_files.append(degs_fn)
    else:
        degs = pd.DataFrame()

    print("The following user files were loaded successfully:")
    for fn in loaded_user_files:
        print(fn)

    return counts, raw_meta, meta, microenvs, degs

def get_user_file(user_dir_root, h5ad_fn='test.h5ad') \
        -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """

    Parameters
    ----------
    user_dir_root: str
        The directory in which user stores CellphoneDB files
    h5ad_fn: str
        Path to the user's file in https://broadinstitute.github.io/wot/file_formats/#h5ad format. \
        The file is assumed to contain meta, microenvironments and differentially expressed genes (DEGs) data \
        in the uns portion of the h5ad file, under 'meta', 'microenvs' and 'degs' keys respectively, and each \
        follow the respective format in: \
        https://github.com/ventolab/CellphoneDB/blob/bare-essentials/example_data/test_meta.txt \
        https://github.com/ventolab/CellphoneDB/blob/bare-essentials/example_data/test_microenviroments.txt \
        https://github.com/ventolab/CellphoneDB/blob/bare-essentials/example_data/test_degs.txt

    Returns
    -------
    Tuple
        - counts: pd.DataFrame
        - raw_meta: pd.DataFrame
        - meta: pd.DataFrame
        - microenvs: pd.DataFrame
        - degs: pd.DataFrame

    """
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
    counts = counts_preprocessors.counts_preprocessor(counts, meta)
    ss = subsampler.Subsampler(log=False, num_pc=100, num_cells=None, verbose=False, debug_seed=None)
    if ss is not None:
        counts = ss.subsample(counts)
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
    elif arg == 'rel':
        db_releases_utils.get_remote_database_versions_html()
    elif arg == 'te':
        # Run statistical and deg analyses for endometrium example - for the purpose of comparing
        # results to old CellphoneDB or ones after the new code optiisations
        root_dir = os.path.join(CPDB_ROOT, 'tests', 'data', 'examples')
        dbversion = "v4.0.0"
        interactions, genes, complex_composition, complex_expanded = \
            db_utils.get_interactions_genes_complex(CPDB_ROOT, RELEASED_VERSION)
        adata = utils.read_h5ad(os.path.join(root_dir, 'endometrium_example_counts.h5ad'))
        counts = adata.to_df().T
        raw_meta = utils.read_data_table_from_file(
            os.path.join(root_dir, 'endometrium_example_meta.tsv'))
        meta = method_preprocessors.meta_preprocessor(raw_meta)
        microenvs = utils.read_data_table_from_file(
            os.path.join(root_dir, 'endometrium_example_microenviroments.tsv'))
        ss = subsampler.Subsampler(log=False, num_pc=100, num_cells=None, verbose=False, debug_seed=None)
        if ss is not None:
            counts = ss.subsample(counts)
        t0 = time.time()
        deconvoluted, means, pvalues, significant_means = \
        cpdb_statistical_analysis_method.call(meta,
                                              counts,
                                              'hgnc_symbol',
                                              interactions,
                                              genes,
                                              complex_expanded,
                                              complex_composition,
                                              microenvs=microenvs,
                                              iterations=1000,
                                              threshold=0.1,
                                              threads=4,
                                              debug_seed=-1,
                                              result_precision=3,
                                              pvalue=1,
                                              separator='|',
                                              debug=False)
        print("Statistical method took: ", time.time() - t0, "seconds to complete")
        output_path=os.path.join(root_dir, 'stat_new')
        utils.write_to_file(means, 'means.txt', output_path)
        utils.write_to_file(pvalues, 'pvalues.txt', output_path)
        utils.write_to_file(significant_means, 'significant_means.txt', output_path)
        utils.write_to_file(deconvoluted, 'deconvoluted.txt', output_path)
        """
        degs = utils.read_data_table_from_file(os.path.join(root_dir, 'endometrium_example_DEGs.tsv'))
        output_path = os.path.join(root_dir, 'deg_new')
        deconvoluted, means, relevant_interactions, significant_means = \
            cpdb_degs_analysis_method.call(meta,
                                           counts,
                                           degs,
                                           'hgnc_symbol',
                                           interactions,
                                           genes,
                                           complex_expanded,
                                           complex_composition,
                                           microenvs=microenvs,
                                           iterations=1000,
                                           threshold=0.1,
                                           threads=4,
                                           debug_seed=-1,
                                           result_precision=3,
                                           separator='|',
                                           debug=False,
                                           output_path='')
        utils.write_to_file(means, 'means.txt', output_path)
        utils.write_to_file(relevant_interactions, 'relevant_interactions.txt', output_path)
        utils.write_to_file(significant_means, 'significant_means.txt', output_path)
        utils.write_to_file(deconvoluted, 'deconvoluted.txt', output_path)
        """
    elif arg == 'to':
        # Run statistical and deg analyses for ovarian example - for the purpose of comparing results to old CellphoneDB
        root_dir = os.path.join(CPDB_ROOT, 'tests', 'data', 'bug_2_ovary')
        dbversion = "v4.0.0"
        interactions, genes, complex_composition, complex_expanded = \
            db_utils.get_interactions_genes_complex(CPDB_ROOT, RELEASED_VERSION)
        adata = utils.read_h5ad(os.path.join(root_dir, 'granulosa_normloqTransformed.h5ad'))
        counts = adata.to_df().T
        raw_meta = utils.read_data_table_from_file(os.path.join(root_dir, 'ovarian_meta.tsv'))
        meta = method_preprocessors.meta_preprocessor(raw_meta)
        microenvs = utils.read_data_table_from_file(os.path.join(root_dir, 'ovarian_microenviroment.tsv'))
        deconvoluted, means, pvalues, significant_means = \
            cpdb_statistical_analysis_method.call(meta,
                                                  counts,
                                                  'gene_name',
                                                  interactions,
                                                  genes,
                                                  complex_expanded,
                                                  complex_composition,
                                                  microenvs=microenvs,
                                                  iterations=1000,
                                                  threshold=0.1,
                                                  threads=4,
                                                  debug_seed=-1,
                                                  result_precision=3,
                                                  pvalue=1,
                                                  separator='|',
                                                  debug=False)
        output_path = os.path.join(root_dir, 'stat_new')
        utils.write_to_file(means, 'means.txt', output_path)
        utils.write_to_file(pvalues, 'pvalues.txt', output_path)
        utils.write_to_file(significant_means, 'significant_means.txt', output_path)
        utils.write_to_file(deconvoluted, 'deconvoluted.txt', output_path)
        degs = utils.read_data_table_from_file(root_dir, 'DEGs.tsv')
        output_path = os.path.join(root_dir, 'deg_new')
        deconvoluted, means, relevant_interactions, significant_means = \
            cpdb_degs_analysis_method.call(meta,
                                           counts,
                                           degs,
                                           'gene_name',
                                           interactions,
                                           genes,
                                           complex_expanded,
                                           complex_composition,
                                           microenvs=microenvs,
                                           iterations=1000,
                                           threshold=0.1,
                                           threads=4,
                                           debug_seed=-1,
                                           result_precision=3,
                                           separator='|',
                                           debug=False,
                                           output_path='')
        utils.write_to_file(means, 'means.txt', output_path)
        utils.write_to_file(relevant_interactions, 'relevant_interactions.txt', output_path)
        utils.write_to_file(significant_means, 'significant_means.txt', output_path)
        utils.write_to_file(deconvoluted, 'deconvoluted.txt', output_path)
    else:
        print("Arguments can be a (perform analysis) or db (create database)")

