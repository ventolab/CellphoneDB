import sys
import pandas as pd
import anndata
import os
from cellphonedb.utils import file_utils, generate_input_files, search_utils, db_utils, db_releases_utils
from cellphonedb.src.core.methods import cpdb_analysis_method, cpdb_statistical_analysis_method, cpdb_degs_analysis_method
from cellphonedb.src.core.preprocessors import method_preprocessors, counts_preprocessors
from cellphonedb.src.core.utils import subsampler
import time

RELEASED_VERSION="v4.1.0"
CPDB_ROOT = os.path.join(os.path.expanduser('~'),".cpdb")
USER_FILES_DIR = os.path.join(CPDB_ROOT, "user_files")

KEY2USER_TEST_FILE = {'counts' : os.path.join(USER_FILES_DIR, 'test_counts.txt'), \
                      'meta': os.path.join(USER_FILES_DIR, 'test_meta.txt'), \
                      'microenvs' : os.path.join(USER_FILES_DIR, 'test_microenviroments.txt'), \
                      'degs' : os.path.join(USER_FILES_DIR, 'test_degs.txt')}

def testrun_analyses(cpdb_file_path):
    interactions, genes, complex_composition, complex_expanded = \
        db_utils.get_interactions_genes_complex(cpdb_dir)
    counts, raw_meta, meta, microenvs, degs = file_utils.get_user_files(counts_fp=KEY2USER_TEST_FILE['counts'],
                                                             meta_fp=KEY2USER_TEST_FILE['meta'],
                                                             microenvs_fp=KEY2USER_TEST_FILE['microenvs'],
                                                             degs_fp=KEY2USER_TEST_FILE['degs'])
    # ************ Call analysis method
    means, significant_means, deconvoluted = cpdb_analysis_method.call(
        cpdb_dir,
        meta,
        counts,
        'ensembl',
        microenvs=microenvs,
        # Does not store results files in cellphonedb/out/debug_intermediate.pkl
        debug=False,
        output_path=os.path.join('.','out')
        );
    # dbg("means:", means.index, means.columns, means.info)
    # dbg("significant_means:", significant_means.index, significant_means.columns, significant_means.info)
    # dbg("deconvoluted:", deconvoluted.index, deconvoluted.columns, deconvoluted.info)

    # ************ Call statistical analysis method
    ss = subsampler.Subsampler(log=False, num_pc=100, num_cells=None, verbose=False, debug_seed=None)
    if ss is not None:
        counts = ss.subsample(counts)
    deconvoluted, means, pvalues, significant_means = \
        cpdb_statistical_analysis_method.call(cpdb_dir,
                                              meta,
                                              counts,
                                              'ensembl',
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
        cpdb_degs_analysis_method.call(cpdb_dir,
                                    meta,
                                    counts,
                                    degs,
                                    'ensembl',
                                    microenvs=microenvs,
                                    threshold = 0.1,
                                    debug_seed = -1,
                                    result_precision = 3,
                                    separator = '|',
                                    debug = False,
                                    output_path = '')

    # dbg("means:", means.index, means.columns, means.info)
    # dbg("significant_means:", significant_means.index, significant_means.columns, significant_means.info)
    # dbg("deconvoluted:", deconvoluted.index, deconvoluted.columns, deconvoluted.info)
    # dbg("relevant_interactions:", relevant_interactions.index, relevant_interactions.columns, relevant_interactions.info)


def convert_to_h5ad():
    counts, raw_meta, meta, microenvs, degs = file_utils.get_user_files(counts_fp=KEY2USER_TEST_FILE['counts'],
                                                             meta_fp=KEY2USER_TEST_FILE['meta'],
                                                             microenvs_fp=KEY2USER_TEST_FILE['microenvs'],
                                                             degs_fp='test_degs.txt')
    obs = pd.DataFrame()
    obs.index = counts.columns
    # dataframe for annotating the variables
    var = pd.DataFrame(index=counts.index)
    adata = anndata.AnnData(counts.T.values, obs=obs, var=var, dtype='float64')
    outputPath = os.path.join(USER_FILES_DIR,"test.h5ad")
    adata.write(outputPath)

if __name__ == '__main__':
    cpdb_dir = db_utils.get_db_path(CPDB_ROOT, RELEASED_VERSION)
    cpdb_file_path = os.path.join(cpdb_dir, "cellphonedb.zip")
    arg = sys.argv[1]
    if arg == 'a':
        testrun_analyses(cpdb_file_path)
    elif arg == 'db':
        db_utils.download_released_files(cpdb_dir, RELEASED_VERSION, "cellphonedb.zip")
        db_utils.create_db(cpdb_dir)
    elif arg == 'c':
        convert_to_h5ad(os.path.join(CPDB_ROOT, "user_files"))
    elif arg == 's':
        search_utils.search('ENSG00000134780,integrin_a10b1_complex', cpdb_dir)
    elif arg == 'g':
        generate_input_files.generate_all(cpdb_dir, RELEASED_VERSION, \
                                          user_complex=None, user_interactions=None, user_interactions_only=False)
    elif arg == 'rel':
        db_releases_utils.get_remote_database_versions_html()
    elif arg == 'te':
        # Run statistical and deg analyses for endometrium example - for the purpose of comparing
        # results to old CellphoneDB or ones after the new code optiisations
        root_dir = os.path.join(CPDB_ROOT, 'tests', 'data', 'examples')
        dbversion = "v4.0.0"
        interactions, genes, complex_composition, complex_expanded = \
            db_utils.get_interactions_genes_complex(cpdb_file_path)
        adata = file_utils.read_h5ad(os.path.join(root_dir, 'endometrium_example_counts.h5ad'))
        counts = adata.to_df().T
        raw_meta = file_utils.read_data_table_from_file(
            os.path.join(root_dir, 'endometrium_example_meta.tsv'))
        meta = method_preprocessors.meta_preprocessor(raw_meta)
        microenvs = file_utils.read_data_table_from_file(
            os.path.join(root_dir, 'endometrium_example_microenviroments.tsv'))
        ss = subsampler.Subsampler(log=False, num_pc=100, num_cells=None, verbose=False, debug_seed=None)
        if ss is not None:
            counts = ss.subsample(counts)
        t0 = time.time()
        deconvoluted, means, pvalues, significant_means = \
        cpdb_statistical_analysis_method.call(cpdb_dir,
                                              meta,
                                              counts,
                                              'hgnc_symbol',
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
        file_utils.write_to_file(means, 'means.txt', output_path)
        file_utils.write_to_file(pvalues, 'pvalues.txt', output_path)
        file_utils.write_to_file(significant_means, 'significant_means.txt', output_path)
        file_utils.write_to_file(deconvoluted, 'deconvoluted.txt', output_path)
        """
        degs = file_utils.read_data_table_from_file(os.path.join(root_dir, 'endometrium_example_DEGs.tsv'))
        output_path = os.path.join(root_dir, 'deg_new')
        deconvoluted, means, relevant_interactions, significant_means = \
            cpdb_degs_analysis_method.call(cpdb_dir,
                                           meta,
                                           counts,
                                           degs,
                                           'hgnc_symbol',
                                           microenvs=microenvs,
                                           threshold=0.1,
                                           debug_seed=-1,
                                           result_precision=3,
                                           separator='|',
                                           debug=False,
                                           output_path='')
        file_utils.write_to_file(means, 'means.txt', output_path)
        file_utils.write_to_file(relevant_interactions, 'relevant_interactions.txt', output_path)
        file_utils.write_to_file(significant_means, 'significant_means.txt', output_path)
        file_utils.write_to_file(deconvoluted, 'deconvoluted.txt', output_path)
        """
    elif arg == 'to':
        # Run statistical and deg analyses for ovarian example - for the purpose of comparing results to old CellphoneDB
        root_dir = os.path.join(CPDB_ROOT, 'tests', 'data', 'bug_2_ovary')
        dbversion = "v4.0.0"
        interactions, genes, complex_composition, complex_expanded = \
            db_utils.get_interactions_genes_complex(cpdb_file_path)
        adata = file_utils.read_h5ad(os.path.join(root_dir, 'granulosa_normloqTransformed.h5ad'))
        counts = adata.to_df().T
        raw_meta = file_utils.read_data_table_from_file(os.path.join(root_dir, 'ovarian_meta.tsv'))
        meta = method_preprocessors.meta_preprocessor(raw_meta)
        microenvs = file_utils.read_data_table_from_file(os.path.join(root_dir, 'ovarian_microenviroment.tsv'))
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
        file_utils.write_to_file(means, 'means.txt', output_path)
        file_utils.write_to_file(pvalues, 'pvalues.txt', output_path)
        file_utils.write_to_file(significant_means, 'significant_means.txt', output_path)
        file_utils.write_to_file(deconvoluted, 'deconvoluted.txt', output_path)
        degs = file_utils.read_data_table_from_file(root_dir, 'DEGs.tsv')
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
                                           threshold=0.1,
                                           debug_seed=-1,
                                           result_precision=3,
                                           separator='|',
                                           debug=False,
                                           output_path='')
        file_utils.write_to_file(means, 'means.txt', output_path)
        file_utils.write_to_file(relevant_interactions, 'relevant_interactions.txt', output_path)
        file_utils.write_to_file(significant_means, 'significant_means.txt', output_path)
        file_utils.write_to_file(deconvoluted, 'deconvoluted.txt', output_path)
    else:
        print("Arguments can be a (perform analysis) or db (create database)")

