import unittest
import os
import fnmatch
import sys
import re
import zipfile
from cellphonedb.utils import file_utils, search_utils, db_utils, db_releases_utils
from cellphonedb.src.core.methods import cpdb_analysis_method, cpdb_statistical_analysis_method, cpdb_degs_analysis_method
from cellphonedb.src.core.preprocessors import method_preprocessors
from cellphonedb.src.core.utils import subsampler

RELEASED_VERSION="v4.1.0"
script_dir = os.path.dirname(__file__)
target_dir = os.path.join(script_dir, "test_data")
GENERATED_CPDB_PATTERN="cellphonedb_*.zip"
downloaded_db_dir = os.path.join(target_dir, "downloaded_db")
generated_db_dir = os.path.join(target_dir, "generated_db")
output_dir = os.path.join(target_dir, "out")
CELLPHONEDB_FILES = ['complex_composition_table.csv','gene_synonym_to_gene_name.csv','interaction_table.csv','complex_table.csv','gene_table.csv','multidata_table.csv']
TEST_FILES_PATH = os.path.join(os.path.join("..", os.path.join("..", os.path.join("..","example_data"))))
cpdb_file_path = os.path.join(target_dir, "cellphonedb.zip")
META_FILE_PATH = os.path.join(TEST_FILES_PATH, 'test_meta.txt')
COUNTS_FILE_PATH = os.path.join(TEST_FILES_PATH, 'test.h5ad')
MICROENVS_FILE_PATH = os.path.join(TEST_FILES_PATH, 'test_microenviroments.txt')
DEGS_FILE_PATH = os.path.join(TEST_FILES_PATH, 'test_degs.txt')

class UnitTests(unittest.TestCase):

    def setUp(self):
        if not os.path.isdir(target_dir):
            os.mkdir(target_dir)
        if not os.path.isdir(output_dir):
            os.mkdir(output_dir)
        if not os.path.isdir(downloaded_db_dir):
            os.mkdir(downloaded_db_dir)
        if not os.path.isdir(generated_db_dir):
            os.mkdir(generated_db_dir)

    def test_download_database(self):
        db_utils.download_database(target_dir, RELEASED_VERSION)

    def test_createDatabase(self):
        for f in os.listdir(target_dir):
            if re.search("cellphonedb_.*\\.zip", f):
                os.remove(os.path.join(target_dir, f))
        db_utils.create_db(target_dir)
        generated_db_file = None
        for f in os.listdir(target_dir):
            if fnmatch.fnmatch(f, GENERATED_CPDB_PATTERN):
                generated_db_file = f
                break
        assert generated_db_file is not None

    def test_compare_downloaded_created_dbs(self):
        generated_db_file = None
        for f in os.listdir(target_dir):
            if fnmatch.fnmatch(f, GENERATED_CPDB_PATTERN):
                generated_db_file = f
                break
        with zipfile.ZipFile(os.path.join(target_dir, generated_db_file), 'r') as zip_ref:
            zip_ref.extractall(generated_db_dir)
        with zipfile.ZipFile(os.path.join(target_dir, "cellphonedb.zip"), 'r') as zip_ref:
            zip_ref.extractall(downloaded_db_dir)
        for f in CELLPHONEDB_FILES:
            generated_db_f = os.path.join(target_dir, generated_db_dir, f)
            downloaded_db_f = os.path.join(target_dir, downloaded_db_dir, f)
            with open(generated_db_f, 'r') as fp:
                for generated_count, line in enumerate(fp):
                    pass
            with open(downloaded_db_f, 'r') as fp:
                for downloaded_count, line in enumerate(fp):
                    pass
            print("Comparing {} between generated and downloaded DB".format(f))
            assert generated_count == downloaded_count, "The number of lines in {} differs between generated and downloaded DB".format(f)

    def test_search_db(self):
        cpdb_file_path = os.path.join(target_dir, "cellphonedb.zip")
        results, complex_name2proteins, protein2Info, complex2Info, resource2Complex2Acc, proteinAcc2Name = \
            search_utils.search('ENSG00000134780,integrin_a10b1_complex', cpdb_file_path)
        assert len(results) > 0

    def test_get_remote_release_versions(self):
        result = db_releases_utils.get_remote_database_versions_html(True, 4.0)
        assert 'db_releases_html_table' in result

    def test_basic_method(self):
        assert os.path.exists(TEST_FILES_PATH), "{} does not exist".format(TEST_FILES_PATH)
        analysis_result = \
            cpdb_analysis_method.call(
                cpdb_file_path=cpdb_file_path,
                meta_file_path=META_FILE_PATH,
                counts_file_path=COUNTS_FILE_PATH,
                counts_data='ensembl',
                output_path=output_dir,
                microenvs_file_path=MICROENVS_FILE_PATH,
                separator="|",
                threshold=0.1,
                result_precision=3,
                debug=False,
                output_suffix=None,
                score_interactions=True,
                threads=1)
        assert analysis_result is not None
        assert 'deconvoluted' in analysis_result
        assert 'deconvoluted_percents' in analysis_result
        assert 'means_result' in analysis_result
        assert 'interaction_scores' in analysis_result
        self.assertFalse(analysis_result['deconvoluted'].empty, 'deconvoluted dataframe is empty')
        self.assertFalse(analysis_result['deconvoluted_percents'].empty, 'deconvoluted_percents dataframe is empty')
        self.assertFalse(analysis_result['means_result'].empty, 'means_result dataframe is empty')
        self.assertFalse(analysis_result['interaction_scores'].empty, 'interaction_scores dataframe is empty')
        # TODO: Test CellSign output when test data is available

    def test_statistical_method(self):
        assert os.path.exists(TEST_FILES_PATH), "{} does not exist".format(TEST_FILES_PATH)
        analysis_result = \
            cpdb_statistical_analysis_method.call(
                cpdb_file_path = cpdb_file_path,
                meta_file_path = META_FILE_PATH,
                counts_file_path = COUNTS_FILE_PATH,
                counts_data = 'ensembl',
                output_path = output_dir,
                microenvs_file_path = MICROENVS_FILE_PATH,
                active_tfs_file_path = None,
                iterations = 1000,
                threshold = 0.1,
                threads = 1,
                debug_seed = -1,
                result_precision = 3,
                pvalue = 0.05,
                subsampling = False,
                subsampling_log = False,
                subsampling_num_pc = 100,
                subsampling_num_cells = None,
                separator = '|',
                debug = False,
                output_suffix = None,
                score_interactions = True)
        assert analysis_result is not None
        assert 'deconvoluted' in analysis_result
        assert 'deconvoluted_percents' in analysis_result
        assert 'means' in analysis_result
        assert 'pvalues' in analysis_result
        assert 'significant_means' in analysis_result
        assert 'interaction_scores' in analysis_result
        assert 'CellSign_active_interactions' in analysis_result
        assert 'CellSign_active_interactions_deconvoluted' in analysis_result
        self.assertFalse(analysis_result['deconvoluted'].empty, 'deconvoluted dataframe is empty')
        self.assertFalse(analysis_result['deconvoluted_percents'].empty, 'deconvoluted_percents dataframe is empty')
        self.assertFalse(analysis_result['means'].empty, 'means dataframe is empty')
        self.assertFalse(analysis_result['pvalues'].empty, 'pvalues dataframe is empty')
        self.assertFalse(analysis_result['significant_means'].empty, 'significant_means dataframe is empty')
        self.assertFalse(analysis_result['interaction_scores'].empty, 'interaction_scores dataframe is empty')
        # TODO: Test CellSign output when test data becomes available

    def test_deg_method(self):
        assert os.path.exists(TEST_FILES_PATH), "{} does not exist".format(TEST_FILES_PATH)
        analysis_result = \
            cpdb_degs_analysis_method.call(
                cpdb_file_path=cpdb_file_path,
                meta_file_path=META_FILE_PATH,
                counts_file_path=COUNTS_FILE_PATH,
                degs_file_path=DEGS_FILE_PATH,
                counts_data='ensembl',
                microenvs_file_path=MICROENVS_FILE_PATH,
                active_tfs_file_path=None,
                threshold=0.1,
                result_precision=3,
                separator='|',
                debug=False,
                output_path=output_dir,
                output_suffix=None,
                score_interactions=True,
                threads=1)
        assert analysis_result is not None
        assert 'deconvoluted' in analysis_result
        assert 'deconvoluted_percents' in analysis_result
        assert 'means' in analysis_result
        assert 'relevant_interactions' in analysis_result
        assert 'significant_means' in analysis_result
        assert 'interaction_scores' in analysis_result
        assert 'CellSign_active_interactions' in analysis_result
        assert 'CellSign_active_interactions_deconvoluted' in analysis_result
        self.assertFalse(analysis_result['deconvoluted'].empty, 'deconvoluted dataframe is empty')
        self.assertFalse(analysis_result['deconvoluted_percents'].empty, 'deconvoluted_percents dataframe is empty')
        self.assertFalse(analysis_result['means'].empty, 'means dataframe is empty')
        self.assertFalse(analysis_result['relevant_interactions'].empty, 'relevant_interactions dataframe is empty')
        self.assertFalse(analysis_result['significant_means'].empty, 'significant_means dataframe is empty')
        self.assertFalse(analysis_result['interaction_scores'].empty, 'interaction_scores dataframe is empty')
        # TODO: Test CellSign output when test data becomes available

if __name__ == "__main__":
    unittest.main()