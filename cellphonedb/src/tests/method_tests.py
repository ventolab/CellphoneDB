import unittest
import os
import fnmatch
import re
import zipfile
from cellphonedb.utils import search_utils, db_utils, db_releases_utils
from cellphonedb.src.core.methods import cpdb_analysis_method, cpdb_statistical_analysis_method, cpdb_degs_analysis_method

RELEASED_VERSION = "v4.1.0"
# Ready for v5.0.0 data release: RELEASED_VERSION = "v5.0.0"
script_dir = os.path.dirname(__file__)
target_dir = os.path.join(script_dir, "test_data")
GENERATED_CPDB_PATTERN = "cellphonedb_*.zip"
downloaded_db_dir = os.path.join(target_dir, "downloaded_db")
DOWNLOADED_DB_FILE = os.path.join(downloaded_db_dir, "cellphonedb.zip")
generated_db_dir = os.path.join(target_dir, "generated_db")
output_dir = os.path.join(target_dir, "out")
CELLPHONEDB_FILES = ['complex_composition_table.csv', 'gene_synonym_to_gene_name.csv', 'interaction_table.csv',
                     'complex_table.csv', 'gene_table.csv', 'multidata_table.csv']
TEST_FILES_PATH = os.path.join(os.path.join("..", os.path.join("..", os.path.join("..", "example_data"))))
META_FILE_PATH = os.path.join(TEST_FILES_PATH, 'test_meta.txt')
COUNTS_FILE_PATH = os.path.join(TEST_FILES_PATH, 'test.h5ad')
MICROENVS_FILE_PATH = os.path.join(TEST_FILES_PATH, 'test_microenviroments.txt')
DEGS_FILE_PATH = os.path.join(TEST_FILES_PATH, 'test_degs.txt')
ACTIVE_TFS_PATH = os.path.join(TEST_FILES_PATH, 'test_active_tfs.txt')


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
        if not os.path.exists(DOWNLOADED_DB_FILE):
            db_utils.download_database(downloaded_db_dir, RELEASED_VERSION)
            db_utils.download_database(generated_db_dir, RELEASED_VERSION)

    def test_compare_downloaded_created_dbs(self):
        # First remove any previously generated DB file
        for f in os.listdir(generated_db_dir):
            if re.search("cellphonedb_.*\\.zip", f):
                os.remove(os.path.join(generated_db_dir, f))
        # Now generate a new DB file
        db_utils.create_db(generated_db_dir)
        # Find out the name of the newly generated DB file
        generated_db_file = None
        for f in os.listdir(generated_db_dir):
            if fnmatch.fnmatch(f, GENERATED_CPDB_PATTERN):
                generated_db_file = f
                break
        assert generated_db_file is not None
        # Now compare the number of lines in each file of the generated and downloaded DB in turn
        with zipfile.ZipFile(os.path.join(generated_db_dir, generated_db_file), 'r') as zip_ref:
            zip_ref.extractall(generated_db_dir)
        with zipfile.ZipFile(DOWNLOADED_DB_FILE, 'r') as zip_ref:
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
            assert generated_count == downloaded_count, \
                "The number of lines in {} differs between generated and downloaded DB".format(f)

    def test_search_db(self):
        results, complex_name2proteins, protein2Info, complex2Info, resource2Complex2Acc, proteinAcc2Name = \
            search_utils.search('ENSG00000134780,integrin_a10b1_complex', DOWNLOADED_DB_FILE)
        assert len(results) > 0

    def test_get_remote_release_versions(self):
        result = db_releases_utils.get_remote_database_versions_html(True, 4.0)
        assert 'db_releases_html_table' in result

    def test_basic_method(self):
        assert os.path.exists(TEST_FILES_PATH), "{} does not exist".format(TEST_FILES_PATH)
        analysis_result = \
            cpdb_analysis_method.call(
                cpdb_file_path=DOWNLOADED_DB_FILE,
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

    def test_statistical_method(self):
        assert os.path.exists(TEST_FILES_PATH), "{} does not exist".format(TEST_FILES_PATH)
        analysis_result = \
            cpdb_statistical_analysis_method.call(
                cpdb_file_path=DOWNLOADED_DB_FILE,
                meta_file_path=META_FILE_PATH,
                counts_file_path=COUNTS_FILE_PATH,
                counts_data='ensembl',
                output_path=output_dir,
                microenvs_file_path=MICROENVS_FILE_PATH,
                active_tfs_file_path=None,
                # Ready for v5.0.0 data release: active_tfs_file_path=ACTIVE_TFS_PATH,
                iterations=1000,
                threshold=0.1,
                threads=1,
                debug_seed=-1,
                result_precision=3,
                pvalue=0.05,
                subsampling=False,
                subsampling_log=False,
                subsampling_num_pc=100,
                subsampling_num_cells=None,
                separator='|',
                debug=False,
                output_suffix=None,
                score_interactions=True)
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
        # Ready for v5.0.0 data release: self.assertFalse(analysis_result['CellSign_active_interactions'].empty,
        # 'CellSign_active_interactions dataframe is empty')
        # Ready for v5.0.0 data release: self.assertFalse(analysis_result['CellSign_active_interactions_deconvoluted'].empty,
        # 'CellSign_active_interactions_deconvoluted dataframe is empty')

    def test_deg_method(self):
        assert os.path.exists(TEST_FILES_PATH), "{} does not exist".format(TEST_FILES_PATH)
        analysis_result = \
            cpdb_degs_analysis_method.call(
                cpdb_file_path=DOWNLOADED_DB_FILE,
                meta_file_path=META_FILE_PATH,
                counts_file_path=COUNTS_FILE_PATH,
                degs_file_path=DEGS_FILE_PATH,
                counts_data='ensembl',
                microenvs_file_path=MICROENVS_FILE_PATH,
                active_tfs_file_path=None,
                # Ready for v5.0.0 data release: active_tfs_file_path=ACTIVE_TFS_PATH,
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
        # Ready for v5.0.0 data release: self.assertFalse(analysis_result['CellSign_active_interactions'].empty,
        # 'CellSign_active_interactions dataframe is empty')
        # Ready for v5.0.0 data release: self.assertFalse(analysis_result['CellSign_active_interactions_deconvoluted'].empty,
        # 'CellSign_active_interactions_deconvoluted dataframe is empty')


if __name__ == "__main__":
    unittest.main()
