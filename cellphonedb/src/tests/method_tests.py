import unittest
import os
import fnmatch
import io
import re
import shutil
import tempfile
import zipfile
from contextlib import redirect_stdout
import pandas as pd
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


class CreateDbUnitTests(unittest.TestCase):
    """Unit tests for db_utils.create_db on hand-written *_input.csv files (no database download)."""

    GENE_INPUT = ("gene_name,uniprot,hgnc_symbol,ensembl\n"
                  "LIG1,P00001,LIG1,ENSG00000000001\n"
                  "SUB1,P00003,SUB1,ENSG00000000003\n"
                  "SUB2,P00004,SUB2,ENSG00000000004\n")
    PROTEIN_INPUT = ("uniprot,protein_name,transmembrane,peripheral,secreted,secreted_desc,secreted_highlight,"
                     "receptor,receptor_desc,integrin,other,other_desc,tags,tags_reason,tags_description\n"
                     "P00001,LIG1_HUMAN,False,False,True,,True,False,,False,False,,,,\n"
                     "P00003,SUB1_HUMAN,True,False,False,,False,True,,False,False,,,,\n"
                     "P00004,SUB2_HUMAN,True,False,False,,False,True,,False,False,,,,\n")
    # The released complex_input.csv header, but the only complex is a dimer: uniprot_3 and uniprot_4
    # are empty in every row, so pandas reads them as all-NaN float64 columns (see issue #224).
    COMPLEX_INPUT = ("complex_name,uniprot_1,uniprot_2,uniprot_3,uniprot_4,transmembrane,peripheral,secreted,"
                     "secreted_desc,secreted_highlight,receptor,receptor_desc,integrin,other,other_desc,"
                     "pdb_id,pdb_structure,stoichiometry,comments_complex\n"
                     "RECCPLX,P00003,P00004,,,True,False,False,,False,True,,False,False,,,,,\n")
    INTERACTION_INPUT = ("partner_a,partner_b,protein_name_a,protein_name_b,annotation_strategy,source\n"
                         "P00001,RECCPLX,LIG1_HUMAN,,curated,test\n")

    def setUp(self):
        self.input_dir = tempfile.mkdtemp()
        self.addCleanup(shutil.rmtree, self.input_dir, ignore_errors=True)
        for file_name, content in [("gene_input.csv", self.GENE_INPUT), ("protein_input.csv", self.PROTEIN_INPUT),
                                   ("complex_input.csv", self.COMPLEX_INPUT),
                                   ("interaction_input.csv", self.INTERACTION_INPUT)]:
            with open(os.path.join(self.input_dir, file_name), "w") as f:
                f.write(content)

    def test_create_db_with_unused_subunit_columns(self):
        """complex_input.csv with entirely empty uniprot_N columns must build (used to raise ValueError)."""
        db_utils.create_db(self.input_dir)
        db_files = [f for f in os.listdir(self.input_dir) if fnmatch.fnmatch(f, GENERATED_CPDB_PATTERN)]
        assert len(db_files) == 1, "expected one generated database file, got {}".format(db_files)
        with zipfile.ZipFile(os.path.join(self.input_dir, db_files[0]), 'r') as zip_ref:
            multidata = pd.read_csv(zip_ref.open('multidata_table.csv'))
            complex_composition = pd.read_csv(zip_ref.open('complex_composition_table.csv'))
        assert multidata[multidata['is_complex']]['name'].tolist() == ['RECCPLX']
        assert len(complex_composition) == 2
        assert complex_composition['total_protein'].tolist() == [2, 2]

    def test_unknown_complex_proteins_still_reported(self):
        """The sanity test must still name subunits missing from protein_input.csv, and skip empty slots."""
        protein_df = pd.read_csv(io.StringIO(self.PROTEIN_INPUT))
        complex_df = pd.read_csv(io.StringIO(
            self.COMPLEX_INPUT + "OTHERCPLX,P00003,P99999,,,True,False,False,,False,True,,False,False,,,,,\n"))
        assert str(complex_df['uniprot_4'].dtype) == 'float64'  # the all-NaN column that used to break the merge
        out = io.StringIO()
        with redirect_stdout(out):
            db_utils.sanity_test_report_unknown_proteins(
                protein_df, complex_df, ['uniprot_1', 'uniprot_2', 'uniprot_3', 'uniprot_4'])
        assert 'P99999' in out.getvalue()
        assert 'P00003' not in out.getvalue() and 'nan' not in out.getvalue()


if __name__ == "__main__":
    unittest.main()
