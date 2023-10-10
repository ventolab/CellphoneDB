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
targetDir = os.path.join(script_dir, "test_data")
GENERATED_CPDB_PATTERN="cellphonedb_*.zip"
downloaded_db_dir = os.path.join(targetDir, "downloaded_db")
generated_db_dir = os.path.join(targetDir, "generated_db")
CELLPHONEDB_FILES = ['complex_composition_table.csv','gene_synonym_to_gene_name.csv','interaction_table.csv','complex_table.csv','gene_table.csv','multidata_table.csv']

class UnitTests(unittest.TestCase):

    def setUp(self):
        if not os.path.isdir(targetDir):
            os.mkdir(targetDir)
        if not os.path.isdir(downloaded_db_dir):
            os.mkdir(downloaded_db_dir)
        if not os.path.isdir(generated_db_dir):
            os.mkdir(generated_db_dir)

    def test_download_database(self):
        db_utils.download_database(targetDir, RELEASED_VERSION)

    def test_createDatabase(self):
        for f in os.listdir(targetDir):
            if re.search("cellphonedb_.*\\.zip", f):
                os.remove(os.path.join(targetDir, f))
        db_utils.create_db(targetDir)
        generated_db_file = None
        for f in os.listdir(targetDir):
            if fnmatch.fnmatch(f, GENERATED_CPDB_PATTERN):
                generated_db_file = f
                break
        assert generated_db_file is not None

    def test_compare_downloaded_created_dbs(self):
        generated_db_file = None
        for f in os.listdir(targetDir):
            if fnmatch.fnmatch(f, GENERATED_CPDB_PATTERN):
                generated_db_file = f
                break
        with zipfile.ZipFile(os.path.join(targetDir, generated_db_file), 'r') as zip_ref:
            zip_ref.extractall(generated_db_dir)
        with zipfile.ZipFile(os.path.join(targetDir, "cellphonedb.zip"), 'r') as zip_ref:
            zip_ref.extractall(downloaded_db_dir)
        for f in CELLPHONEDB_FILES:
            generated_db_f = os.path.join(targetDir, generated_db_dir, f)
            downloaded_db_f = os.path.join(targetDir, downloaded_db_dir, f)
            with open(generated_db_f, 'r') as fp:
                for generated_count, line in enumerate(fp):
                    pass
            with open(downloaded_db_f, 'r') as fp:
                for downloaded_count, line in enumerate(fp):
                    pass
            print("Comparing {} between generated and downloaded DB".format(f))
            assert generated_count == downloaded_count, "The number of lines in {} differs between generated and downloaded DB".format(f)

if __name__ == "__main__":
    unittest.main()