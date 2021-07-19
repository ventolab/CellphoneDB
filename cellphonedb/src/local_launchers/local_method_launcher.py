import os
from typing import Optional, Tuple

import numpy as np
import pandas as pd

from cellphonedb.src.app.app_logger import app_logger
from cellphonedb.src.app.cellphonedb_app import output_dir
from cellphonedb.src.core.utils.subsampler import Subsampler
from cellphonedb.utils import utils
from cellphonedb.utils.utils import write_to_file


class LocalMethodLauncher(object):
    def __init__(self, cellphonedb_app):

        self.cellphonedb_app = cellphonedb_app

    def __getattribute__(self, name):
        method = object.__getattribute__(self, name)
        if hasattr(method, '__call__'):
            app_logger.info('Launching Method {}'.format(name))

        return method

    def cpdb_statistical_analysis_local_method_launcher(self, meta_filename: str,
                                                        counts_filename: str,
                                                        counts_data: str,
                                                        microenvs_filename: str = '',
                                                        project_name: str = '',
                                                        iterations: int = 1000,
                                                        threshold: float = 0.1,
                                                        output_path: str = '',
                                                        output_format: Optional[str] = None,
                                                        means_filename: str = 'means',
                                                        pvalues_filename: str = 'pvalues',
                                                        significant_means_filename: str = 'significant_means',
                                                        deconvoluted_filename='deconvoluted',
                                                        debug_seed: int = -1,
                                                        threads: int = -1,
                                                        result_precision: int = 3,
                                                        pvalue: float = 0.05,
                                                        subsampler: Subsampler = None,
                                                        debug: bool = False
                                                        ) -> None:
        output_path = self._set_paths(output_path, project_name)

        debug_seed = int(debug_seed)
        iterations = int(iterations)
        threads = int(threads)
        threshold = float(threshold)
        result_precision = int(result_precision)

        counts, meta = self._load_meta_counts(counts_filename, meta_filename)
        self._check_counts_data(counts, counts_data)
        if not microenvs_filename:
            microenvs = pd.DataFrame()
        else:
            microenvs = self._load_microenvs(microenvs_filename, meta)
        
        pvalues_simple, means_simple, significant_means_simple, deconvoluted_simple = \
            self.cellphonedb_app.method.cpdb_statistical_analysis_launcher(
                meta,
                counts,
                counts_data,
                microenvs,
                iterations,
                threshold,
                threads,
                debug_seed,
                result_precision,
                pvalue,
                subsampler,
                debug,
                output_path
            )

        write_to_file(means_simple, means_filename, output_path, output_format)
        write_to_file(pvalues_simple, pvalues_filename, output_path, output_format)
        write_to_file(significant_means_simple, significant_means_filename, output_path, output_format)
        write_to_file(deconvoluted_simple, deconvoluted_filename, output_path, output_format)

    def cpdb_analysis_local_method_launcher(self, meta_filename: str,
                                            counts_filename: str,
                                            counts_data: str,
                                            microenvs_filename: str = '',
                                            project_name: str = '',
                                            threshold: float = 0.1,
                                            output_path: str = '',
                                            output_format: Optional[str] = None,
                                            means_filename: str = 'means',
                                            significant_means_filename: str = 'significant_means',
                                            deconvoluted_filename='deconvoluted',
                                            result_precision: int = 3,
                                            subsampler: Subsampler = None,
                                            debug: bool = False,
                                            ) -> None:
        output_path = self._set_paths(output_path, project_name)

        result_precision = int(result_precision)
        threshold = float(threshold)

        counts, meta = self._load_meta_counts(counts_filename, meta_filename)
        self._check_counts_data(counts, counts_data)
        if not microenvs_filename:
            microenvs = pd.DataFrame()
        else:
            microenvs = self._load_microenvs(microenvs_filename, meta)

        means, significant_means, deconvoluted = \
            self.cellphonedb_app.method.cpdb_method_analysis_launcher(meta,
                                                                      counts,
                                                                      counts_data,
                                                                      microenvs,
                                                                      threshold,
                                                                      result_precision,
                                                                      subsampler,
                                                                      debug,
                                                                      output_path)

        write_to_file(means, means_filename, output_path, output_format)
        write_to_file(significant_means, significant_means_filename, output_path, output_format)
        write_to_file(deconvoluted, deconvoluted_filename, output_path, output_format)


    def cpdb_degs_analysis_local_method_launcher(self, meta_filename: str,
                                                        counts_filename: str,
                                                        degs_filename: str,
                                                        counts_data: str,
                                                        microenvs_filename: str,
                                                        project_name: str = '',
                                                        iterations: int = 1000,
                                                        threshold: float = 0.1,
                                                        output_path: str = '',
                                                        output_format: Optional[str] = None,
                                                        means_filename: str = 'means',
                                                        relevant_interactions_filename: str = 'relevant_interactions',
                                                        significant_means_filename: str = 'significant_means',
                                                        deconvoluted_filename='deconvoluted',
                                                        debug_seed: int = -1,
                                                        threads: int = -1,
                                                        result_precision: int = 3,
                                                        subsampler: Subsampler = None,
                                                        debug: bool = True,
                                                        ) -> None:
        output_path = self._set_paths(output_path, project_name)

        debug_seed = int(debug_seed)
        iterations = int(iterations)
        threads = int(threads)
        threshold = float(threshold)
        result_precision = int(result_precision)

        counts, meta = self._load_meta_counts(counts_filename, meta_filename)
        self._check_counts_data(counts, counts_data)
        degs = self._load_degs(degs_filename, meta)
        if not microenvs_filename:
            microenvs = pd.DataFrame()
        else:
            microenvs = self._load_microenvs(microenvs_filename, meta)
        
        relevant_interactions, means_simple, significant_means_simple, deconvoluted_simple = \
            self.cellphonedb_app.method.cpdb_degs_analysis_launcher(
                meta,
                counts,
                degs,
                counts_data,
                microenvs,
                iterations,
                threshold,
                threads,
                debug_seed,
                result_precision,
                subsampler,
                debug,
                output_path
            )

        write_to_file(means_simple, means_filename, output_path, output_format)
        write_to_file(relevant_interactions, relevant_interactions_filename, output_path, output_format)
        write_to_file(significant_means_simple, significant_means_filename, output_path, output_format)
        write_to_file(deconvoluted_simple, deconvoluted_filename, output_path, output_format)


    @staticmethod
    def _path_is_empty(path):
        return bool([f for f in os.listdir(path) if not f.startswith('.')])

    @staticmethod
    def _set_paths(output_path, project_name):
        """
        Set and create the output path.
        
        Defines output path based on the output_path and project_name options
        and then creates the required folder if it doesn't exist already.
        """
        if not output_path:
            output_path = output_dir
        if project_name:
            output_path = os.path.realpath(os.path.expanduser('{}/{}'.format(output_path, project_name)))
        os.makedirs(output_path, exist_ok=True)
        if LocalMethodLauncher._path_is_empty(output_path):
            app_logger.warning(
                'Output directory ({}) exist and is not empty. Result can overwrite old results'.format(output_path))
        return output_path

    @staticmethod
    def _load_meta_counts(counts_filename: str, meta_filename: str) -> Tuple[pd.DataFrame, pd.DataFrame]:
        """Load meta and counts file

        This methdos reads both meta and counts input files required
        for all methods. Also runs validations and cleanup on both files.

        Parameters
        ----------
        counts_filename: str
            Path to the counts file.
        meta_filename: str
            Path to the meta file.

        Returns
        -------
        Tuple: A tuple containing:
            - counts loaded as a DataFrame.
            - meta loaded as a DataFrame.
        
        Raises
        ------
        ParseMetaException
            Error when parsing Meta file.
        """
        meta = utils.read_data_table_from_file(os.path.realpath(meta_filename))
        counts = utils.read_data_table_from_file(os.path.realpath(counts_filename), index_column_first=True)

        return counts, meta


    @staticmethod
    def _load_microenvs(microenvs_filename: str, meta: pd.DataFrame) -> pd.DataFrame:
        """Load microenvironment file

        This methdos reads a microenvironment file into a DataFrame. 
        Runs validations to make sure the file has enough columns and
        that all the clusters in the microenvironment are included in meta.

        Parameters
        ----------
        microenvs_filename
            Path to the microenvironments file.
        meta
            Meta DataFrame.

        Returns
        -------
        pd.DataFrame
            Microenvrionments as a DataFrame read fro the input file.
        """
        CELL_TYPE = "cell_type"
        MICRO_ENVIRONMENT = "microenvironment"
        microenvs = utils.read_data_table_from_file(os.path.realpath(microenvs_filename))
        microenvs.drop_duplicates(inplace = True)
        len_columns = len(microenvs.columns)
        if len_columns<2:
            raise Exception(f"Missing columns in microenvironments: 2 required but {len_columns} provieded")
        elif len_columns>2:
            app_logger.warn(f"Microenvrionemnts expects 2 columns and got {len_columns}. Droppoing extra columns.")
        microenvs = microenvs.iloc[:, 0:2]
        if any(~microenvs.iloc[:,0].isin(meta.iloc[:,1])):
            raise Exception("Some clusters/cell_types in microenvironments are not present in meta")
        microenvs.columns = [CELL_TYPE,MICRO_ENVIRONMENT]
        return microenvs


    @staticmethod
    def _load_degs(degs_filename: str, meta: pd.DataFrame) -> pd.DataFrame:
        """Load DEGs file

        This methdos reads a DEGs file into a DataFrame. Runs validations
        to make sure the file has enough columns and that all the clusters
        in DEGs are included in meta.

        Parameters
        ----------
        degs_filename
            Path to the DEGs file.
        meta
            Meta DataFrame.

        Returns
        -------
        DataFrame
            DEGs as a DataFrame with cluster and genes columns.
        """
        CLUSTER = "cluster"
        GENE = "gene"
        degs_filename = os.path.realpath(degs_filename)
        degs = utils.read_data_table_from_file(degs_filename)
        len_columns = len(degs.columns)
        if len_columns<2:
            raise Exception(f"Missing columns in DEGs: 2 required but {len_columns} provieded")
        #elif len_columns>2:
        #    app_logger.warn(f"DEGs expects 2 columns and got {len_columns}. Droppoing extra columns.")
        degs = degs.iloc[:, 0:2]

        if any(~degs.iloc[:,0].isin(meta.iloc[:,1])):
            raise Exception("Some clusters/cell_types in DEGs are not present in meta")
        degs.columns = [CLUSTER,GENE]
        degs.drop_duplicates(inplace=True)
        return degs

    @staticmethod
    def _check_counts_data(counts: pd.DataFrame, counts_data: str) -> None:
        """Naive check count data agains counts gene names. 

        This methdos quickly checks if count_data matches the all gene names and
        gives a comprehensive warning.

        Parameters
        ----------
        counts: pd.DataFrame
            Counts data
        counts_data: str
            Gene format expected in counts data
        """
        if ~np.all(counts.index.str.startswith("ENSG0")) and counts_data=="ensembl":
            app_logger.warn(f"Gene format missmatch. Using gene type '{counts_data}' expects gene names to start with 'ENSG' but some genes seem to be in another format. "
            "Try using '--counts-data hgnc_symbol' if all counts are filtered.")
