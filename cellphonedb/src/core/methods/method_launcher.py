from typing import Tuple
import pandas as pd
import numpy as np

from cellphonedb.src.core.core_logger import core_logger
from cellphonedb.src.core.database import DatabaseManager
from cellphonedb.src.core.exceptions.ThresholdValueException import ThresholdValueException
from cellphonedb.src.core.methods import cpdb_analysis_method, cpdb_statistical_analysis_method, cpdb_degs_analysis_method
from cellphonedb.src.core.preprocessors import method_preprocessors
from cellphonedb.src.core.utils.subsampler import Subsampler
from cellphonedb.src.exceptions.ParseCountsException import ParseCountsException


class MethodLauncher:
    def __init__(self, database_manager: DatabaseManager, default_threads: int, separator: str = '|'):
        self.database_manager = database_manager
        self.default_threads = default_threads
        self.separator = separator

    def __getattribute__(self, name):
        method = object.__getattribute__(self, name)
        if hasattr(method, '__call__'):
            core_logger.info('Launching Method {}'.format(name))

        return method

    def get_multidatas_from_string(self, string: str) -> pd.DataFrame:
        multidatas = self.database_manager.get_repository('multidata').get_multidatas_from_string(string)
        return multidatas

    def get_interactions_genes_complex(self) -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
        """
        Get basic data from the database for all analyses.

        Rerieves interactions, genes and complex data from the database required
        for all analyses. This data is integrated and includes 'multidata'.

        Returns
        -------
        Tuple: A tuple containing:
            - interactions: interactions joined with multidata from the database
            - genes: genes from the database
            - complex_composition: complex_compositions from the database
            - complex_expanded: complex joined with multidata from the database

        """
        # get data form database 
        interactions = self.database_manager.get_repository('interaction').get_all_expanded(include_gene=False)
        genes = self.database_manager.get_repository('gene').get_all_expanded()
        complex_composition = self.database_manager.get_repository('complex').get_all_compositions()
        complex_expanded = self.database_manager.get_repository('complex').get_all_expanded()

        # index interactions and complex data frames
        interactions.set_index('id_interaction', drop=True, inplace=True)
        complex_composition.set_index('id_complex_composition', inplace=True, drop=True)

        return interactions, genes, complex_composition, complex_expanded

    def cpdb_statistical_analysis_launcher(self,
                                           raw_meta: pd.DataFrame,
                                           counts: pd.DataFrame,
                                           counts_data: str,
                                           microenvs: pd.DataFrame,
                                           iterations: int,
                                           threshold: float,
                                           threads: int,
                                           debug_seed: int,
                                           result_precision: int,
                                           pvalue: float,
                                           subsampler: Subsampler = None,
                                           debug: bool = False,
                                           output_path: str = ''
                                           ) -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
        """
        Invokes Statistical analysis with the arguments obtained from the command line.
        """

        if threads < 1:
            core_logger.info('Using Default thread number: %s' % self.default_threads)
            threads = self.default_threads

        if threshold < 0 or threshold > 1:
            raise ThresholdValueException(threshold)

        meta = method_preprocessors.meta_preprocessor(raw_meta)
        counts = self._counts_validations(counts, meta)

        if subsampler is not None:
            counts = subsampler.subsample(counts)
            meta = meta.filter(items=(list(counts)), axis=0)

        interactions, genes, complex_composition, complex_expanded = self.get_interactions_genes_complex()

        deconvoluted, means, pvalues, significant_means = \
            cpdb_statistical_analysis_method.call(meta,
                                                  counts,
                                                  counts_data,
                                                  interactions,
                                                  genes,
                                                  complex_expanded,
                                                  complex_composition,
                                                  microenvs,
                                                  iterations,
                                                  threshold,
                                                  threads,
                                                  debug_seed,
                                                  result_precision,
                                                  pvalue,
                                                  self.separator,
                                                  debug,
                                                  output_path)

        return pvalues, means, significant_means, deconvoluted

    def cpdb_method_analysis_launcher(self,
                                      raw_meta: pd.DataFrame,
                                      counts: pd.DataFrame,
                                      counts_data: str,
                                      microenvs: pd.DataFrame,
                                      threshold: float,
                                      result_precision: int,
                                      subsampler: Subsampler = None,
                                      debug: bool =False,
                                      output_path: str = ''
                                      ) -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
        """
        Invokes Non-statistical analysis with the arguments obtained from the command line.
        """

        if threshold < 0 or threshold > 1:
            raise ThresholdValueException(threshold)
        meta = method_preprocessors.meta_preprocessor(raw_meta)

        counts = self._counts_validations(counts, meta)

        if subsampler is not None:
            counts = subsampler.subsample(counts)
            meta = meta.filter(items=list(counts), axis=0)

        interactions, genes, complex_composition, complex_expanded = self.get_interactions_genes_complex()

        means, significant_means, deconvoluted = cpdb_analysis_method.call(
            meta,
            counts,
            counts_data,
            interactions,
            genes,
            complex_expanded,
            complex_composition,
            microenvs,
            self.separator,
            threshold,
            result_precision,
            debug,
            output_path)

        return means, significant_means, deconvoluted


    def cpdb_degs_analysis_launcher(self,
                                    raw_meta: pd.DataFrame,
                                    counts: pd.DataFrame,
                                    degs: pd.DataFrame,
                                    counts_data: str,
                                    microenvs: pd.DataFrame,
                                    iterations: int,
                                    threshold: float,
                                    threads: int,
                                    debug_seed: int,
                                    result_precision: int,
                                    subsampler: Subsampler = None,
                                    debug: bool =False,
                                    output_path: str = "",
                                    ) -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
        """
        Invokes DEGs analysis with the arguments obtained from the command line.
        """

        if threads < 1:
            core_logger.info('Using Default thread number: %s' % self.default_threads)
            threads = self.default_threads

        if threshold < 0 or threshold > 1:
            raise ThresholdValueException(threshold)

        meta = method_preprocessors.meta_preprocessor(raw_meta)
        counts = self._counts_validations(counts, meta)

        if subsampler is not None:
            counts = subsampler.subsample(counts)
            meta = meta.filter(items=(list(counts)), axis=0)

        interactions, genes, complex_composition, complex_expanded = self.get_interactions_genes_complex()

        deconvoluted, means, relevant_interactions, significant_means = \
            cpdb_degs_analysis_method.call(meta,
                                            counts,
                                            degs,
                                            counts_data,
                                            interactions,
                                            genes,
                                            complex_expanded,
                                            complex_composition,
                                            microenvs,
                                            self.separator,
                                            iterations,
                                            threshold,
                                            threads,
                                            debug_seed,
                                            result_precision,
                                            debug,
                                            output_path)

        return relevant_interactions, means, significant_means, deconvoluted

    @staticmethod
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
