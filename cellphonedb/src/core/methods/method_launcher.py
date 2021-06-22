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
                                           ) -> (pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame):

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

        interactions = self.database_manager.get_repository('interaction').get_all_expanded(include_gene=False)
        genes = self.database_manager.get_repository('gene').get_all_expanded()
        complex_composition = self.database_manager.get_repository('complex').get_all_compositions()
        complex_expanded = self.database_manager.get_repository('complex').get_all_expanded()

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
                                                  self.separator)

        return pvalues, means, significant_means, deconvoluted

    def cpdb_method_analysis_launcher(self,
                                      raw_meta: pd.DataFrame,
                                      counts: pd.DataFrame,
                                      counts_data: str,
                                      microenvs: pd.DataFrame,
                                      threshold: float,
                                      result_precision: int,
                                      subsampler: Subsampler = None,
                                      ) -> (pd.DataFrame, pd.DataFrame, pd.DataFrame):

        if threshold < 0 or threshold > 1:
            raise ThresholdValueException(threshold)
        meta = method_preprocessors.meta_preprocessor(raw_meta)

        counts = self._counts_validations(counts, meta)

        if subsampler is not None:
            counts = subsampler.subsample(counts)
            meta = meta.filter(items=list(counts), axis=0)

        interactions = self.database_manager.get_repository('interaction').get_all_expanded(include_gene=False)
        genes = self.database_manager.get_repository('gene').get_all_expanded()
        complex_composition = self.database_manager.get_repository('complex').get_all_compositions()
        complex_expanded = self.database_manager.get_repository('complex').get_all_expanded()

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
            result_precision)

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
                                    pvalue: float,
                                    subsampler: Subsampler = None,
                                    ) -> (pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame):

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

        interactions = self.database_manager.get_repository('interaction').get_all_expanded(include_gene=False)
        genes = self.database_manager.get_repository('gene').get_all_expanded()
        complex_composition = self.database_manager.get_repository('complex').get_all_compositions()
        complex_expanded = self.database_manager.get_repository('complex').get_all_expanded()

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
                                            iterations,
                                            threshold,
                                            threads,
                                            debug_seed,
                                            result_precision,
                                            pvalue,
                                            self.separator)

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
            raise ParseCountsException('Some cells in meta didnt exist in counts',
                                       'Maybe incorrect file format')
        if np.any(~counts.columns.isin(meta.index)):
            counts = counts.loc[:, counts.columns.isin(meta.index)].copy()
        return counts
