import pandas as pd

from cellphonedb.src.core.methods import cpdb_analysis_complex_method


def call(meta: pd.DataFrame,
         counts: pd.DataFrame,
         counts_data: str,
         interactions: pd.DataFrame,
         genes: pd.DataFrame,
         complexes: pd.DataFrame,
         complex_compositions: pd.DataFrame,
         microenvs: pd.DataFrame,
         separator: str,
         threshold: float = 0.1,
         result_precision: int = 3) -> (pd.DataFrame, pd.DataFrame, pd.DataFrame):
    means, significant_means, deconvoluted = \
        cpdb_analysis_complex_method.call(meta.copy(),
                                          counts.copy(),
                                          counts_data,
                                          interactions.copy(),
                                          genes,
                                          complexes,
                                          complex_compositions,
                                          microenvs,
                                          separator,
                                          threshold,
                                          result_precision)
    max_rank = significant_means['rank'].max()
    significant_means['rank'] = significant_means['rank'].apply(lambda rank: rank if rank != 0 else (1 + max_rank))
    significant_means.sort_values('rank', inplace=True)

    return means, significant_means, deconvoluted
