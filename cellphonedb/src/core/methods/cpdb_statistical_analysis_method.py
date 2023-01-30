from typing import Tuple
import pandas as pd

from cellphonedb.src.core.methods import cpdb_statistical_analysis_complex_method
from cellphonedb.utils import db_utils, file_utils

def call(cpdb_file_path: str,
         meta: pd.DataFrame,
         count: pd.DataFrame,
         counts_data: str,
         output_path: str,
         microenvs: pd.DataFrame,
         iterations: int,
         threshold: float,
         threads: int,
         debug_seed: int,
         result_precision: int,
         pvalue: float,
         separator: str = '|',
         debug: bool = False,
         output_suffix: str = None
         ) -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Statistical method for analysis

     This methods calculates the mean and percent for the cluster interactions
     and for each gene interaction. No shuffling nor DEGs are involved.

     Parameters
     ----------
    cpdb_file_path: str
        CellphoneDB database file path
     meta: str
         Meta data.
     counts: str
         Counts data.
     counts_data: str
         Type of gene identifiers in the counts data: "ensembl", "gene_name", "hgnc_symbol"
     output_path: str
        Output path used to store the analysis results (and to store intermediate files when debugging)
     microenvs: pd.DataFrame
         Micro-environment data to limit cluster interactions
     iterations: int
        Number of times cell type labels will be shuffled across cells in order to
        determine statistically significant ligand/receptor expression means.
     threshold: float
         Percentage of cells expressing the specific ligand/receptor [0.0 - 1.0]
     threads: int
        Number of threads to be used during the shuffling of clusters/cell types across cells
     debug_seed: int
        This parameter is used for testing only (and only in single-threaded mode
        only - see: https://stackoverflow.com/questions/21494489/what-does-numpy-random-seed0-do).
     result_precision: int
         Number of decimal digits in results.
     pvalue: float
         A p-value below which a ligand/receptor expression mean is considered to be
         statistically significant.
     separator: str
         Separator for pairs of genes (gene1|gene2) and clusters (cluster1|cluster2).
     debug: bool
         Storge intermediate data as pickle file (debug_intermediate.pkl).
     output_suffix: str, optional
         Suffix to append to the result file's name (if not provided, timestamp will be used)

     Returns
     -------
     Tuple
         - pvalues
         - means_result
         - significant_means
         - deconvoluted_result
     """

    # Load into memory CellphoneDB data
    interactions, genes, complex_composition, complex_expanded = \
        db_utils.get_interactions_genes_complex(cpdb_file_path)
    
    pvalues, means, significant_means, deconvoluted = \
        cpdb_statistical_analysis_complex_method.call(meta.copy(),
                                                      count.copy(),
                                                      counts_data,
                                                      interactions,
                                                      genes,
                                                      complex_expanded,
                                                      complex_composition,
                                                      microenvs,
                                                      pvalue,
                                                      separator,
                                                      iterations,
                                                      threshold,
                                                      threads,
                                                      debug_seed,
                                                      result_precision,
                                                      debug,
                                                      output_path
                                                      )


    max_rank = significant_means['rank'].max()
    significant_means['rank'] = significant_means['rank'].apply(lambda rank: rank if rank != 0 else (1 + max_rank))
    significant_means.sort_values('rank', inplace=True)

    file_utils.save_dfs_as_csv(output_path, output_suffix, "statistical_analysis", \
                            {"deconvoluted" : deconvoluted, \
                            "means" : means, \
                            "pvalues" : pvalues, \
                            "significant_means" : significant_means} )

    return deconvoluted, means, pvalues, significant_means
