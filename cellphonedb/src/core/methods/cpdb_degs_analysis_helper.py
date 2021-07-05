from typing import Tuple
import warnings

warnings.simplefilter("ignore", UserWarning)

import pandas as pd
import numpy as np
import numpy_groupies as npg
from cellphonedb.src.core.core_logger import core_logger
from cellphonedb.src.core.models.complex import complex_helper


def get_significant_means(real_mean_analysis: pd.DataFrame,
                          relevant_interactions: pd.DataFrame) -> pd.DataFrame:
    """
    If the result_percent value is > min_significant_mean, sets the value to NaN, else, sets the mean value.

    EXAMPLE:
    INPUT:

    real mean
                cluster1    cluster2    cluster
    ensembl1    0.1         1.0         2.0
    ensembl2    2.0         0.1         0.2
    ensembl3    0.3         0.0         0.5

    result percent

                cluster1    cluster2    cluster
    ensembl1    0.0         1.0         1.0
    ensembl2    0.04        0.03        0.62
    ensembl3    0.3         0.55        0.02

    min_significant_mean = 0.05

    RESULT:

                cluster1    cluster2    cluster
    ensembl1    0.1         NaN         NaN
    ensembl2    2.0         0.1         NaN
    ensembl3    NaN         NaN         0.5
    """
    significant_means = real_mean_analysis.values.copy()
    mask = relevant_interactions == 0
    significant_means[mask] = np.nan
    return pd.DataFrame(significant_means,
                        index=real_mean_analysis.index,
                        columns=real_mean_analysis.columns)


def build_significant_means(real_mean_analysis: pd.DataFrame,
                            relevant_interactions: pd.DataFrame) -> Tuple[pd.Series, pd.DataFrame]:
    """
    Calculates the significant means and add rank (number of non empty entries divided by total entries)
    """
    significant_means = get_significant_means(real_mean_analysis, relevant_interactions)
    significant_mean_rank = significant_means.count(axis=1)  # type: pd.Series
    number_of_clusters = len(significant_means.columns)
    significant_mean_rank = significant_mean_rank.apply(lambda rank: rank / number_of_clusters)
    significant_mean_rank = significant_mean_rank.round(3)
    significant_mean_rank.name = 'rank'
    return significant_mean_rank, significant_means


def degs_analysis(degs: pd.DataFrame,
                genes: pd.DataFrame,
                interactions: pd.DataFrame,
                cluster_interactions: list,
                counts_data:pd.DataFrame,
                separator: str) -> pd.DataFrame:
    """
    Calculates the filter matrix based on the DEGs for the list of interactions and for each cluster
    Sets 1 if one of both is 1.

    EXAMPLE:
        DEGs (input deg file is transformed into matrix)
        
               cluster1   cluster2    cluster3
        ensembl1     0           1           0
        ensembl2     1           1           0
        ensembl3     0           0           1

        interactions:

        ensembl1,ensembl2
        ensembl2,ensembl3

        RESULT:
                            cluster1_cluster1   cluster1_cluster2   ...   cluster3_cluster2   cluster3_cluster3
        ensembl1_ensembl2      (0,1) -> 1          (0,1) -> 1                (0,1) -> 1             (0,0) -> 0
        ensembl2_ensembl3      (1,0) -> 1          (1,0) -> 1                (0,0) -> 0             (0,1) -> 1
        
    """    

    #Prepare DEGs matrix from input file
    degs["deg"]=1
    degs = pd.pivot_table(degs, values="deg", index="gene", columns="cluster", fill_value=0)
    degs = degs.merge(genes[['id_multidata', 'ensembl', 'gene_name', 'hgnc_symbol']], left_on="gene", right_on=counts_data)
    degs.set_index("id_multidata", inplace=True)
    degs = degs[~degs.index.duplicated(keep='first')]

    GENE_ID1 = 'multidata_1_id'
    GENE_ID2 = 'multidata_2_id'

    cluster1_names = cluster_interactions[:, 0]
    cluster2_names = cluster_interactions[:, 1]
    gene1_ids = interactions[GENE_ID1].values
    gene2_ids = interactions[GENE_ID2].values

    x = degs.reindex(index=gene1_ids, columns=cluster1_names, fill_value=0).values
    y = degs.reindex(index=gene2_ids, columns=cluster2_names, fill_value=0).values

    result = pd.DataFrame(
        ((x + y) > 0).astype(int),
        index=interactions.index,
        columns=(pd.Series(cluster1_names) + separator + pd.Series(cluster2_names)).values)

    return result