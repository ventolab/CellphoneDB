import warnings
warnings.simplefilter("ignore", UserWarning)

import itertools
from functools import partial
from multiprocessing.pool import Pool

import pandas as pd
import numpy as np
import numpy_groupies as npg
from cellphonedb.src.core.core_logger import core_logger
from cellphonedb.src.core.models.complex import complex_helper


def get_significant_means(real_mean_analysis: pd.DataFrame,
                          result_percent: pd.DataFrame,
                          min_significant_mean: float) -> pd.DataFrame:
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
    mask = result_percent > min_significant_mean
    significant_means[mask] = np.nan
    return pd.DataFrame(significant_means,
                        index=real_mean_analysis.index,
                        columns=real_mean_analysis.columns)


def shuffle_meta(meta: pd.DataFrame) -> pd.DataFrame:
    """
    Permutates the meta values aleatory generating a new meta file
    """
    meta_copy = meta.copy()
    np.random.shuffle(meta_copy['cell_type'].values)

    return meta_copy


def build_clusters(meta: pd.DataFrame,
                   counts: pd.DataFrame,
                   complex_composition: pd.DataFrame,
                   skip_percent: bool) -> dict:
    """
    Builds a cluster structure and calculates the means values
    """
    CELL_TYPE = 'cell_type'
    COMPLEX_ID = 'complex_multidata_id'
    PROTEIN_ID = 'protein_multidata_id'

    meta[CELL_TYPE] = meta[CELL_TYPE].astype('category')
    cluster_names = meta[CELL_TYPE].cat.categories

    # Simple genes cluster counts
    cluster_means = pd.DataFrame(
        npg.aggregate(meta[CELL_TYPE].cat.codes, counts.values, func='mean', axis=1),
        index=counts.index,
        columns=cluster_names.to_list()
    )
    if not skip_percent:
        cluster_pcts = pd.DataFrame(
            npg.aggregate(meta[CELL_TYPE].cat.codes, (counts > 0).astype(int).values, func='mean', axis=1),
            index=counts.index,
            columns=cluster_names.to_list()
        )
    else:
        cluster_pcts = pd.DataFrame(index=counts.index, columns=cluster_names.to_list())

    # Complex genes cluster counts
    if not complex_composition.empty:
        complexes = complex_composition.groupby(COMPLEX_ID).apply(lambda x: x[PROTEIN_ID].values).to_dict()
        complex_cluster_means = pd.DataFrame(
            {complex_id: cluster_means.loc[protein_ids].min(axis=0).values
             for complex_id, protein_ids in complexes.items()},
            index=cluster_means.columns
        ).T
        cluster_means = cluster_means.append(complex_cluster_means)
        if not skip_percent:
            complex_cluster_pcts = pd.DataFrame(
                {complex_id: cluster_pcts.loc[protein_ids].min(axis=0).values
             for complex_id, protein_ids in complexes.items()},
                index=cluster_pcts.columns
            ).T
            cluster_pcts = cluster_pcts.append(complex_cluster_pcts)

    return {'names': cluster_names, 'means': cluster_means, 'percents': cluster_pcts}


def filter_counts_by_interactions(counts: pd.DataFrame,
                                  interactions: pd.DataFrame) -> pd.DataFrame:
    """
    Removes counts if is not defined in interactions components
    """
    multidata_genes_ids = interactions['multidata_1_id'].append(
        interactions['multidata_2_id']).drop_duplicates().tolist()

    counts_filtered = counts.filter(multidata_genes_ids, axis=0)

    return counts_filtered


def filter_empty_cluster_counts(counts: pd.DataFrame) -> pd.DataFrame:
    """
    Remove count with all values to zero
    """
    if counts.empty:
        return counts

    filtered_counts = counts[counts.apply(lambda row: row.sum() > 0, axis=1)]
    return filtered_counts


def mean_pvalue_result_build(real_mean_analysis: pd.DataFrame, result_percent: pd.DataFrame,
                             interactions_data_result: pd.DataFrame) -> pd.DataFrame:
    """
    Merges the pvalues and means in one table
    """
    mean_pvalue_result = pd.DataFrame(index=real_mean_analysis.index)
    for interaction_cluster in real_mean_analysis.columns.values:
        mean_pvalue_result[interaction_cluster] = real_mean_analysis[interaction_cluster].astype(str).str.cat(
            result_percent[interaction_cluster].astype(str), sep=' | ')

    mean_pvalue_result = pd.concat([interactions_data_result, mean_pvalue_result], axis=1, join='inner', sort=False)

    return mean_pvalue_result


def get_cluster_combinations(cluster_names: np.array, microenvs: pd.DataFrame = pd.DataFrame()) -> np.array:
    """
    Calculates and sort combinations including itself

    ie

    INPUT
    cluster_names = ['cluster1', 'cluster2', 'cluster3']

    RESULT
    [('cluster1','cluster1'),('cluster1','cluster2'),('cluster1','cluster3'),
     ('cluster2','cluster1'),('cluster2','cluster2'),('cluster2','cluster3'),
     ('cluster3','cluster1'),('cluster3','cluster2'),('cluster3','cluster3')]

    if microenvironments are provided combinations are performed only within each microenv

    INPUT
    cluster_names = ['cluster1', 'cluster2', 'cluster3']
    microenvs = [
        ('cluster1', 'env1'),
        ('cluster2', 'env1'),
        ('cluster3', 'env2')]

    RESULT
    [('cluster1','cluster1'),('cluster1','cluster2'),
     ('cluster2','cluster1'),('cluster2','cluster2'),
     ('cluster3','cluster3')]

    """
    result = np.array([])
    if microenvs.empty:
        result = np.array(np.meshgrid(cluster_names.values, cluster_names.values)).T.reshape(-1, 2)
    else:
        core_logger.info('Limiting cluster combinations using microenvironments')
        cluster_combinations = []
        for me in microenvs["microenvironment"].unique():
            me_cell_types = microenvs[microenvs["microenvironment"]==me]["cell_type"]
            combinations = np.array(np.meshgrid(me_cell_types, me_cell_types))
            cluster_combinations.extend(combinations.T.reshape(-1, 2))
        result = pd.DataFrame(cluster_combinations).drop_duplicates().to_numpy()
    core_logger.debug(f'Using {len(result)} cluster combinations for analysis')
    return result

def build_result_matrix(interactions: pd.DataFrame, cluster_interactions: list, separator: str) -> pd.DataFrame:
    """
    builds an empty cluster matrix to fill it later
    """
    columns = []

    for cluster_interaction in cluster_interactions:
        columns.append('{}{}{}'.format(cluster_interaction[0], separator, cluster_interaction[1]))

    result = pd.DataFrame(index=interactions.index, columns=columns, dtype=float)

    return result


def mean_analysis(interactions: pd.DataFrame,
                  clusters: dict,
                  cluster_interactions: list,
                  separator: str) -> pd.DataFrame:
    """
    Calculates the mean for the list of interactions and for each cluster

    sets 0 if one of both is 0

    EXAMPLE:
        cluster_means
                   cluster1    cluster2    cluster3
        ensembl1     0.0         0.2         0.3
        ensembl2     0.4         0.5         0.6
        ensembl3     0.7         0.0         0.9

        interactions:

        ensembl1,ensembl2
        ensembl2,ensembl3

        RESULT:
                              cluster1_cluster1   cluster1_cluster2   ...   cluster3_cluster2   cluster3_cluster3
        ensembl1_ensembl2     mean(0.0,0.4)*      mean(0.0,0.5)*            mean(0.3,0.5)       mean(0.3,0.6)
        ensembl2_ensembl3     mean(0.4,0.7)       mean(0.4,0.0)*            mean(0.6,0.0)*      mean(0.6,0.9)


        results with * are 0 because one of both components is 0.
    """
    GENE_ID1 = 'multidata_1_id'
    GENE_ID2 = 'multidata_2_id'

    cluster1_names = cluster_interactions[:, 0]
    cluster2_names = cluster_interactions[:, 1]
    gene1_ids = interactions[GENE_ID1].values
    gene2_ids = interactions[GENE_ID2].values

    x = clusters['means'].loc[gene1_ids, cluster1_names].values
    y = clusters['means'].loc[gene2_ids, cluster2_names].values

    result = pd.DataFrame(
        (x > 0) * (y > 0) * (x + y) / 2,
        index=interactions.index,
        columns=(pd.Series(cluster1_names) + separator + pd.Series(cluster2_names)).values)

    return result


def percent_analysis(clusters: dict,
                     threshold: float,
                     interactions: pd.DataFrame,
                     cluster_interactions: list,
                     separator: str) -> pd.DataFrame:
    """
    Calculates the percents for cluster interactions and foreach gene interaction

    If one of both is not 0 sets the value to 0. Else sets 1

    EXAMPLE:
        INPUT:

        threshold = 0.1
        cluster1 = cell1,cell2
        cluster2 = cell3

                     cell1       cell2      cell3
        ensembl1     0.0         0.6         0.3
        ensembl2     0.1         0.05         0.06
        ensembl3     0.0         0.0         0.9

        interactions:

        ensembl1,ensembl2
        ensembl1,ensembl3


        (after percents calculation)

                     cluster1    cluster2
        ensembl1     0           0
        ensembl2     1           1
        ensembl3     1           0

        RESULT:
                            cluster1_cluster1   cluster1_cluster2   cluster2_cluster1   cluster2_cluster2
        ensembl1_ensembl2   (0,1)-> 0           (0,1)-> 0           (0,1)->0            (0,1)->0
        ensembl1_ensembl3   (0,1)-> 0           (0,0)-> 1           (0,1)->0            (0,0)->1


    """
    GENE_ID1 = 'multidata_1_id'
    GENE_ID2 = 'multidata_2_id'

    cluster1_names = cluster_interactions[:, 0]
    cluster2_names = cluster_interactions[:, 1]
    gene1_ids = interactions[GENE_ID1].values
    gene2_ids = interactions[GENE_ID2].values

    x = clusters['percents'].loc[gene1_ids, cluster1_names].values
    y = clusters['percents'].loc[gene2_ids, cluster2_names].values

    result = pd.DataFrame(
        ((x > threshold) * (y > threshold)).astype(int),
        index=interactions.index,
        columns=(pd.Series(cluster1_names) + separator + pd.Series(cluster2_names)).values)

    return result


def shuffled_analysis(iterations: int,
                      meta: pd.DataFrame,
                      counts: pd.DataFrame,
                      interactions: pd.DataFrame,
                      cluster_interactions: list,
                      complex_composition: pd.DataFrame,
                      real_mean_analysis: pd.DataFrame,
                      base_result: pd.DataFrame,
                      threads: int,
                      separator: str) -> list:
    """
    Shuffles meta and calculates the means for each and saves it in a list.

    Runs it in a multiple threads to run it faster
    """
    with Pool(processes=threads) as pool:
        statistical_analysis_thread = partial(_statistical_analysis,
                                              base_result,
                                              cluster_interactions,
                                              counts,
                                              interactions,
                                              meta,
                                              complex_composition,
                                              separator,
                                              real_mean_analysis)
        results = pool.map(statistical_analysis_thread, range(iterations))

    return results


def _statistical_analysis(base_result,
                          cluster_interactions,
                          counts,
                          interactions,
                          meta,
                          complex_composition: pd.DataFrame,
                          separator,
                          real_mean_analysis: pd.DataFrame,
                          iteration_number) -> pd.DataFrame:
    """
    Shuffles meta dataset and calculates calculates the means
    """
    shuffled_meta = shuffle_meta(meta)
    shuffled_clusters = build_clusters(shuffled_meta,
                                       counts,
                                       complex_composition,
                                       skip_percent=True)

    shuffled_mean_analysis = mean_analysis(interactions,
                                           shuffled_clusters,
                                           cluster_interactions,
                                           separator)

    result_mean_analysis = shuffled_greater_than_real(shuffled_mean_analysis,
                                                    real_mean_analysis)
    return result_mean_analysis


def shuffled_greater_than_real(shuffled_mean_analysis: pd.DataFrame,
                             real_mean_analysis: pd.DataFrame):
    return np.packbits(shuffled_mean_analysis.values > real_mean_analysis.values, axis=None)


def build_percent_result(real_mean_analysis: pd.DataFrame, real_percents_analysis: pd.DataFrame,
                         statistical_mean_analysis: list, interactions: pd.DataFrame, cluster_interactions: list,
                         base_result: pd.DataFrame, separator: str) -> pd.DataFrame:
    """
    Calculates the pvalues after statistical analysis.

    If real_percent or real_mean are zero, result_percent is 1

    If not:
    Calculates how many shuffled means are bigger than real mean and divides it for the number of
    the total iterations

    EXAMPLE:
        INPUT:

        real_mean_analysis:
                      cluster1_cluster1   cluster1_cluster2 ...
        interaction1  0.5                 0.4
        interaction2  0.0                 0.2


        real_percents_analysis:
                      cluster1_cluster1   cluster1_cluster2 ...
        interaction1  1                   0
        interaction2  0                   1

        statistical means:
        [
                        cluster1_cluster1   cluster1_cluster2 ...
        interaction1    0.6                 0.1
        interaction2    0.0                 0.2

        ,
                      cluster1_cluster1   cluster1_cluster2 ...
        interaction1  0.5                 0.4
        interaction2  0.0                 0.6
        ]

        iterations = 2


        RESULT:

                        cluster1_cluster1   cluster1_cluster2 ...
        interaction1    1                   1
        interaction2    1                   0.5


    """
    core_logger.info('Building Pvalues result')
    percent_result = np.zeros(real_mean_analysis.shape)
    result_size = percent_result.size
    result_shape = percent_result.shape

    for statistical_mean in statistical_mean_analysis:
        percent_result += np.unpackbits(statistical_mean, axis=None)[:result_size].reshape(result_shape)
    percent_result /= len(statistical_mean_analysis)

    mask = (real_mean_analysis.values == 0) | (real_percents_analysis == 0)

    percent_result[mask] = 1

    return pd.DataFrame(percent_result, index=base_result.index, columns=base_result.columns)


def interacting_pair_build(interactions: pd.DataFrame) -> pd.Series:
    """
    Returns the interaction result formated with prefixes
    """

    def get_interactor_name(interaction: pd.Series, suffix: str) -> str:
        if interaction['is_complex{}'.format(suffix)]:
            return interaction['name{}'.format(suffix)]

        return interaction['gene_name{}'.format(suffix)]

    interacting_pair = interactions.apply(
        lambda interaction: '{}_{}'.format(get_interactor_name(interaction, '_1'),
                                           get_interactor_name(interaction, '_2')), axis=1)

    interacting_pair.rename('interacting_pair', inplace=True)

    return interacting_pair


def build_significant_means(real_mean_analysis: pd.DataFrame,
                            result_percent: pd.DataFrame,
                            min_significant_mean: float) -> (pd.Series, pd.DataFrame):
    """
    Calculates the significant means and add rank (number of non empty entries divided by total entries)
    """
    significant_means = get_significant_means(real_mean_analysis, result_percent, min_significant_mean)
    significant_mean_rank = significant_means.count(axis=1)  # type: pd.Series
    number_of_clusters = len(significant_means.columns)
    significant_mean_rank = significant_mean_rank.apply(lambda rank: rank / number_of_clusters)
    significant_mean_rank = significant_mean_rank.round(3)
    significant_mean_rank.name = 'rank'
    return significant_mean_rank, significant_means


def filter_interactions_by_counts(interactions: pd.DataFrame,
                                  counts: pd.DataFrame,
                                  complex_composition: pd.DataFrame) -> pd.DataFrame:
    multidatas = list(counts.index)

    if not complex_composition.empty:
        multidatas += complex_composition['complex_multidata_id'].to_list() + complex_composition[
            'protein_multidata_id'].to_list()

    multidatas = list(set(multidatas))

    def filter_interactions(interaction: pd.Series) -> bool:
        if interaction['multidata_1_id'] in multidatas and interaction['multidata_2_id'] in multidatas:
            return True
        return False

    interactions_filtered = interactions[interactions.apply(filter_interactions, axis=1)]

    return interactions_filtered


def prefilters(interactions: pd.DataFrame,
               counts: pd.DataFrame,
               complexes: pd.DataFrame,
               complex_compositions: pd.DataFrame):
    """
    - Finds the complex defined in counts and calculates their counts values
    - Remove interactions if the simple component ensembl is not in the counts list
    - Remove interactions if the complex component is not in the calculated complex list
    - Remove undefined simple counts
    - Merge simple filtered counts and calculated complex counts
    - Remove duplicated counts
    """
    counts_filtered = filter_empty_cluster_counts(counts)
    complex_composition_filtered, counts_complex = get_involved_complex_from_counts(counts_filtered,
                                                                                    complex_compositions)

    interactions_filtered = filter_interactions_by_counts(interactions,
                                                          counts_filtered,
                                                          complex_composition_filtered)

    counts_simple = filter_counts_by_interactions(counts_filtered, interactions_filtered)

    counts_filtered = counts_simple.append(counts_complex, sort=False)
    counts_filtered = counts_filtered[~counts_filtered.index.duplicated()]

    return interactions_filtered, counts_filtered, complex_composition_filtered

def get_involved_complex_from_counts(counts: pd.DataFrame,
                                     complex_composition: pd.DataFrame) -> (
        pd.DataFrame, pd.DataFrame):
    """
    Finds the complexes defined in counts and calculates the counts values
    """
    proteins_in_complexes = complex_composition['protein_multidata_id'].drop_duplicates().tolist()

    # Remove counts that can't be part of a complex
    counts_filtered = counts[
        counts.apply(lambda count: count.name in proteins_in_complexes, axis=1)]

    # Find complexes with all components defined in counts
    complex_composition_filtered = complex_helper.get_involved_complex_composition_from_protein(counts_filtered,
                                                                                                complex_composition)

    if complex_composition_filtered.empty:
        return complex_composition_filtered, pd.DataFrame(columns=counts.columns)

    available_complex_proteins = complex_composition_filtered['protein_multidata_id'].drop_duplicates().to_list()

    # Remove counts that are not defined in selected complexes
    counts_filtered = counts_filtered[
        counts_filtered.apply(lambda count: count.name in available_complex_proteins, axis=1)]

    return complex_composition_filtered, counts_filtered
