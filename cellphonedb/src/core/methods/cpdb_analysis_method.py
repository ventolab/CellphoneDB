
from typing import Tuple
import pandas as pd
import numpy as np
import pickle

from cellphonedb.src.core.core_logger import core_logger
from cellphonedb.src.core.exceptions.AllCountsFilteredException import AllCountsFilteredException
from cellphonedb.src.core.exceptions.NoInteractionsFound import NoInteractionsFound
from cellphonedb.src.core.methods import cpdb_statistical_analysis_complex_method
from cellphonedb.src.core.methods import cpdb_statistical_analysis_helper


def call(meta: pd.DataFrame,
         counts: pd.DataFrame,
         counts_data: str,
         interactions: pd.DataFrame,
         genes: pd.DataFrame,
         complexes: pd.DataFrame,
         complex_compositions: pd.DataFrame,
         microenvs: pd.DataFrame,
         separator: str = "|",
         threshold: float = 0.1,
         result_precision: int = 3,
         debug: bool = False,
         output_path: str = "",
         ) -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Non-statistical method for analysis

    This methods calculates the mean and percent for the cluster interactions
    and for each gene interaction. No shuffling nor DEGs are involved.

    Parameters
    ----------
    meta: str
        Meta data.
    counts: str
        Counts data.
    counts_data: str
        Type of gene identifiers in the counts data: "ensembl", "gene_name", "hgnc_symbol"
    interactions: pd.DataFrame
        Interactions from CellPhoneDB database
    genes: pd.DataFrame
        Genes from CellPhoneDB database
    complexes: pd.DataFrame
        Complex and Multidata joined from CellPhoneDB database
    complex_compositions: pd.DataFrame
        ComplexComposition from CellPhoneDB database
    microenvs: pd.DataFrame
        Microenvironment data to limit cluster interactions
    separator: str
        Separator for pairs of genes (gene1|gene2) and clustes (cluster1|cluster2).
    threshold: float
        Percentage of cells expressing the specific ligand/receptor [0.0 - 1.0]
    result_precision: int
        Number of decimal digits in results.
    debug: bool
        Storge intermediate data as pickle file (debug_intermediate.pkl).
    output_path: str
        Output path used to store intermediate files when debuging.

    Returns
    -------
    Tuple: A tuple containing:
        - means_result: Result of the means
        - significant_means: Result of the significant means
        - deconvoluted_result: Result of the deconvoluted complex
    """
    core_logger.info(
        '[Non Statistical Method] Threshold:{} Precision:{}'.format(threshold,
                                                                    result_precision))

    # get reduced interactions (drop duplicateds)
    interactions_reduced = interactions[['multidata_1_id', 'multidata_2_id']].drop_duplicates()

    # add id multidata and means to counts input
    counts, counts_relations = cpdb_statistical_analysis_helper.add_multidata_and_means_to_counts(
        counts, genes, counts_data)

    if counts.empty:
        raise AllCountsFilteredException(hint='Are you using human data?')

    interactions_filtered, counts_filtered, complex_composition_filtered = \
        cpdb_statistical_analysis_helper.prefilters(interactions_reduced,
                                                    counts,
                                                    complexes,
                                                    complex_compositions)
    if interactions_filtered.empty:
        raise NoInteractionsFound()

    meta = meta.loc[counts.columns]

    clusters = cpdb_statistical_analysis_helper.build_clusters(meta,
                                                               counts_filtered,
                                                               complex_composition_filtered,
                                                               skip_percent=False)
    core_logger.info('Running Real Analysis')

    cluster_interactions = cpdb_statistical_analysis_helper.get_cluster_combinations(clusters['names'], microenvs)

    base_result = cpdb_statistical_analysis_helper.build_result_matrix(interactions_filtered,
                                                                       cluster_interactions,
                                                                       separator)

    mean_analysis = cpdb_statistical_analysis_helper.mean_analysis(interactions_filtered,
                                                                   clusters,
                                                                   cluster_interactions,
                                                                   separator)

    percent_analysis = cpdb_statistical_analysis_helper.percent_analysis(clusters,
                                                             threshold,
                                                             interactions_filtered,
                                                             cluster_interactions,
                                                             separator)

    if debug:
        with open(f"{output_path}/debug_intermediate.pkl", "wb") as fh:
            pickle.dump({
                "genes": genes,
                "interactions": interactions,
                "interactions_filtered": interactions_filtered,
                "interactions_reduced": interactions_reduced,
                "complex_compositions": complex_compositions,
                "counts": counts,
                "counts_relations": counts_relations,
                "clusters_means_percents": clusters,
                "cluster_interactions": cluster_interactions,
                "base_result": base_result,
                "mean_analysis": mean_analysis,
                "percent_analysis": percent_analysis}, fh)

    means_result, significant_means, deconvoluted_result = build_results(
        interactions_filtered,
        interactions,
        counts_relations,
        mean_analysis,
        percent_analysis,
        clusters['means'],
        complex_composition_filtered,
        counts,
        genes,
        result_precision,
        counts_data
    )

    max_rank = significant_means['rank'].max()
    significant_means['rank'] = significant_means['rank'].apply(lambda rank: rank if rank != 0 else (1 + max_rank))
    significant_means.sort_values('rank', inplace=True)

    return means_result, significant_means, deconvoluted_result


def build_results(interactions: pd.DataFrame,
                  interactions_original: pd.DataFrame,
                  counts_relations: pd.DataFrame,
                  mean_analysis: pd.DataFrame,
                  percent_analysis: pd.DataFrame,
                  clusters_means: pd.DataFrame,
                  complex_compositions: pd.DataFrame,
                  counts: pd.DataFrame,
                  genes: pd.DataFrame,
                  result_precision: int,
                  counts_data: str) -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """
    Sets the results data structure from method generated data. 
    
    Results documents are defined by specs.

    Returns
    -------
    Tuple: A tuple containing the results for:
        - means
        - significant_means
        - deconvoluted
    """
    core_logger.info('Building results')
    interactions: pd.DataFrame = interactions_original.loc[interactions.index]
    interactions['interaction_index'] = interactions.index
    interactions = interactions.merge(counts_relations, how='left', left_on='multidata_1_id', right_on='id_multidata', )
    interactions = interactions.merge(counts_relations, how='left', left_on='multidata_2_id', right_on='id_multidata',
                                      suffixes=('_1', '_2'))
    interactions.set_index('interaction_index', inplace=True, drop=True)

    interacting_pair = cpdb_statistical_analysis_helper.interacting_pair_build(interactions)

    def simple_complex_indicator(interaction: pd.Series, suffix: str) -> str:
        """
        Add simple/complex prefixes to interaction components
        """
        if interaction['is_complex{}'.format(suffix)]:
            return 'complex:{}'.format(interaction['name{}'.format(suffix)])

        return 'simple:{}'.format(interaction['name{}'.format(suffix)])

    interactions['partner_a'] = interactions.apply(lambda interaction: simple_complex_indicator(interaction, '_1'),
                                                   axis=1)
    interactions['partner_b'] = interactions.apply(lambda interaction: simple_complex_indicator(interaction, '_2'),
                                                   axis=1)

    significant_mean_rank, significant_means = cpdb_statistical_analysis_helper.build_significant_means(
        mean_analysis,
        percent_analysis)
    significant_means = significant_means.round(result_precision)

    gene_columns = ['{}_{}'.format(counts_data, suffix) for suffix in ('1', '2')]
    gene_renames = {column: 'gene_{}'.format(suffix) for column, suffix in zip(gene_columns, ['a', 'b'])}

    # Remove useless columns
    interactions_data_result = pd.DataFrame(
        interactions[['id_cp_interaction', 'partner_a', 'partner_b', 'receptor_1', 'receptor_2', *gene_columns,
                      'annotation_strategy']].copy())

    interactions_data_result = pd.concat([interacting_pair, interactions_data_result], axis=1, sort=False)

    interactions_data_result['secreted'] = (interactions['secreted_1'] | interactions['secreted_2'])
    interactions_data_result['is_integrin'] = (interactions['integrin_1'] | interactions['integrin_2'])

    interactions_data_result.rename(
        columns={**gene_renames, 'receptor_1': 'receptor_a', 'receptor_2': 'receptor_b'},
        inplace=True)

    # Dedupe rows and filter only desired columns
    interactions_data_result.drop_duplicates(inplace=True)

    means_columns = ['id_cp_interaction', 'interacting_pair', 'partner_a', 'partner_b', 'gene_a', 'gene_b', 'secreted',
                     'receptor_a', 'receptor_b', 'annotation_strategy', 'is_integrin']

    interactions_data_result = interactions_data_result[means_columns]

    mean_analysis = mean_analysis.round(result_precision)

    # Round result decimals
    for key, cluster_means in clusters_means.items():
        clusters_means[key] = cluster_means.round(result_precision)

    # Document 2
    means_result = pd.concat([interactions_data_result, mean_analysis], axis=1, join='inner', sort=False)

    # Document 3
    significant_means_result = pd.concat([interactions_data_result, significant_mean_rank, significant_means], axis=1,
                                         join='inner', sort=False)

    # Document 5
    deconvoluted_result = cpdb_statistical_analysis_complex_method.deconvoluted_complex_result_build(
        clusters_means,
        interactions,
        complex_compositions,
        counts,
        genes,
        counts_data)

    return means_result, significant_means_result, deconvoluted_result


def deconvoluted_complex_result_build(clusters_means: dict, interactions: pd.DataFrame,
                                      complex_compositions: pd.DataFrame, counts: pd.DataFrame,
                                      genes: pd.DataFrame, counts_data: str) -> pd.DataFrame:
    genes_counts = list(counts.index)
    genes_filtered = genes[genes[counts_data].apply(lambda gene: gene in genes_counts)]

    deconvoluted_complex_result_1 = deconvolute_complex_interaction_component(complex_compositions, genes_filtered,
                                                                              interactions, '_1', counts_data)
    deconvoluted_simple_result_1 = deconvolute_interaction_component(interactions, '_1', counts_data)

    deconvoluted_complex_result_2 = deconvolute_complex_interaction_component(complex_compositions, genes_filtered,
                                                                              interactions, '_2', counts_data)
    deconvoluted_simple_result_2 = deconvolute_interaction_component(interactions, '_2', counts_data)

    deconvoluted_result = deconvoluted_complex_result_1.append(
        [deconvoluted_simple_result_1, deconvoluted_complex_result_2, deconvoluted_simple_result_2], sort=False)

    deconvoluted_result.set_index('gene', inplace=True)

    cluster_counts = pd.DataFrame(index=deconvoluted_result.index)

    for key, cluster_means in clusters_means.items():
        cluster_counts[key] = cluster_means

    cluster_counts = cluster_counts.reindex(sorted(cluster_counts.columns), axis=1)

    # Here we sort and filter unwanted columns
    deconvoluted_columns = ['gene_name', 'name', 'is_complex', 'protein_name', 'complex_name', 'id_cp_interaction']

    deconvoluted_result = deconvoluted_result[deconvoluted_columns]
    deconvoluted_result.rename({'name': 'uniprot'}, axis=1, inplace=True)

    deconvoluted_result = pd.concat([deconvoluted_result, cluster_counts], axis=1, join='inner', sort=False)

    deconvoluted_result.reset_index(inplace=True)
    deconvoluted_result.drop(columns='gene', inplace=True)
    return deconvoluted_result


def deconvolute_interaction_component(interactions, suffix, counts_data):
    interactions = interactions[~interactions['is_complex{}'.format(suffix)]]
    deconvoluted_result = pd.DataFrame()
    deconvoluted_result['gene'] = interactions['{}{}'.format(counts_data, suffix)]

    deconvoluted_result[['protein_name', 'gene_name', 'name', 'is_complex', 'id_cp_interaction', 'receptor']] = \
        interactions[['protein_name{}'.format(suffix), 'gene_name{}'.format(suffix), 'name{}'.format(suffix),
                      'is_complex{}'.format(suffix), 'id_cp_interaction', 'receptor{}'.format(suffix)]]
    deconvoluted_result['complex_name'] = np.nan

    return deconvoluted_result


def deconvolute_complex_interaction_component(complex_compositions, genes_filtered, interactions, suffix, counts_data):
    deconvoluted_result = pd.DataFrame()
    component = pd.DataFrame()
    component[counts_data] = interactions['{}{}'.format(counts_data, suffix)]
    component[['protein_name', 'gene_name', 'name', 'is_complex', 'id_cp_interaction', 'id_multidata', 'receptor']] = \
        interactions[
            ['protein_name{}'.format(suffix), 'gene_name{}'.format(suffix),
             'name{}'.format(suffix), 'is_complex{}'.format(suffix), 'id_cp_interaction',
             'id_multidata{}'.format(suffix), 'receptor{}'.format(suffix)]]

    deconvolution_complex = pd.merge(complex_compositions, component, left_on='complex_multidata_id',
                                     right_on='id_multidata')
    deconvolution_complex = pd.merge(deconvolution_complex, genes_filtered, left_on='protein_multidata_id',
                                     right_on='protein_multidata_id', suffixes=['_complex', '_simple'])

    deconvoluted_result['gene'] = deconvolution_complex['{}_simple'.format(counts_data)]

    deconvoluted_result[
        ['protein_name', 'gene_name', 'name', 'is_complex', 'id_cp_interaction', 'receptor', 'complex_name']] = \
        deconvolution_complex[
            ['protein_name_simple', 'gene_name_simple', 'name_simple',
             'is_complex_complex', 'id_cp_interaction', 'receptor_simple', 'name_complex']]

    return deconvoluted_result
