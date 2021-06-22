from functools import partial
import os
import pandas as pd
import numpy as np
import pickle

from cellphonedb.src.core.core_logger import core_logger
from cellphonedb.src.core.exceptions.AllCountsFilteredException import AllCountsFilteredException
from cellphonedb.src.core.exceptions.NoInteractionsFound import NoInteractionsFound
from cellphonedb.src.core.methods import cpdb_degs_analysis_helper


def call(meta: pd.DataFrame,
         counts: pd.DataFrame,
         degs: pd.DataFrame,
         counts_data: str,
         interactions: pd.DataFrame,
         genes: pd.DataFrame,
         complexes: pd.DataFrame,
         complex_compositions: pd.DataFrame,
         microenvs: pd.DataFrame,
         pvalue: float,
         separator: str,
         iterations: int = 1000,
         threshold: float = 0.1,
         threads: int = 4,
         debug_seed: int = -1,
         result_precision: int = 3,
         ) -> (pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame):
    core_logger.info(
        '[Cluster DEGs Analysis] '
        'Threshold:{} Iterations:{} Debug-seed:{} Threads:{} Precision:{}'.format(threshold,
                                                                                  iterations,
                                                                                  debug_seed,
                                                                                  threads,
                                                                                  result_precision))
    if debug_seed >= 0:
        np.random.seed(debug_seed)
        core_logger.warning('Debug random seed enabled. Setted to {}'.format(debug_seed))
    cells_names = sorted(counts.columns)
    
    interactions.set_index('id_interaction', drop=True, inplace=True)
    
    interactions_reduced = interactions[['multidata_1_id', 'multidata_2_id']].drop_duplicates()

    complex_compositions.set_index('id_complex_composition', inplace=True, drop=True)
    
    # Add id multidata to counts input
    counts: pd.DataFrame = counts.merge(genes[['id_multidata', 'ensembl', 'gene_name', 'hgnc_symbol']],
                                        left_index=True, right_on=counts_data)
    counts_relations = counts[['id_multidata', 'ensembl', 'gene_name', 'hgnc_symbol']].copy()
    
    counts.set_index('id_multidata', inplace=True, drop=True)
    counts = counts[cells_names]
    if np.any(counts.dtypes.values != np.dtype('float32')):
        counts = counts.astype(np.float32)
    counts = counts.groupby(counts.index).mean()

    if counts.empty:
        raise AllCountsFilteredException(hint='Are you using human data?')
    # End add id multidata


    interactions_filtered, counts_filtered, complex_composition_filtered = \
        cpdb_degs_analysis_helper.prefilters(interactions_reduced,
                                                    counts,
                                                    complexes,
                                                    complex_compositions)

    if interactions_filtered.empty:
        raise NoInteractionsFound()

    meta = meta.loc[counts.columns]

    clusters = cpdb_degs_analysis_helper.build_clusters(meta, counts_filtered, complex_composition_filtered, skip_percent=False)
    
    core_logger.info('Running Real Analysis')
    core_logger.debug('Generating cluster combinations')
    cluster_interactions = cpdb_degs_analysis_helper.get_cluster_combinations(clusters['names'], microenvs)
    core_logger.debug(f'Got {len(cluster_interactions)} cluster interactions combinations')
    
    
    core_logger.debug('Build empty result matrix')
    base_result = cpdb_degs_analysis_helper.build_result_matrix(interactions_filtered,
                                                                cluster_interactions,
                                                                separator)

    core_logger.debug('Run Mean Analyisis (real mean)')
    real_mean_analysis = cpdb_degs_analysis_helper.mean_analysis(interactions_filtered,
                                                                clusters,
                                                                cluster_interactions,
                                                                separator)

    # filter real_mean_analysis using the threshold value:
    #   if mean > threashold then value of [gene1|gene2][cluster1|cluster2] = 0; otherwise = 1
    # in this case 0 is 'True' and 1 is 'False'
    core_logger.debug('Run Percent Analyisis (real percent)')
    real_percents_analysis = cpdb_degs_analysis_helper.percent_analysis(clusters,
                                                                        threshold,
                                                                        interactions_filtered,
                                                                        cluster_interactions,
                                                                        separator)
    # interexchange 0s and 1s
    real_percents_analysis = real_percents_analysis.replace(0,2).replace(1,0).replace(2,1)

    save_intermediate_results = True
    if save_intermediate_results:
        with open("debug_intermediate.pkl","wb") as fh:
            pickle.dump( {
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
                "real_mean_analysis": real_mean_analysis,
                "real_percents_analysis": real_percents_analysis}, fh )

        # os.makedirs("_debug", exist_ok=True)
        # genes.to_csv("_debug/genes.csv")
        # interactions.to_csv("_debug/interactions.csv")
        # interactions_reduced.to_csv("_debug/interactions_reduced.csv")
        # complex_compositions.to_csv("_debug/complex_compositions.csv")    
        # counts.to_csv("_debug/counts.csv")
        # interactions_filtered.to_csv("_debug/interactions_filtered.csv")
        # clusters["means"].to_csv("_debug/clusters__means.csv")
        # clusters["percents"].to_csv("_debug/clusters__percents.csv")
        # np.savetxt("_debug/cluster_interactions.csv",cluster_interactions, fmt="%s",delimiter=",")
        # base_result.to_csv("_debug/base_result.csv")
        # real_mean_analysis.to_csv("_debug/real_mean_analysis.csv")
        # real_percents_analysis.to_csv("_debug/real_percents_analysis.csv")


    core_logger.info('Running DEGs-based Analysis')

    core_logger.debug('Run DEGs Analyisis')
    degs_analysis = cpdb_degs_analysis_helper.degs_analysis(degs,genes,
                                                            interactions_filtered,
                                                            cluster_interactions,
                                                            counts_data,
                                                            separator)

    core_logger.debug('Building relevant interactions (merge percent & DEGs analyisis)')
    relevant_interactions =  pd.DataFrame(real_percents_analysis.values&degs_analysis.values,
        columns=real_percents_analysis.columns, index=real_percents_analysis.index)

    core_logger.debug('Building results')
    relevant_interactions_result, means_result, significant_means, deconvoluted_result = build_results(
        interactions_filtered,
        interactions,
        counts_relations,
        real_mean_analysis,
        relevant_interactions,
        clusters['means'],
        complex_composition_filtered,
        counts,
        genes,
        result_precision,
        pvalue,
        counts_data
    )
    
    if save_intermediate_results:
        with open("debug_results.pkl","wb") as fh:
                    pickle.dump( {
                        "relevant_interactions_result": relevant_interactions_result,
                        "means_result": means_result,
                        "significant_means": significant_means,
                        "deconvoluted_result": deconvoluted_result,
                        "complex_compositions": complex_compositions,
                        }, fh )

    return relevant_interactions_result, means_result, significant_means, deconvoluted_result


def build_results(interactions: pd.DataFrame,
                  interactions_original: pd.DataFrame,
                  counts_relations: pd.DataFrame,
                  real_mean_analysis: pd.DataFrame,
                  relevant_interactions: pd.DataFrame,
                  clusters_means: pd.DataFrame,
                  complex_compositions: pd.DataFrame,
                  counts: pd.DataFrame,
                  genes: pd.DataFrame,
                  result_precision: int,
                  pvalue: float,
                  counts_data: str
                  ) -> (pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame):
    """
    Sets the results data structure from method generated data. Results documents are defined by specs.
    """
    core_logger.info('Building results')
    interactions: pd.DataFrame = interactions_original.loc[interactions.index]
    interactions['interaction_index'] = interactions.index
    interactions = interactions.merge(counts_relations, how='left', left_on='multidata_1_id', right_on='id_multidata', )
    interactions = interactions.merge(counts_relations, how='left', left_on='multidata_2_id', right_on='id_multidata',
                                      suffixes=('_1', '_2'))
    interactions.set_index('interaction_index', inplace=True, drop=True)

    interacting_pair = cpdb_degs_analysis_helper.interacting_pair_build(interactions)

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

    significant_mean_rank, significant_means = cpdb_degs_analysis_helper.build_significant_means(
        real_mean_analysis, relevant_interactions)
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

    real_mean_analysis = real_mean_analysis.round(result_precision)
    significant_means = significant_means.round(result_precision)

    # Round result decimals
    for key, cluster_means in clusters_means.items():
        clusters_means[key] = cluster_means.round(result_precision)

    # Document 1
    relevant_interactions_result = pd.concat([interactions_data_result, relevant_interactions], axis=1, join='inner', sort=False)

    # Document 2
    means_result = pd.concat([interactions_data_result, real_mean_analysis], axis=1, join='inner', sort=False)

    # Document 3
    significant_means_result = pd.concat([interactions_data_result, significant_mean_rank, significant_means], axis=1,
                                         join='inner', sort=False)

    # Document 5
    deconvoluted_result = deconvoluted_complex_result_build(clusters_means,
                                                            interactions,
                                                            complex_compositions,
                                                            counts,
                                                            genes,
                                                            counts_data)

    return relevant_interactions_result, means_result, significant_means_result, deconvoluted_result


def deconvoluted_complex_result_build(clusters_means: pd.DataFrame,
                                      interactions: pd.DataFrame,
                                      complex_compositions: pd.DataFrame,
                                      counts: pd.DataFrame,
                                      genes: pd.DataFrame,
                                      counts_data: str) -> pd.DataFrame:
    genes_counts = list(counts.index)
    genes_filtered = genes[genes['id_multidata'].apply(lambda gene: gene in genes_counts)]

    deconvoluted_complex_result_1 = deconvolute_complex_interaction_component(complex_compositions,
                                                                              genes_filtered,
                                                                              interactions,
                                                                              '_1',
                                                                              counts_data)
    deconvoluted_simple_result_1 = deconvolute_interaction_component(interactions,
                                                                     '_1',
                                                                     counts_data)

    deconvoluted_complex_result_2 = deconvolute_complex_interaction_component(complex_compositions,
                                                                              genes_filtered,
                                                                              interactions,
                                                                              '_2',
                                                                              counts_data)
    deconvoluted_simple_result_2 = deconvolute_interaction_component(interactions,
                                                                     '_2',
                                                                     counts_data)

    deconvoluted_result = deconvoluted_complex_result_1.append(
        [deconvoluted_simple_result_1, deconvoluted_complex_result_2, deconvoluted_simple_result_2], sort=False)

    deconvoluted_result.set_index('multidata_id', inplace=True, drop=True)

    deconvoluted_columns = ['gene_name', 'name', 'is_complex', 'protein_name', 'complex_name', 'id_cp_interaction',
                            'gene']

    deconvoluted_result = deconvoluted_result[deconvoluted_columns]
    deconvoluted_result.rename({'name': 'uniprot'}, axis=1, inplace=True)
    deconvoluted_result = pd.concat([deconvoluted_result, clusters_means.reindex(deconvoluted_result.index)], axis=1, join='inner', sort=False)
    deconvoluted_result.set_index('gene', inplace=True, drop=True)
    deconvoluted_result.drop_duplicates(inplace=True)

    return deconvoluted_result


def deconvolute_interaction_component(interactions, suffix, counts_data):
    interactions = interactions[~interactions['is_complex{}'.format(suffix)]]
    deconvoluted_result = pd.DataFrame()
    deconvoluted_result['gene'] = interactions['{}{}'.format(counts_data, suffix)]

    deconvoluted_result[
        ['multidata_id', 'protein_name', 'gene_name', 'name', 'is_complex', 'id_cp_interaction', 'receptor']] = \
        interactions[
            ['multidata{}_id'.format(suffix), 'protein_name{}'.format(suffix), 'gene_name{}'.format(suffix),
             'name{}'.format(suffix),
             'is_complex{}'.format(suffix), 'id_cp_interaction', 'receptor{}'.format(suffix)]]
    deconvoluted_result['complex_name'] = np.nan

    return deconvoluted_result


def deconvolute_complex_interaction_component(complex_compositions,
                                              genes_filtered,
                                              interactions,
                                              suffix,
                                              counts_data):
    return_properties = [counts_data, 'protein_name', 'gene_name', 'name', 'is_complex', 'id_cp_interaction',
                         'receptor', 'complex_name']
    if complex_compositions.empty:
        return pd.DataFrame(
            columns=return_properties)

    deconvoluted_result = pd.DataFrame()
    component = pd.DataFrame()
    component[counts_data] = interactions['{}{}'.format(counts_data, suffix)]
    component[[counts_data, 'protein_name', 'gene_name', 'name', 'is_complex', 'id_cp_interaction', 'id_multidata',
               'receptor']] = \
        interactions[['{}{}'.format(counts_data, suffix), 'protein_name{}'.format(suffix), 'gene_name{}'.format(suffix),
                      'name{}'.format(suffix), 'is_complex{}'.format(suffix), 'id_cp_interaction',
                      'multidata{}_id'.format(suffix), 'receptor{}'.format(suffix)]]

    deconvolution_complex = pd.merge(complex_compositions,
                                     component,
                                     left_on='complex_multidata_id',
                                     right_on='id_multidata')
    deconvolution_complex = pd.merge(deconvolution_complex,
                                     genes_filtered,
                                     left_on='protein_multidata_id',
                                     right_on='protein_multidata_id',
                                     suffixes=['_complex', '_simple'])

    deconvoluted_result['gene'] = deconvolution_complex['{}_simple'.format(counts_data)]

    deconvoluted_result[
        ['multidata_id', 'protein_name', 'gene_name', 'name', 'is_complex', 'id_cp_interaction', 'receptor',
         'complex_name']] = \
        deconvolution_complex[
            ['complex_multidata_id', 'protein_name_simple', 'gene_name_simple', 'name_simple',
             'is_complex_complex', 'id_cp_interaction', 'receptor_simple', 'name_complex']]

    return deconvoluted_result

