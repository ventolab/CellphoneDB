import pandas as pd
import numpy as np
import pickle
from cellphonedb.src.core.core_logger import core_logger
from cellphonedb.src.core.methods import cpdb_statistical_analysis_helper
from cellphonedb.src.core.models.complex import complex_helper
from cellphonedb.src.core.utils import cellsign


def call(meta: pd.DataFrame,
         counts: pd.DataFrame,
         counts_relations: pd.DataFrame,
         counts_data: str,
         active_tf2cell_types: dict,
         interactions: pd.DataFrame,
         genes: pd.DataFrame,
         complexes: pd.DataFrame,
         complex_compositions: pd.DataFrame,
         microenvs: pd.DataFrame,
         receptor2tfs: dict,
         pvalue: float,
         separator: str = '|',
         iterations: int = 1000,
         threshold: float = 0.1,
         threads: int = 4,
         debug_seed: int = -1,
         result_precision: int = 3,
         debug: bool = False,
         output_path: str = '',
         ) -> dict:
    analysis_result = {}
    core_logger.info(
        '[Cluster Statistical Analysis] '
        'Threshold:{} Iterations:{} Debug-seed:{} Threads:{} Precision:{}'.format(threshold,
                                                                                  iterations,
                                                                                  debug_seed,
                                                                                  threads,
                                                                                  result_precision))
    if debug_seed >= 0:
        np.random.seed(debug_seed)
        core_logger.warning('Debug random seed enabled. Set to {}'.format(debug_seed))

    # get reduced interactions (drop duplicates)
    interactions_reduced = interactions[['multidata_1_id', 'multidata_2_id']].drop_duplicates()

    interactions_filtered, counts_filtered, complex_composition_filtered = \
        cpdb_statistical_analysis_helper.prefilters(interactions_reduced,
                                                    counts,
                                                    complexes,
                                                    complex_compositions)

    if interactions_filtered.empty:
        core_logger.info('No CellphoneDB interactions found in this input.')
        return analysis_result

    meta = meta.loc[counts.columns]
    # Make sure all cell types are strings
    meta['cell_type'] = meta['cell_type'].apply(str)
    if not microenvs.empty:
        microenvs['cell_type'] = microenvs['cell_type'].apply(str)

    complex_to_protein_row_ids = complex_helper.map_complex_to_protein_row_ids(complex_composition_filtered,
                                                                               counts_filtered)
    clusters = cpdb_statistical_analysis_helper.build_clusters(meta, counts_filtered,
                                                               complex_to_protein_row_ids,
                                                               skip_percent=False)
    core_logger.info('Running Real Analysis')
    cluster_combinations = cpdb_statistical_analysis_helper.get_cluster_combinations(clusters['names'], microenvs)
    base_result = cpdb_statistical_analysis_helper.build_result_matrix(interactions_filtered,
                                                                       cluster_combinations,
                                                                       separator)

    real_mean_analysis = cpdb_statistical_analysis_helper.mean_analysis(interactions_filtered,
                                                                        clusters,
                                                                        cluster_combinations,
                                                                        separator)

    real_percents_analysis = cpdb_statistical_analysis_helper.percent_analysis(clusters,
                                                                               threshold,
                                                                               interactions_filtered,
                                                                               cluster_combinations,
                                                                               separator)

    core_logger.info('Running Statistical Analysis')
    statistical_mean_analysis = cpdb_statistical_analysis_helper.shuffled_analysis(iterations,
                                                                                   meta,
                                                                                   counts_filtered,
                                                                                   interactions_filtered,
                                                                                   cluster_combinations,
                                                                                   complex_to_protein_row_ids,
                                                                                   real_mean_analysis,
                                                                                   threads,
                                                                                   separator)

    result_percent = cpdb_statistical_analysis_helper.build_percent_result(real_mean_analysis,
                                                                           real_percents_analysis,
                                                                           statistical_mean_analysis,
                                                                           interactions_filtered,
                                                                           cluster_combinations,
                                                                           base_result,
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
                "cluster_combinations": cluster_combinations,
                "base_result": base_result,
                "real_mean_analysis": real_mean_analysis,
                "real_percent_analysis": real_percents_analysis,
                "statistical_mean_analysis": statistical_mean_analysis,
                "result_percent": result_percent}, fh)

    analysis_result = build_results(
        interactions_filtered,
        interactions,
        counts_relations,
        real_mean_analysis,
        result_percent,
        clusters['means'],
        clusters['percents'],
        complex_composition_filtered,
        counts,
        genes,
        result_precision,
        pvalue,
        counts_data,
        separator,
        active_tf2cell_types,
        receptor2tfs
    )
    return analysis_result


def build_results(interactions: pd.DataFrame,
                  interactions_original: pd.DataFrame,
                  counts_relations: pd.DataFrame,
                  real_mean_analysis: pd.DataFrame,
                  result_percent: pd.DataFrame,
                  clusters_means: pd.DataFrame,
                  clusters_percents: pd.DataFrame,
                  complex_compositions: pd.DataFrame,
                  counts: pd.DataFrame,
                  genes: pd.DataFrame,
                  result_precision: int,
                  pvalue: float,
                  counts_data: str,
                  separator: str,
                  active_tf2cell_types: dict,
                  receptor2tfs: dict
                  ) -> dict:
    """
    Sets the results data structure from method generated data. Results documents are defined by specs.
    """
    core_logger.info('Building results')
    interactions: pd.DataFrame = interactions_original.loc[interactions.index]
    interactions['interaction_index'] = interactions.index
    interactions = interactions.merge(counts_relations, how='left', left_on='multidata_1_id', right_on='id_multidata', )
    # The column drop below prevents: 'FutureWarning: Passing 'suffixes' which cause duplicate columns {'id_multidata_1'}
    # in the result is deprecated and will raise a MergeError in a future version.'
    interactions = interactions.drop('id_multidata', axis=1)
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
        real_mean_analysis, result_percent, pvalue)
    significant_means = significant_means.round(result_precision)

    gene_columns = ['{}_{}'.format('gene_name', suffix) for suffix in ('1', '2')]
    gene_renames = {column: 'gene_{}'.format(suffix) for column, suffix in zip(gene_columns, ['a', 'b'])}
    # Cater for DB version-dependent column names
    interaction_columns = []
    if 'directionality' in interactions.columns:
        interaction_columns = ['directionality', 'classification']
    # Remove superfluous columns
    interactions_data_result = pd.DataFrame(
        interactions[['id_cp_interaction', 'partner_a', 'partner_b', 'receptor_1', 'receptor_2', *gene_columns,
                      'annotation_strategy'] + interaction_columns].copy())

    interactions_data_result = pd.concat([interacting_pair, interactions_data_result], axis=1, sort=False)

    interactions_data_result['secreted'] = (interactions['secreted_1'] | interactions['secreted_2'])
    interactions_data_result['is_integrin'] = (interactions['integrin_1'] | interactions['integrin_2'])
    interactions_data_result.rename(
        columns={**gene_renames, 'receptor_1': 'receptor_a', 'receptor_2': 'receptor_b'},
        inplace=True)

    # Dedupe rows and filter only desired columns
    interactions_data_result.drop_duplicates(inplace=True)

    means_columns = ['id_cp_interaction', 'interacting_pair', 'partner_a', 'partner_b', 'gene_a', 'gene_b', 'secreted',
                     'receptor_a', 'receptor_b', 'annotation_strategy', 'is_integrin'] + interaction_columns

    interactions_data_result = interactions_data_result[means_columns]
    real_mean_analysis = real_mean_analysis.round(result_precision)
    significant_means = significant_means.round(result_precision)

    # Round result decimals
    for key, cluster_means in clusters_means.items():
        clusters_means[key] = cluster_means.round(result_precision)
    for key, cluster_percents in clusters_percents.items():
        clusters_percents[key] = cluster_percents.round(result_precision)

    # Document 1
    pvalues_result = pd.concat([interactions_data_result, result_percent], axis=1, join='inner', sort=False)

    # Document 2
    means_result = pd.concat([interactions_data_result, real_mean_analysis], axis=1, join='inner', sort=False)

    # Document 3
    significant_means_result = pd.concat([interactions_data_result, significant_mean_rank, significant_means], axis=1,
                                         join='inner', sort=False)

    # Perform CellSign analysis to identify active interactions
    active_interactions, active_interactions_deconvoluted = \
        cellsign.find_active_interactions(significant_means_result, receptor2tfs, active_tf2cell_types, separator)

    # Document 5
    deconvoluted_result, deconvoluted_percents = deconvoluted_complex_result_build(clusters_means,
                                                                                   clusters_percents,
                                                                                   interactions,
                                                                                   complex_compositions,
                                                                                   counts,
                                                                                   genes,
                                                                                   counts_data)
    analysis_result = {'deconvoluted': deconvoluted_result,
                       'deconvoluted_percents': deconvoluted_percents,
                       'means': means_result,
                       'pvalues': pvalues_result,
                       'significant_means': significant_means_result,
                       'CellSign_active_interactions': active_interactions,
                       'CellSign_active_interactions_deconvoluted': active_interactions_deconvoluted}

    return analysis_result


def deconvoluted_complex_result_build(clusters_means: pd.DataFrame,
                                      clusters_percents: pd.DataFrame,
                                      interactions: pd.DataFrame,
                                      complex_compositions: pd.DataFrame,
                                      counts: pd.DataFrame,
                                      genes: pd.DataFrame,
                                      counts_data: str) -> (pd.DataFrame, pd.DataFrame):
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

    deconvoluted_result = pd.concat([deconvoluted_complex_result_1, deconvoluted_simple_result_1,
                                     deconvoluted_complex_result_2, deconvoluted_simple_result_2], sort=False)

    deconvoluted_result.set_index('multidata_id', inplace=True, drop=True)

    deconvoluted_columns = ['gene_name', 'name', 'is_complex', 'protein_name', 'complex_name', 'id_cp_interaction',
                            'gene']

    deconvoluted_result = deconvoluted_result[deconvoluted_columns]
    deconvoluted_result.rename({'name': 'uniprot'}, axis=1, inplace=True)
    deconvoluted_result4pcts = deconvoluted_result.copy(deep=True)
    deconvoluted_result = pd.concat([deconvoluted_result, clusters_means.reindex(deconvoluted_result.index)],
                                    axis=1, join='inner', sort=False)
    deconvoluted_percents = pd.concat([deconvoluted_result4pcts,
                                       clusters_percents.reindex(deconvoluted_result4pcts.index)],
                                      axis=1, join='inner', sort=False)

    deconvoluted_result.drop_duplicates(inplace=True)
    deconvoluted_percents.drop_duplicates(inplace=True)

    return deconvoluted_result, deconvoluted_percents


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

    # The two operations below remove duplicates - this prevents the latest version of Pandas throwing
    # 'ValueError: Columns must be same length as key' when counts_data == 'gene_name'
    component_columns = list(dict.fromkeys([counts_data, 'protein_name', 'gene_name', 'name', 'is_complex',
                                            'id_cp_interaction', 'id_multidata', 'receptor']))
    interactions_columns = list(dict.fromkeys(['{}{}'.format(counts_data, suffix), 'protein_name{}'
                                              .format(suffix), 'gene_name{}'.format(suffix),
                                               'name{}'.format(suffix), 'is_complex{}'.format(suffix),
                                               'id_cp_interaction', 'multidata{}_id'
                                              .format(suffix), 'receptor{}'.format(suffix)]))
    component[component_columns] = interactions[interactions_columns]

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
