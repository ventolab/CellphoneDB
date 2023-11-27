import pickle
import pandas as pd

from cellphonedb.src.core.core_logger import core_logger
from cellphonedb.src.core.exceptions.AllCountsFilteredException import AllCountsFilteredException
from cellphonedb.src.core.exceptions.MissingRequiredArgumentsException import MissingRequiredArgumentsException
from cellphonedb.src.core.methods import cpdb_statistical_analysis_helper, cpdb_statistical_analysis_complex_method
from cellphonedb.src.core.models.complex import complex_helper
from cellphonedb.utils import db_utils, file_utils, scoring_utils
from cellphonedb.src.core.utils import cellsign


def call(cpdb_file_path: str = None,
         meta_file_path: str = None,
         counts_file_path=None,
         degs_file_path: str = None,
         counts_data: str = None,
         output_path: str = None,
         microenvs_file_path: str = None,
         active_tfs_file_path: str = None,
         separator: str = "|",
         threshold: float = 0.1,
         result_precision: int = 3,
         debug: bool = False,
         output_suffix: str = None,
         score_interactions: bool = False,
         threads: int = 4
         ) -> dict:
    """Differentially Expressed Genes (DEGs) analysis

    This analysis bypass previous statistical analysis where mean's pvalues are
    computed using a permutation approach. Instead of deriving pvalues from a
    re-shufling strategy to identify relevant means (aka, mean expression of
    ligand/receptor in a cell-cell pair), relevant interactions are identified
    from a list of differentially expressed genes (DEGs) provided by the user
    and computed from their counts matrix.

    Parameters
    ----------
    cpdb_file_path: str
        CellphoneDB database file path
    meta_file_path: str
        Path to metadata csv file
    counts_file_path:
        Path to counts csv file, or an in-memory AnnData object
    degs_file_path: str
        Path to differential expression csv file
    counts_data: str
        Type of gene identifiers in the counts data: "ensembl", "gene_name", "hgnc_symbol"
    output_path: str
        Output path used to store the analysis results (and to store intermediate files when debugging)
    microenvs_file_path: str, optional
        Path to Micro-environment file. Its content is used to limit cluster interactions
    active_tfs_file_path: str, optional
         Path to active TFs. Its content is used to limit cluster interactions.
    separator: str, optional
        Separator for pairs of genes (gene1|gene2) and clusters (cluster1|cluster2).
    threshold: float, optional
        Percentage of cells expressing the specific ligand/receptor [0.0 - 1.0]
    result_precision: int, optional
        Number of decimal digits in results.
    debug: bool, optional
        Storge intermediate data as pickle file (debug_intermediate.pkl).
    output_suffix: str, optional
        Suffix to append to the result file's name (if not provided, timestamp will be used)
    score_interactions: bool
        If True, CellphoneDB interactions will be scored per cell type pair, and returned in interaction_scores_dict
    threads: int
        Number of threads to be used when scoring interactions

    Return
    -------
    Dict
         - deconvoluted_result
         - deconvoluted_percents
         - means_result
         - relevant_interactions_result
         - significant_means
         - interaction_scores_dict
    """
    analysis_result = {}
    core_logger.info(
        '[Cluster DEGs Analysis] '
        'Threshold:{} Precision:{}'.format(threshold, result_precision))

    # Report error unless the required arguments have been provided
    required_arguments = [cpdb_file_path, meta_file_path, degs_file_path, counts_data, output_path]
    if None in required_arguments or '' in required_arguments or not counts_file_path:
        raise MissingRequiredArgumentsException(description="All of the following arguments need to be provided: {}".format(
            "cpdb_file_path, meta_file_path, counts_file_path, degs_file_path, counts_data, output_path"))

    # Load into memory CellphoneDB data
    interactions, genes, complex_compositions, complexes, gene_synonym2gene_name, receptor2tfs = \
        db_utils.get_interactions_genes_complex(cpdb_file_path)

    # Load user files into memory
    counts, meta, microenvs, degs, active_tf2cell_types = file_utils.get_user_files(
        counts=counts_file_path, meta_fp=meta_file_path, microenvs_fp=microenvs_file_path, degs_fp=degs_file_path,
        active_tfs_fp=active_tfs_file_path,
        gene_synonym2gene_name=gene_synonym2gene_name, counts_data=counts_data)

    # get reduced interactions (drop duplicates)
    interactions_reduced = interactions[['multidata_1_id', 'multidata_2_id']].drop_duplicates()

    # add multidata id and means to counts
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
        core_logger.info('No CellphoneDB interactions found in this input.')
        return analysis_result

    meta = meta.loc[counts.columns]
    # Make sure all cell types are strings
    meta['cell_type'] = meta['cell_type'].apply(str)
    if not microenvs.empty:
        microenvs['cell_type'] = microenvs['cell_type'].apply(str)
    degs['cluster'] = degs['cluster'].apply(str)

    complex_to_protein_row_ids = complex_helper.map_complex_to_protein_row_ids(complex_composition_filtered, counts_filtered)
    clusters = cpdb_statistical_analysis_helper.build_clusters(
                                                            meta,
                                                            counts_filtered,
                                                            complex_to_protein_row_ids,
                                                            skip_percent=False)

    core_logger.info('Running Real Analysis')
    core_logger.debug('Generating cluster combinations')
    cluster_interactions = cpdb_statistical_analysis_helper.get_cluster_combinations(clusters['names'], microenvs)

    core_logger.debug('Build empty result matrix')
    base_result = cpdb_statistical_analysis_helper.build_result_matrix(interactions_filtered,
                                                                       cluster_interactions,
                                                                       separator)

    core_logger.debug('Run Mean Analysis (real mean)')
    real_mean_analysis = cpdb_statistical_analysis_helper.mean_analysis(interactions_filtered,
                                                                        clusters,
                                                                        cluster_interactions,
                                                                        separator)
    # do percent_analysis using the threshold value:
    core_logger.debug('Run Percent Analysis (real percent)')
    real_percents_analysis = cpdb_statistical_analysis_helper.percent_analysis(clusters,
                                                                               threshold,
                                                                               interactions_filtered,
                                                                               cluster_interactions,
                                                                               separator)

    core_logger.info('Running DEGs-based Analysis')

    # Prepare DEGs matrix from input file
    degs = build_degs_matrix(degs,
                             genes,
                             meta,
                             counts_filtered,
                             complex_composition_filtered,
                             counts_data)

    real_degs_analysis = degs_analysis(degs,
                                       interactions_filtered,
                                       cluster_interactions,
                                       separator)

    core_logger.debug('Building relevant interactions (intersect percent & DEGs analysis)')
    # get relevant interactions by intersecting percent_analysis and DEGs result
    relevant_interactions = pd.DataFrame(real_percents_analysis.values & real_degs_analysis.values,
                                         columns=real_degs_analysis.columns,
                                         index=real_degs_analysis.index)

    if debug:
        core_logger.info('Saving intermediate data to file debug_intermediate.pkl')
        with open(f"{output_path}/debug_intermediate.pkl", "wb") as fh:
            pickle.dump({
                "meta": meta,
                "genes": genes,
                "interactions": interactions,
                "interactions_filtered": interactions_filtered,
                "interactions_reduced": interactions_reduced,
                "complex_compositions": complex_compositions,
                "complex_composition_filtered": complex_composition_filtered,
                "counts": counts,
                "counts_filtered": counts_filtered,
                "counts_relations": counts_relations,
                "clusters": clusters,
                "cluster_interactions": cluster_interactions,
                "base_result": base_result,
                "real_mean_analysis": real_mean_analysis,
                "real_percent_analysis": real_percents_analysis,
                "degs": degs,
                "real_degs_analysis": real_degs_analysis,
                "relevant_interactions": relevant_interactions}, fh)

    core_logger.debug('Building results')

    analysis_result = build_results(
        interactions_filtered,
        interactions,
        counts_relations,
        real_mean_analysis,
        relevant_interactions,
        clusters['means'],
        clusters['percents'],
        complex_composition_filtered,
        counts,
        genes,
        result_precision,
        counts_data,
        separator,
        active_tf2cell_types,
        receptor2tfs
    )

    significant_means = analysis_result['significant_means']
    max_rank = significant_means['rank'].max()
    significant_means['rank'] = significant_means['rank'].apply(lambda rank: rank if rank != 0 else (1 + max_rank))
    significant_means.sort_values('rank', inplace=True)
    means_result = analysis_result['means']

    if score_interactions:
        interaction_scores = scoring_utils.score_interactions_based_on_participant_expressions_product(
            cpdb_file_path, counts, means_result.copy(), separator, meta, threshold, "cell_type", threads)
        analysis_result['interaction_scores'] = interaction_scores

    file_utils.save_dfs_as_tsv(output_path, output_suffix, "degs_analysis", analysis_result)
    return analysis_result


def build_results(interactions: pd.DataFrame,
                  interactions_original: pd.DataFrame,
                  counts_relations: pd.DataFrame,
                  real_mean_analysis: pd.DataFrame,
                  relevant_interactions: pd.DataFrame,
                  clusters_means: pd.DataFrame,
                  clusters_percents: pd.DataFrame,
                  complex_compositions: pd.DataFrame,
                  counts: pd.DataFrame,
                  genes: pd.DataFrame,
                  result_precision: int,
                  counts_data: str,
                  separator: str,
                  active_tf2cell_types: dict,
                  receptor2tfs: dict
                  ) -> dict:
    """
    Sets the results data structure from method generated data. Results documents are defined by specs.
    """
    core_logger.info('Building results')
    interactions = interactions_original.loc[interactions.index]
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
        real_mean_analysis, relevant_interactions)
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

    # Document 1: relevant_intearcitons.txt
    # drop irrelevant interactions (all zeros)
    relevant_interactions = relevant_interactions.loc[
        (relevant_interactions != 0).any(axis=1)]
    # drop irrelevant clusters (columns with all zeros)
    relevant_interactions = relevant_interactions.loc[
        :, (relevant_interactions != 0).any(axis=0)]
    # concat interactions data and relevant interactions data
    relevant_interactions_result = pd.concat(
        [interactions_data_result, relevant_interactions], axis=1, join='inner', sort=False)

    # Perform CellSign analysis to identify active interactions
    active_interactions, active_interactions_deconvoluted = \
        cellsign.find_active_interactions(relevant_interactions_result, receptor2tfs, active_tf2cell_types, separator)

    # Document 2: means.txt
    means_result = pd.concat(
        [interactions_data_result, real_mean_analysis], axis=1, join='inner', sort=False)

    # Document 3: significant_means.txt
    significant_means_result = pd.concat(
        [interactions_data_result, significant_mean_rank, significant_means], axis=1, join='inner', sort=False)

    # Document 4: deconvoluted.txt
    deconvoluted_result, deconvoluted_percents = \
        cpdb_statistical_analysis_complex_method.deconvoluted_complex_result_build(clusters_means,
                                                                                   clusters_percents,
                                                                                   interactions,
                                                                                   complex_compositions,
                                                                                   counts,
                                                                                   genes,
                                                                                   counts_data)

    analysis_result = {'deconvoluted': deconvoluted_result,
                       'deconvoluted_percents': deconvoluted_percents,
                       'means': means_result,
                       'relevant_interactions': relevant_interactions_result,
                       'significant_means': significant_means_result,
                       'CellSign_active_interactions': active_interactions,
                       'CellSign_active_interactions_deconvoluted': active_interactions_deconvoluted}
    return analysis_result


def build_degs_matrix(degs: pd.DataFrame,
                      genes: pd.DataFrame,
                      meta: pd.DataFrame,
                      counts: pd.DataFrame,
                      complex_composition: pd.DataFrame,
                      counts_data: str) -> pd.DataFrame:
    """
    Reshape the input DEGs into a matrix with genes (id_metadata) as index
    and clusters (cell_types) as columns.

    The prepare process integrates id_multidata from genes into the DEGs
    using the provided gene (2nd column from the input file) and the count_data.
    Then drops possible duplicated genes and adds missing genes to match
    the index. Also adds the complexes calculations.

    Parameters
    ----------
    degs: pd.DataFrame
        DEGs data.
    genes: pd.DataFrame
        Genes from CellPhoneDB database required for DEGs-CellPhoneDB gene mapping
    meta: pd.DataFrame
        Meta data.
    counts: pd.DataFrame
        Multidata required for indexing
    complex_composition: list
        List of cluster names to be used as column names
    counts_data: str
        Type of gene identifiers in the counts data: "ensembl", "gene_name", "hgnc_symbol"
    """

    CELL_TYPE = 'cell_type'
    COMPLEX_ID = 'complex_multidata_id'
    PROTEIN_ID = 'protein_multidata_id'

    meta[CELL_TYPE] = meta[CELL_TYPE].astype('category')
    cluster_names = meta[CELL_TYPE].cat.categories

    # mark provided DEGs as active
    d = degs.copy()
    d["deg"] = 1

    # create matrix from input and add id_multidata
    d = pd.pivot_table(d, values="deg", index="gene", columns="cluster", fill_value=0)
    d = d.merge(genes[['id_multidata', 'ensembl', 'gene_name', 'hgnc_symbol']], left_index=True, right_on=counts_data)
    d.set_index('id_multidata', inplace=True, drop=True)
    d = d.groupby(d.index).max()

    # sort cell names
    cells_names = sorted(d.columns)
    d = d[cells_names]

    # Simple DEGs cluster count
    cluster_degs = pd.DataFrame(0, index=counts.index, columns=cluster_names.to_list())
    cluster_degs.update(d)

    # add complexes and check if they are DEGs (max)
    complexes = complex_composition.groupby(COMPLEX_ID).apply(lambda x: x[PROTEIN_ID].values).to_dict()
    complex_degs = pd.DataFrame(
        {complex_id: cluster_degs.loc[protein_ids].max(axis=0).values
         for complex_id, protein_ids in complexes.items()},
        index=cluster_degs.columns
    ).T
    cluster_degs = pd.concat([cluster_degs, complex_degs])
    return cluster_degs


def degs_analysis(degs: pd.DataFrame,
                  interactions: pd.DataFrame,
                  cluster_interactions: list,
                  separator: str) -> pd.DataFrame:
    """
    Calculates the filter matrix based on the DEGs for the list of interactions and for each cluster
    Sets 1 if one of both is 1.

    Example
    -------
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

    Returns
    ----------
    pd.DataFrame:
        A DataFrame where each column is a cluster combination (cluster1|cluster2)
        and each row represents an interaction (gene1|gene2). Values are the binary:
            - 1 means the interaction is relevant for the cluster because either
              gene1-cluster1 or gene2-cluster2 are in the DEGs input.
            - 0 means there is no relevant interaction pair (missing from DEGs input).
    """

    GENE_ID1 = 'multidata_1_id'
    GENE_ID2 = 'multidata_2_id'

    cluster1_names = cluster_interactions[:, 0]
    cluster2_names = cluster_interactions[:, 1]
    gene1_ids = interactions[GENE_ID1].values
    gene2_ids = interactions[GENE_ID2].values

    x = degs.loc[gene1_ids, cluster1_names].values
    y = degs.loc[gene2_ids, cluster2_names].values

    result = pd.DataFrame(
        ((x + y) > 0).astype(int),
        index=interactions.index,
        columns=(pd.Series(cluster1_names) + separator + pd.Series(cluster2_names)).values)

    return result
