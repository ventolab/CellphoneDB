import pickle
from typing import Tuple
import pandas as pd
import numpy as np

from cellphonedb.src.core.core_logger import core_logger
from cellphonedb.src.core.exceptions.AllCountsFilteredException import AllCountsFilteredException
from cellphonedb.src.core.exceptions.NoInteractionsFound import NoInteractionsFound
from cellphonedb.src.core.methods import cpdb_statistical_analysis_helper, cpdb_statistical_analysis_complex_method

def call(meta: pd.DataFrame,
         counts: pd.DataFrame,
         degs: pd.DataFrame,
         counts_data: str,
         interactions: pd.DataFrame,
         genes: pd.DataFrame,
         complexes: pd.DataFrame,
         complex_compositions: pd.DataFrame,
         microenvs: pd.DataFrame = None,
         separator: str = "|",
         iterations: int = 1000,
         threshold: float = 0.1,
         threads: int = 4,
         debug_seed: int= -1,
         result_precision: int = 3,
         debug: bool = False,
         output_path: str = ''
         ) -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Differentially Expressed Genes (DEGs) analysis

    This analysis bypass previous statistical analysis where mean's pvalues are
    computed using a permutation approach. Instead of deriving pvalues from a
    re-shufling strategy to identify relevant means (aka, mean expression of
    ligand/receptor in a cell-cell pair), relevant interactions are identified
    from a list of differentially expressed genes (DEGs) provided by the user
    and computed from their counts matrix.

    Parameters
    ----------
    meta: pd.DataFrame
        Meta data.
    counts: pd.DataFrame
        Counts data.
    degs: pd.DataFrame
        DEGs data.
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
    microenvs: pd.DataFrame, optional
        Microenvironment data to limit cluster interactions
    separator: str, optional
        Separator for pairs of genes (gene1|gene2) and clustes (cluster1|cluster2).
    threshold: float, optional
        Percentage of cells expressing the specific ligand/receptor [0.0 - 1.0]
    result_precision: int, optional
        Number of decimal digits in results.
    debug: bool, optional
        Storge intermediate data as pickle file (debug_intermediate.pkl).
    output_path: str, optional
        Output path used to store intermediate files when debuging.

    Returns
    -------
    means_result DataFrame 
        Result of the means
    significant_means DataFrame
        Result of the significant means
    deconvoluted_result DataFrame
        Result of the deconvoluted complex
    """
    core_logger.info(
        '[Cluster DEGs Analysis] '
        'Threshold:{} Iterations:{} Debug-seed:{} Threads:{} Precision:{}'.format(threshold,
                                                                                  iterations,
                                                                                  debug_seed,
                                                                                  threads,
                                                                                  result_precision))
    core_logger.warning("""
***********************************
DEGs ANALYSIS IS AN EXPERIMENTAL METHOD STILL UNDER DEVELOPMENT!
***********************************""")

    if debug_seed >= 0:
        np.random.seed(debug_seed)
        core_logger.warning('Debug random seed enabled. Set to {}'.format(debug_seed))

    # get reduced interactions (drop duplicateds)
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
        raise NoInteractionsFound()

    meta = meta.loc[counts.columns]

    clusters = cpdb_statistical_analysis_helper.build_clusters(
                                                            meta,
                                                            counts_filtered,
                                                            complex_composition_filtered,
                                                            skip_percent=False)

    core_logger.info('Running Real Analysis')
    core_logger.debug('Generating cluster combinations')
    cluster_interactions = cpdb_statistical_analysis_helper.get_cluster_combinations(clusters['names'], microenvs)
    
    
    core_logger.debug('Build empty result matrix')
    base_result = cpdb_statistical_analysis_helper.build_result_matrix(interactions_filtered,
                                                                cluster_interactions,
                                                                separator)

    core_logger.debug('Run Mean Analyisis (real mean)')
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
    raw_degs = degs.copy()
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
        counts_data
    )

    max_rank = significant_means['rank'].max()
    significant_means['rank'] = significant_means['rank'].apply(lambda rank: rank if rank != 0 else (1 + max_rank))
    significant_means.sort_values('rank', inplace=True)

    return deconvoluted_result, means_result, relevant_interactions_result, significant_means


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
                  counts_data: str) -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """
    Sets the results data structure from method generated data. Results documents are defined by specs.
    """
    core_logger.info('Building results')
    interactions = interactions_original.loc[interactions.index]
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

    # Document 1: relevant_intearcitons.txt
    # drop irrelevant interactions (all zeros)
    relevant_interactions = relevant_interactions.loc[
        (relevant_interactions!=0).any(axis=1)]
    # drop irrelvant clusters (columns with all zeros)
    relevant_interactions = relevant_interactions.loc[
        :,(relevant_interactions!=0).any(axis=0)]
    # concat interactions data and relevant interactions data
    relevant_interactions_result = pd.concat(
        [interactions_data_result, relevant_interactions], axis=1, join='inner', sort=False)

    # Document 2: means.txt
    means_result = pd.concat(
        [interactions_data_result, real_mean_analysis], axis=1, join='inner', sort=False)

    # Document 3: significant_means.txt
    significant_means_result = pd.concat(
        [interactions_data_result, significant_mean_rank, significant_means], axis=1,join='inner', sort=False)

    # Document 4: deconvoluted.txt
    deconvoluted_result = cpdb_statistical_analysis_complex_method.deconvoluted_complex_result_build(clusters_means,
                                                            interactions,
                                                            complex_compositions,
                                                            counts,
                                                            genes,
                                                            counts_data)
    

    return relevant_interactions_result, means_result, significant_means_result, deconvoluted_result

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
    d["deg"]=1

    # create matrix from input and add id_multidata
    d = pd.pivot_table(d, values="deg", index="gene", columns="cluster", fill_value=0)
    d = d.merge(genes[['id_multidata', 'ensembl', 'gene_name', 'hgnc_symbol']],left_index=True,right_on=counts_data)
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
    cluster_degs = cluster_degs.append(complex_degs)
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
