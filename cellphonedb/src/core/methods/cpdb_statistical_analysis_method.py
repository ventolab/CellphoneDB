from cellphonedb.src.core.exceptions.AllCountsFilteredException import AllCountsFilteredException
from cellphonedb.src.core.methods import cpdb_statistical_analysis_complex_method, cpdb_statistical_analysis_helper
from cellphonedb.src.core.exceptions.MissingRequiredArgumentsException import MissingRequiredArgumentsException
from cellphonedb.utils import db_utils, file_utils, scoring_utils
from cellphonedb.src.core.utils import subsampler


def call(cpdb_file_path: str = None,
         meta_file_path: str = None,
         counts_file_path=None,
         counts_data: str = None,
         output_path: str = None,
         microenvs_file_path: str = None,
         active_tfs_file_path: str = None,
         iterations: int = 1000,
         threshold: float = 0.1,
         threads: int = 4,
         debug_seed: int = -1,
         result_precision: int = 3,
         pvalue: float = 0.05,
         subsampling=False,
         subsampling_log=False,
         subsampling_num_pc=100,
         subsampling_num_cells=None,
         separator: str = '|',
         debug: bool = False,
         output_suffix: str = None,
         score_interactions: bool = False
         ) -> dict:
    """Statistical method for analysis

     This method calculates the mean and percent for the cluster interactions
     and for each gene interaction. No shuffling nor DEGs are involved.

     Parameters
     ----------
     cpdb_file_path: str
        CellphoneDB database file path
     meta_file_path: str
         Path to metadata csv file
     counts_file_path:
         Path to counts csv file, or an in-memory AnnData object
     counts_data: str
         Type of gene identifiers in the counts data: "ensembl", "gene_name", "hgnc_symbol"
     output_path: str
        Output path used to store the analysis results (and to store intermediate files when debugging)
     microenvs_file_path: str, optional
         Path to Micro-environment file. Its content is used to limit cluster interactions
     active_tfs_file_path: str, optional
         Path to active TFs. Its content is used to limit cluster interactions.
     iterations: int
        Number of times cell type labels will be shuffled across cells in order to
        determine statistically significant ligand/receptor expression means.
     threshold: float
         Percentage of cells expressing the specific ligand/receptor [0.0 - 1.0]
     threads: int
        Number of threads to be used during the shuffling of clusters/cell types across cells, and in scoring interactions -
        if score_interactions argument was set to True
     debug_seed: int
        This parameter is used for testing only (and only in single-threaded mode
        only - see: https://stackoverflow.com/questions/21494489/what-does-numpy-random-seed0-do).
     result_precision: int
         Number of decimal digits in results.
     pvalue: float
         A p-value below which a ligand/receptor expression mean is considered to be
         statistically significant.
     subsampling: bool
        Enable subsampling
     subsampling_log: bool,
        Enable subsampling log1p for non log-transformed data inputs !!mandatory!!
     subsampling_num_pc: int,
        Subsampling NumPC argument (number of PCs to use) [100]
     subsampling_num_cells: int
        Number of cells to subsample to [1/3 of cells]
     separator: str
         Separator for pairs of genes (gene1|gene2) and clusters (cluster1|cluster2).
     debug: bool
         Storge intermediate data as pickle file (debug_intermediate.pkl).
     output_suffix: str, optional
         Suffix to append to the result file's name (if not provided, timestamp will be used)
    score_interactions: bool
        If True, CellphoneDB interactions will be scored per cell type pair, and returned in interaction_scores_dict
     Returns
     -------
     Dict with the following keys:
         - deconvoluted
         - deconvoluted_percents,
         - means
         - pvalues
         - significant_means
         - interaction_scores
     """
    # Report error unless the required arguments have been provided
    required_arguments = [cpdb_file_path, meta_file_path, counts_data, output_path]
    if None in required_arguments or '' in required_arguments or not counts_file_path:
        raise MissingRequiredArgumentsException(description="All of the following arguments need to be provided: {}"
                                                .format("cpdb_file_path, meta_file_path, counts_file_path, " +
                                                        "counts_data, output_path"))

    # Load into memory CellphoneDB data
    interactions, genes, complex_composition, complex_expanded, gene_synonym2gene_name, receptor2tfs = \
        db_utils.get_interactions_genes_complex(cpdb_file_path)

    # Load user files into memory
    counts, meta, microenvs, degs, active_tf2cell_types = file_utils.get_user_files(
        counts=counts_file_path, meta_fp=meta_file_path, microenvs_fp=microenvs_file_path,
        active_tfs_fp=active_tfs_file_path,
        gene_synonym2gene_name=gene_synonym2gene_name, counts_data=counts_data)

    # add multidata id and means to counts
    counts, counts_relations = cpdb_statistical_analysis_helper.add_multidata_and_means_to_counts(
        counts, genes, counts_data)

    if counts.empty:
        raise AllCountsFilteredException(hint='Are you using human data?')

    # Note that for consistency with the other analysis methods, interaction scoring is done on all cells
    # whether subsampling takes place or not
    counts4scoring = counts.copy()
    # Subsample counts data, if required
    if subsampling:
        ss = subsampler.Subsampler(log=subsampling_log, num_pc=subsampling_num_pc,
                                   num_cells=subsampling_num_cells, verbose=False, debug_seed=None)
        counts = ss.subsample(counts)

    analysis_result = cpdb_statistical_analysis_complex_method.call(meta.copy(),
                                                                    counts,
                                                                    counts_relations,
                                                                    counts_data,
                                                                    active_tf2cell_types,
                                                                    interactions,
                                                                    genes,
                                                                    complex_expanded,
                                                                    complex_composition,
                                                                    microenvs,
                                                                    receptor2tfs,
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

    significant_means = analysis_result['significant_means']
    max_rank = significant_means['rank'].max()
    significant_means['rank'] = significant_means['rank'].apply(lambda rank: rank if rank != 0 else (1 + max_rank))
    significant_means.sort_values('rank', inplace=True)
    means_result = analysis_result['means']

    if score_interactions:
        # Make sure all cell types are strings
        meta['cell_type'] = meta['cell_type'].apply(str)
        interaction_scores = scoring_utils.score_interactions_based_on_participant_expressions_product(
            cpdb_file_path, counts4scoring, means_result.copy(), separator, meta, threshold, "cell_type", threads)
        analysis_result['interaction_scores'] = interaction_scores

    file_utils.save_dfs_as_tsv(output_path, output_suffix, "statistical_analysis", analysis_result)

    return analysis_result
