import pandas as pd
import numpy as np
from cellphonedb.utils import db_utils
from sklearn.preprocessing import MinMaxScaler
from cellphonedb.src.core.core_logger import core_logger
from collections import ChainMap

from functools import partial
from multiprocessing.pool import Pool
from tqdm.std import tqdm


def filter_genes_per_cell_type(
        matrix: pd.DataFrame,
        metadata: pd.DataFrame,
        min_pct_cell: float,
        cell_column_name: str) -> pd.DataFrame:
    """
    This function takes as input a normalized count matrix and for each gene calculates the
    percentage of cells expressing it (anything but 0). Then sets to 0 the expression of a given
    gene for all cells of a specific cell type if this gene is expressed in less than min_pct_cell.

    Parameters
    ----------
    matrix: Normalized gene expression matrix (genes x barcodes).
    metadata: Index contains the barcode id and a single column named 'cell_type' indicating the group/cell type which
              the barcode belongs to.
    min_pct_cell : Percentage of cells required to express a given gene.
    cell_column_name: Name of the column containing cell types

    Returns
    -------
    pd.DataFrame
        matrix with expression of a gene set to 0 for all cells - if that gene is expressed in less than
        min_pct_cell of cells
    """
    core_logger.info("Scoring interactions: Filtering genes per cell type..")
    matrix = matrix.copy()

    # Cell types present in metadata
    cell_type_data = set(metadata[cell_column_name])

    for cell_type in tqdm(cell_type_data):
        # Obtain the barcode of the cells annotated under cell_type
        idx = metadata[cell_column_name] == cell_type
        cell_barcode = list(metadata.index[idx])

        # Calculate percentage of cells expressing the gene (expression != 0)
        gene_expr_pct = (matrix[cell_barcode] != 0).sum(axis=1) / len(cell_barcode)

        # List of genes lowly expressed
        gene_lowly_expr = list(gene_expr_pct.index[gene_expr_pct < min_pct_cell])

        # Set expression to zero for genes expressed in a given cell type below the
        # user defined min_pct_cell
        matrix.loc[gene_lowly_expr, cell_barcode] = 0

    # Return filtered matrix
    return matrix


def mean_expression_per_cell_type(matrix: pd.DataFrame, metadata: pd.DataFrame, cell_column_name: str) -> pd.DataFrame:
    """
    This functions calculates the mean expression of each gene per group/cell type.

    Parameters
    ----------
    matrix: Normalized gene expression matrix (genes x barcodes).
    metadata: Index contains the barcode id and a single column named 'cell_type' indicating the group/cell type which
              the barcode belongs to.
    cell_column_name: Name of the column containing cell types

    Returns
    -------
    pd.DataFrame
        (genes x cell types) containing mean expression of a gene in a given cell type

    """
    core_logger.info("Scoring interactions: Calculating mean expression of each gene per group/cell type..")
    matrix = matrix.copy()
    out_dict = {}

    # Cell types present in Metadata
    cell_type_data = set(metadata[cell_column_name])

    for cell_type in tqdm(cell_type_data):
        # Obtain the barcode of the cells annotated with cell_type
        idx = metadata[cell_column_name] == cell_type
        cell_barcode = list(metadata.index[idx])

        # Calculate mean expression per cell type
        out_dict[cell_type] = matrix[cell_barcode].mean(axis=1)

    # Convert the dictionary to a dataframe
    matrix_mean_expr = pd.DataFrame.from_dict(out_dict)

    return matrix_mean_expr


def _geometric_mean(x):
    sub_values = list(x)
    sub_prod = np.prod(sub_values)
    geom = np.power(sub_prod, 1 / len(sub_values))
    return (geom)


def heteromer_geometric_expression_per_cell_type(
        matrix: pd.DataFrame,
        counts_data: str,
        genes: pd.DataFrame,
        complex_composition: pd.DataFrame,
        complex_expanded: pd.DataFrame,
        id2name: dict) -> pd.DataFrame:
    """
    Parameters
    ----------
    matrix: (genes x cell types) mean scaled expression matrix
    counts_data: Type of gene identifiers in the counts data: "ensembl", "gene_name", "hgnc_symbol"
    genes: CellphoneDB genes table
    complex_composition: CellphoneDB complex_composition table
    complex_expanded: CellphoneDB complex_expanded table
    id2name: A mapping between multidata_id and genes[counts_data]

    Returns
    -------
    pd.DataFrame
        (genes/complexes by cell type) in which the expression of a complex in a given cell type is
        a geometric mean of expressions of its constituents. Note that the only complexes included
        are the ones for which all the constituent genes are present in matrix
    """

    matrix = matrix.copy()

    # If the existing index does not contain genes[counts_data], replace it with index containing genes[counts_data] -
    if matrix.index.intersection(genes[counts_data]).empty:
        index_name = matrix.index.name
        matrix = matrix.reset_index()
        matrix[index_name].replace(to_replace=id2name, inplace=True)
        matrix.set_index(index_name, inplace=True)
        matrix.index.name = None

    # Subset the mean expression matrix to keep only the genes in CellphoneDB
    idx = [gene in list(genes[counts_data]) for gene in matrix.index]
    matrix = matrix.loc[idx]

    # Map complex name to its constituents/subunits
    complex_composition = pd.merge(complex_composition,
                                   complex_expanded[['complex_multidata_id', 'name']], on='complex_multidata_id')
    complex_composition = pd.merge(complex_composition,
                                   genes[[counts_data, 'protein_id']], left_on='protein_multidata_id',
                                   right_on='protein_id')
    d = complex_composition.groupby('name')[counts_data].apply(list).reset_index(name='subunits')
    complex_name_2_subunits = dict(zip(d['name'], d['subunits']))

    # Iterate over the complexes to calculate the geometric mean of expressions of their constituents
    # If any member of the subunit is not present in the means matrix then
    # the geometric mean expression is not calculated
    # The geometric mean is always lower than the arithmetic mean due to a compounding effect. Note that this is different
    # than in CellphoneDB methods themselves where the minimum expression across all parts of a heteromer is taken
    # as its overall expression.
    complex_geom_mean = dict()
    for complex_id in complex_name_2_subunits:

        # set used below to eliminate duplicate gene names
        # (in cases where the same gene appears more than once in genes table)
        subunits_set = set(complex_name_2_subunits[complex_id])

        # Remove nan values from heteromers
        subunits_list = [i for i in subunits_set if str(i) != 'nan']

        # Test if all the members of subunits_list are present in matrix
        # if true then calculate the geometric mean of the complex
        check_subunit = all([sub in matrix.index for sub in subunits_list])
        if check_subunit:
            complex_geom_mean[complex_id] = matrix.loc[subunits_list,].apply(_geometric_mean, axis=0)

    # Convert geometric mean dictionary to dataframe
    complex_geom_mean_df = pd.DataFrame.from_dict(complex_geom_mean,
                                                  orient='index')

    # Detect and remove genes that have the same name as a complex, e.g. OSMR, LIFR, IL2
    # Otherwise two rows assigned to the same gene would appear in final_df
    idx = [i not in complex_geom_mean_df.index for i in matrix.index]
    matrix = matrix.loc[idx]

    final_df = pd.concat([matrix, complex_geom_mean_df])

    return final_df


def scale_expression(matrix: pd.DataFrame, upper_range: int) -> pd.DataFrame:
    """
    Scale (up to upper_range) the expression of genes across all cell types in matrix

    Parameters
    ----------
    matrix: (genes x cell types) mean expression matrix
    upper_range: 0-upper_range is the range to which the expression of genes should be scaled

    Returns
    -------
    pd.DataFrame
        (genes x cell types) in which, for each gene in a given cell type, the expression was scaled to 0-upper_range
    """

    # Transpose matrix to apply scaling per row (i.e. scale across cell types)
    scaler = MinMaxScaler(feature_range=(0, upper_range), clip=True).fit(matrix.T)
    matrix_scaled = scaler.transform(matrix.T).T

    matrix_scaled = pd.DataFrame(matrix_scaled,
                                 index=matrix.index,
                                 columns=matrix.columns)
    return matrix_scaled


def _get_lr_scores(matrix, cpdb_set_all_lrs, separator, cell_type_pair) -> tuple:
    cell_type_A = cell_type_pair.split(separator)[0]
    cell_type_B = cell_type_pair.split(separator)[1]
    # Each cell in lr_outer is an arithmetic product:
    # gene in row's mean expression in cell_type_A and gene in column's mean expression in cell_type_B
    lr_outer = pd.DataFrame(np.outer(list(matrix[cell_type_A]),
                                     list(matrix[cell_type_B])),
                            index=matrix.index,
                            columns=matrix.index)
    # Round scores to 3 decimal places
    lr_outer = np.around(lr_outer, 3)
    df = lr_outer.stack().reset_index()
    lr_outer_long = pd.DataFrame({
        'interacting_pair': df.iloc[:, 0] + "_" + df.iloc[:, 1],
        "score": df.iloc[:, 2]})
    # Filtering interactions to only those in CellphoneDB
    idx_interactions = [i in cpdb_set_all_lrs for i in lr_outer_long['interacting_pair'].values]
    lr_outer_long_filtered = lr_outer_long.loc[idx_interactions]
    return (cell_type_pair, lr_outer_long_filtered)


def score_product(matrix: pd.DataFrame,
                  interactions: pd.DataFrame,
                  means: pd.DataFrame,
                  separator: str,
                  id2name: dict,
                  threads: int) -> dict:
    """
    For each interaction in CellphoneDB and a pair of cell types it calculates a score based on
    an arithmetic product of expressions of its participants in matrix

    Parameters
    ----------
    matrix: (genes/complexes by cell type) mean scaled expression matrix,
        where complex means are geometric means across their constituents(=subunits)
    counts_data: Type of gene identifiers in the counts data: "ensembl", "gene_name", "hgnc_symbol"
    genes: CellphoneDB genes table
    complex_expanded: CellphoneDB complex_expanded table
    interactions: CellphoneDB interactions table
    means: A copy of one of the result DataFrames of CellphoneDB that interaction scores will be populated into
    separator: separator character used in cell type pairs
    id2name: A mapping between multidata_id and genes[counts_data]
    threads: Number of threads to be used for parallel processing

    Returns
    -------
    dict
        Cell type pair identifier -> DataFrame containing the score annotated with each cell type and
        CellphoneDB interaction id.
    """

    core_logger.info("Scoring interactions: Calculating scores for all interactions and cell types..")
    matrix = matrix.copy()
    # All means in interaction_scores will be overwritten with scores
    interaction_scores = means

    interactions_df = interactions[['id_cp_interaction', 'multidata_1_id', 'multidata_2_id']].copy()
    interactions_df.replace(to_replace=id2name, inplace=True)
    interactions_df.rename(columns={'multidata_1_id': 'partner_a', 'multidata_2_id': 'partner_b'}, inplace=True)

    # Create lists of the interacting LR (and the reverse) to keep in lr_outer_long only those entries
    # that are described in the cpdb interaction file
    cpdb_list_a = list(interactions_df['partner_a']+'_'+interactions_df['partner_b'])
    cpdb_list_b = list(interactions_df['partner_b']+'_'+interactions_df['partner_a'])
    cpdb_set_all_lrs = set(cpdb_list_a + cpdb_list_b)

    cell_type_pairs = [c for c in interaction_scores.columns if separator in c]
    results = []

    with Pool(processes=threads) as pool:
        _get_lr_scores_thread = partial(_get_lr_scores, matrix, cpdb_set_all_lrs, separator)
        for tp in tqdm(pool.imap(_get_lr_scores_thread, cell_type_pairs),
                       total=len(cell_type_pairs)):
            results.append(tp)

    for ct_pair, lr_scores_filtered in results:
        interacting_pair2score = dict(zip(lr_scores_filtered['interacting_pair'], lr_scores_filtered['score']))
        interaction_scores[ct_pair] = [interacting_pair2score[id] for id in interaction_scores['interacting_pair']]

    return interaction_scores


# For each interaction in CellphoneDB and a pair of cell types it calculates a score based on
# an arithmetic product of expressions of its participants in counts matrix
def score_interactions_based_on_participant_expressions_product(
        cpdb_file_path: str,
        counts: pd.DataFrame,
        means: pd.DataFrame,
        separator: str,
        metadata: pd.DataFrame,
        threshold: float,
        cell_type_col_name: str,
        threads: int) -> pd.DataFrame:

    # Even if in the original analysis was done with counts_data = 'ensembl', we set it to 'gene_name' (= 'hgnc_symbol')
    # This is because id2name below relies on one-to-one mapping between multidata_id and genes[counts_data]
    counts_data = 'gene_name'

    # Get DB files
    interactions, genes, complex_composition, complex_expanded, _, _ = \
        db_utils.get_interactions_genes_complex(cpdb_file_path)

    #  Get mapping between multidata_id and genes[counts_data]
    id2name = dict(ChainMap(
        dict(zip(genes.protein_id, genes[counts_data])),
        dict(zip(complex_expanded.complex_multidata_id, complex_expanded.name))))

    # Step 1: Filter genes expressed in less than min_pct_cell of cells in a given cell type.
    cpdb_f = filter_genes_per_cell_type(matrix=counts,
                                        metadata=metadata,
                                        min_pct_cell=threshold,
                                        cell_column_name=cell_type_col_name)

    # Step 2: Calculate the gene's mean expression per cell type
    cpdb_fm = mean_expression_per_cell_type(matrix=cpdb_f,
                                            metadata=metadata,
                                            cell_column_name=cell_type_col_name)

    # Step 3: Calculate geometric expression mean per heteromer
    cpdb_fmsh = heteromer_geometric_expression_per_cell_type(matrix=cpdb_fm,
                                                             counts_data=counts_data,
                                                             genes=genes,
                                                             complex_composition=complex_composition,
                                                             complex_expanded=complex_expanded,
                                                             id2name=id2name)

    # Step 4: Scale the gene's mean expression across cell types.
    cpdb_fms = scale_expression(cpdb_fmsh,
                                upper_range=10)

    # Step 5: calculate the ligand-receptor score.
    interaction_scores = score_product(matrix=cpdb_fms,
                                       means=means,
                                       separator=separator,
                                       interactions=interactions,
                                       id2name=id2name,
                                       threads=threads)
    return interaction_scores
