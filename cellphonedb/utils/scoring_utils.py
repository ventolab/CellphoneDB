import pandas as pd
import numpy as np
from cellphonedb.utils import db_utils
from scipy.stats.mstats import gmean
from sklearn.preprocessing import MinMaxScaler
from itertools import combinations_with_replacement

from functools import partial
from multiprocessing.pool import Pool
from tqdm.std import tqdm

import time

def filter_genes_cluster(
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
    metadata: Index contains the barcode id and a single column named 'cell_type' indicating the group/cell type which the barcode belongs to.
    min_pct_cell : Percentage of cells required to express a given gene.
    cell_column_name: Name of column containing the cluster name

    Returns
    -------
    pd.DataFrame
        matrix with expression of a gene for all cells set to 0 - if that gene is expressed in less than min_pct_cell

    """
    matrix = matrix.copy()

    # Cell types present in metadata
    cell_type_data = set(metadata[cell_column_name])

    for cell_type in cell_type_data:
        # Obtain the barcode of the cells annotated under cell_type
        idx = metadata[cell_column_name] == cell_type
        cell_barcode = list(metadata.index[idx])

        # Calculate percentage of cells expressing the gene (expression != 0)
        gene_expr_pct = (matrix[cell_barcode] != 0).sum(axis=1) / len(cell_barcode)

        # List of genes lowly expressed
        gene_lowly_expr = list(gene_expr_pct.index[gene_expr_pct < min_pct_cell])

        # Set to zero those genes that are expressed in a given cell type below the 
        # user defined min_pct_cell
        matrix.loc[gene_lowly_expr, cell_barcode] = 0

    # Return filtered matrix
    return matrix


def mean_expression_cluster(matrix: pd.DataFrame, metadata: pd.DataFrame, cell_column_name: str) -> pd.DataFrame:
    """
    This functions calculates the mean expression of each gene per group/cell type.

    Parameters
    ----------
    matrix: Normalized gene expression matrix (genes x barcodes).
    metadata: Index contains the barcode id and a single column named 'cell_type' indicating the group/cell type which the barcode belongs to.
    cell_column_name: Name of column containing the cluster name

    Returns
    -------
    pd.DataFrame
        (genes x cell types) containing mean expression of a gene in a given cell type

    """
    matrix = matrix.copy()
    out_dict = {}

    # Cell types present in Metadata
    cell_type_data = set(metadata[cell_column_name])

    for cell_type in cell_type_data:
        # Obtain the barcode of the cells annotated with cell_type
        idx = metadata[cell_column_name] == cell_type
        cell_barcode = list(metadata.index[idx])

        # Calculate mean expression per cluster
        out_dict[cell_type] = matrix[cell_barcode].mean(axis=1)

    # Convert the dictionary to a dataframe
    matrix_mean_expr = pd.DataFrame.from_dict(out_dict)

    return matrix_mean_expr


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
    scaler = MinMaxScaler(feature_range=(0, upper_range)).fit(matrix.T)
    matrix_scaled = scaler.transform(matrix.T).T

    matrix_scaled = pd.DataFrame(matrix_scaled,
                                 index=matrix.index,
                                 columns=matrix.columns)
    return matrix_scaled

def _geometric_mean(x):
    sub_values = list(x)
    sub_prod = np.prod(sub_values)
    geom = np.power(sub_prod, 1 / len(sub_values))
    return (geom)

def heteromer_geometric_expression(matrix: pd.DataFrame, cpdb_file_path: str) -> pd.DataFrame:
    """
    Parameters
    ----------
    matrix: (genes x cell types) mean scaled expression matrix
    cpdb_file_path: A full path of CellphoneDB database file

    Returns
    -------
    pd.DataFrame
        (genes/complexes by cell type) in which the expression of a complex in a given cell type is
        a geometric mean of expressions of its constituents. Note that the only complexes included
        are the ones for which all the constituent genes are in matrix
    """
    interactions, genes, complex_composition, complex_expanded, gene_synonym2gene_name = \
        db_utils.get_interactions_genes_complex(cpdb_file_path)

    matrix = matrix.copy()

    # Subset the mean expression matrix to keep only the genes in CellphoneDB
    print(matrix.shape)
    idx = [gene in list(genes['gene_name']) for gene in matrix.index]
    matrix = matrix.loc[idx]
    print(matrix.shape)

    # Map complex name to its constituents/subunits
    complex_composition = pd.merge(complex_composition,
                                   complex_expanded[['complex_multidata_id', 'name']], on='complex_multidata_id')
    complex_composition = pd.merge(complex_composition, \
                                   genes[['gene_name', 'protein_id']], left_on='protein_multidata_id',
                                   right_on='protein_id')
    d = complex_composition.groupby('name')['gene_name'].apply(list).reset_index(name='subunits')
    complex_name_2_subunits = dict(zip(d['name'], d['subunits']))

    # Iterate over the complexes to calculate the geometric mean of expressions of their constituents
    # If any member of the subunit is not present in the means matrix
    # the geometric mean expression is not calculated
    complex_geom_mean = dict()
    for complex_id in complex_name_2_subunits:

        subunits_list = complex_name_2_subunits[complex_id]

        # Remove nan values from heteromers
        subunits_list = [i for i in subunits_list if str(i) != 'nan']

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

def _get_lr_outer_long(matrix, cell_type_tuple):
    cell_A = cell_type_tuple[0]
    cell_B = cell_type_tuple[1]
    # outer product
    # rows: first element outer function
    # cols: second element outer function
    lr_outer = pd.DataFrame(np.outer(list(matrix[cell_A]),
                                     list(matrix[cell_B])),
                            index=matrix.index,
                            columns=matrix.index)

    upper_idx = np.triu(np.ones(lr_outer.shape)).astype(bool)
    upper_tri = lr_outer.where(upper_idx)

    lower_idx = np.tril(np.ones(lr_outer.shape)).astype(bool)
    lower_tri = lr_outer.where(lower_idx)

    # Upper triangle (mtx_up) represents communication between cell B to cell A
    # Lower triangle (mtx_df) represents communication between cell A to cell B
    mtx_up = upper_tri.stack().reset_index()
    mtx_dn = lower_tri.stack().reset_index()

    mtx_up.columns = [cell_A, cell_B, 'Score']
    mtx_dn.columns = [cell_A, cell_B, 'Score']

    lr_outer_long = pd.concat([mtx_up, mtx_dn])
    return (cell_type_tuple, lr_outer_long)

def _get_score(idx_interactions, interactions_df, cell_type_tuple, lr_outer_long):
    cell_A = cell_type_tuple[0]
    cell_B = cell_type_tuple[1]
    lr_outer_long_filt = lr_outer_long.loc[idx_interactions].copy()

    # Add interaction id
    lr_outer_sorted = ['-'.join(sorted([a, b])) for a, b in
                       zip(lr_outer_long_filt.iloc[:, 0], lr_outer_long_filt.iloc[:, 1])]
    cpdb_sorted = ['-'.join(sorted([a, b])) for a, b in zip(interactions_df['partner_a'], interactions_df['partner_b'])]

    dict_cpdb_id = dict(zip(cpdb_sorted, interactions_df['id_cp_interaction']))
    lr_outer_long_filt['id_cp_interaction'] = [dict_cpdb_id[i] for i in lr_outer_sorted]

    # assign cell names
    # first column is the sender
    # second column the receiver
    # third column the score
    lr_outer_long_filt.columns = [cell_A, cell_B, 'Score', 'id_cp_interaction']
    lr_outer_long_filt = lr_outer_long_filt.reset_index(drop=True)

    return ('|'.join(sorted([cell_A, cell_B])), lr_outer_long_filt)


def score_product(matrix: pd.DataFrame, cpdb_file_path: str, threads: int) -> dict:
    """
    TODO
    Parameters
    ----------
    matrix: (genes/complexes by cell type) mean scaled expression matrix,
        where complex means are geometric means across their constituents/subunits
    cpdb_file_path: A full path of CellphoneDB database file
    threads: Number of threads to be used for parallel processing

    Returns
    -------
    dict
        TODO

    """
    matrix = matrix.copy()

    interactions, genes, complex_composition, complex_expanded, gene_synonym2gene_name = \
        db_utils.get_interactions_genes_complex(cpdb_file_path)

    id2name = dict(zip(genes.protein_id, genes.gene_name))
    id2name = id2name | dict(zip(complex_expanded.complex_multidata_id, complex_expanded.name))
    interactions_df = interactions[['id_cp_interaction', 'multidata_1_id', 'multidata_2_id']].copy()
    interactions_df.replace(to_replace=id2name, inplace=True)
    interactions_df.rename(columns={'multidata_1_id': 'partner_a', 'multidata_2_id': 'partner_b'}, inplace=True)

    # Create lists of the interacting LR (and the reverse) to keep in the lr_outer_long only those entries
    # that are described in the cpdb interaction file
    cpdb_list_a = list(interactions_df['partner_a']+'|'+interactions_df['partner_b'])
    cpdb_list_b = list(interactions_df['partner_b']+'|'+interactions_df['partner_a'])
    cpdb_list_all = cpdb_list_a + cpdb_list_b

    # Calculate all cell type combinations
    # Initialize score dictionary
    combinations_cell_types = list(combinations_with_replacement(matrix.columns, 2))
    score_dict = {}

    t0 = time.time()
    with Pool(processes=threads) as pool:
        _get_lr_outer_long_thread = partial(_get_lr_outer_long, matrix)
        results = tqdm(pool.imap(_get_lr_outer_long_thread, combinations_cell_types),
                       total=len(combinations_cell_types))
        tuple(results)  # fetch the lazy results
    print(time.time() - t0, " got lr_outer_longs")

    cnt = 0
    for cell_type_tuple, lr_outer_long in results:
        if cnt == 0:
            print(cell_type_tuple + " : " + lr_outer_long)
            lr_list = list(lr_outer_long.iloc[:, 0] + '|' + lr_outer_long.iloc[:, 1])
            idx_interactions = [i in cpdb_list_all for i in lr_list]
        tp = _get_score(idx_interactions, interactions_df, cell_type_tuple, lr_outer_long)
        score_dict[tp[0]] = tp[1]
        cnt += 1


    return score_dict
