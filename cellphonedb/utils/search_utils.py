import re
import time
from cellphonedb.utils import db_utils
from cellphonedb.utils.file_utils import dbg
import pandas as pd

SIMPLE_PFX="simple:"
COMPLEX_PFX = 'complex:'
ENS_PFX = "ENS"
INTERACTION_COLUMNS = ['interacting_pair', 'partner_a', 'partner_b', 'gene_a', 'gene_b']

def populate_proteins_for_complex(complex_name, complex_name2proteins, genes, complex_expanded, complex_composition):
    constituent_proteins = []
    if complex_name not in complex_name2proteins:
        complex_multidata_id = complex_expanded['complex_multidata_id'] \
            [complex_expanded[['name']].apply(lambda row: row.astype(str).eq(complex_name).any(), axis=1)].to_list()[0]
        for protein_multidata_id in complex_composition['protein_multidata_id'] \
                [complex_composition['complex_multidata_id'] == complex_multidata_id].to_list():
            constituent_proteins.append(genes['name'][genes['protein_multidata_id'] == protein_multidata_id].to_list()[0])
        complex_name2proteins[complex_name] = constituent_proteins

def search(query_str: str = "",
           cpdb_file_path: str = None)->(list, map):
    """
    Searches CellphoneDB interactions for genes/proteins/complexes in query_str

    Parameters
    ----------
    query_str: str
        A comma- or space-separated list of Ensembl ID (e.g. ENSG00000165029), Gene name (e.g. ABCA1), \
        UniProt ID (e.g. KLRG2_HUMAN), UniProt Accession (e.g. A4D1S0) or Complex name \
        (e.g. 12oxoLeukotrieneB4_byPTGR1)
    cpdb_file_path: str
        CellphoneDB database file path

    Returns
    -------
    Tuple
        - A list of sub-lists, where each sub-list represents details of a single interaction
        - For all complexes that participate in the interactions found, a mapping to their constituent proteins
    """
    start = time.time()
    results = []
    interactions, genes, complex_composition, complex_expanded, gene_synonym2gene_name = \
        db_utils.get_interactions_genes_complex(cpdb_file_path)


    complex_name2proteins = {}
    # Assemble a list of multidata_ids to search interactions DF with
    multidata_ids = []
    for token in re.split(',\s*| ', query_str):

        if token in gene_synonym2gene_name:
            # Map any gene synonyms not in gene_input to gene names in gene_input
            token = gene_synonym2gene_name[token]

        complex_multidata_ids = []
        # Attempt to find token in genes (N.B. genes contains protein information also)
        gene_protein_data_list = genes['protein_multidata_id'] \
            [genes[['ensembl','gene_name','name','protein_name']].apply(lambda row: row.astype(str).eq(token).any(), axis=1)].to_list()
        if (len(gene_protein_data_list) > 0):
            multidata_ids += gene_protein_data_list
            for protein_multidata_id in gene_protein_data_list:
                complex_multidata_ids = \
                    complex_composition['complex_multidata_id'] \
                        [complex_composition['protein_multidata_id'] == protein_multidata_id].to_list()
                multidata_ids += complex_multidata_ids
        else:
            # No match in genes - attempt to find token in complex_expanded
            complex_multidata_ids += complex_expanded['complex_multidata_id'] \
                [complex_expanded[['name']].apply(lambda row: row.astype(str).eq(token).any(), axis=1)].to_list()
            multidata_ids += complex_multidata_ids
    # Now search for all multidata_ids in interactions
    duration = time.time() - start
    dbg("Output for query '{}':".format(query_str))
    # Output header
    results.append(['CellphoneDB interaction ', 'Partner A', 'Partner B', 'Gene name A', 'Gene name B', ' Ensembl ID A', 'Ensembl ID B', 'Annotation strategy', 'Curator','Source','Is PPI'])
    for multidata_id in multidata_ids:
         interactions_data_list = interactions[[ \
             'id_cp_interaction','multidata_1_id', 'multidata_2_id', \
             'name_1', 'name_2','is_complex_1', 'is_complex_2', 'annotation_strategy','curator','source','is_ppi']] \
            [interactions[['multidata_1_id', 'multidata_2_id']].apply(lambda row: row.astype(int).eq(multidata_id).any(),
                                                                     axis=1)].values.tolist()
         for interaction in interactions_data_list:
             output_row = [interaction[0]]
             if interaction[5] == False:
                 data_1 = generate_output(interaction[1], genes)
                 prefix_1 = SIMPLE_PFX
             else:
                 data_1 = ['','']
                 prefix_1 = COMPLEX_PFX
                 populate_proteins_for_complex(interaction[3], complex_name2proteins,
                                               genes, complex_expanded, complex_composition)
             if interaction[6] == False:
                 data_2 = generate_output(interaction[2], genes)
                 prefix_2 = SIMPLE_PFX
             else:
                 data_2 = ['', '']
                 prefix_2 = COMPLEX_PFX
                 populate_proteins_for_complex(interaction[4], complex_name2proteins,
                                               genes, complex_expanded, complex_composition)
             output_row += [
                 prefix_1+interaction[3],
                 prefix_2+interaction[4],
                 data_1[0],
                 data_2[0],
                 data_1[1],
                 data_2[1],
                 interaction[7],
                 interaction[8],
                 interaction[9],
                 interaction[10]
             ]
             results.append(output_row)
    dbg("Total search time: " + str(round(duration, 2)) + "s")
    return results, complex_name2proteins


def generate_output(multidata_id, genes):
    protein_data_list = genes[['gene_name','ensembl']][genes['protein_multidata_id'] == multidata_id].values.tolist()
    if len(protein_data_list) == 0:
        return None
    # Expect only a single result
    return protein_data_list[0]

def get_uniprot_url(uniprot_accessions) -> str:
    url_prefix = "https://www.uniprot.org/uniprotkb?query="
    url_query = ""
    # (accession:Q15825)%20OR%20(accession:P43681)%20OR%20(accession:Q05901)
    for acc in uniprot_accessions:
        if url_query != "":
            url_query += "%20OR%20"
        url_query += "(accession:{})".format(acc)
    return url_prefix + url_query

def get_html_table(data, complex_name2proteins) -> str:
    """
    Parameters
    ----------
    data: list
        A list of sub-lists, where each sub-list represents details of a single interaction
    complex_name2proteins: map
        For all complexes that participate in the interactions found, a mapping to their constituent proteins
    Returns
    -------
    str
        A string containing html representation of interactions in data and complex_name2proteins. \
        The mapping in the latter is used to populate tooltips showing proteins that are a part of a complex.
    """
    html = "<table id=\"cpdb_search_results\" class=\"display compact\">"
    first_row = True
    for row in data:
        if first_row:
            html += "<thead>"
        html += "<tr>"
        for field in row:
            # NB. the two last elements of row contain the mouseover text containing the complex's constituent proteins
            if first_row:
                html += "<th style=\"text-align:left\">{}</th>".format(field)
            else:
                if field.startswith(COMPLEX_PFX):
                    name = field.split(":")[1]
                    complex_mouseover = "Contains proteins: " + ', '.join(complex_name2proteins[name])
                    multi_protein_uniprot_url = get_uniprot_url(complex_name2proteins[name])
                    html += "<td style=\"text-align:left\"><a class=\"teal-text\" target=\"_blank\" title=\"{}\" href=\"{}\">{}</a></td>".format(complex_mouseover, multi_protein_uniprot_url, name)
                elif field.startswith(SIMPLE_PFX):
                    name = field.split(":")[1]
                    html += "<td style=\"text-align:left\"><a class=\"teal-text\" target=\"_blank\" href=\"https://www.uniprot.org/uniprotkb/{}/entry\">{}</a></td>" \
                        .format(name, name)
                elif field.startswith(ENS_PFX):
                    html += "<td style=\"text-align:left\"><a class=\"teal-text\" target=\"_blank\" href=\"https://www.ensembl.org/id/{}\">{}</a></td>" \
                        .format(field, field)
                else:
                    html += "<td style=\"text-align:left\">{}</td>".format(field)
        if first_row:
            html += "</thead><tbody>"
        html += "</tr>"
        first_row = False
    html += "</tbody></table>"
    return html

def autocomplete_query(genes: pd.DataFrame, interactions: pd.DataFrame, partial_element: str) -> pd.DataFrame:
    values = _partial_filter(genes, 'ensembl', partial_element)
    by_protein_name = _partial_filter(genes, 'protein_name', partial_element)
    by_gene_name = _partial_filter(genes, 'gene_name', partial_element)
    with_hgnc_symbol = genes.dropna(subset=['hgnc_symbol'])
    by_hgnc_symbol = _partial_filter(with_hgnc_symbol, 'hgnc_symbol', partial_element)
    by_name_1 = _partial_filter(interactions, 'name_1', partial_element)
    by_name_2 = _partial_filter(interactions, 'name_2', partial_element)

    values = values.append(by_protein_name, ignore_index=True)
    values = values.append(by_gene_name, ignore_index=True)
    values = values.append(by_hgnc_symbol, ignore_index=True)
    values = values.append(by_name_1, ignore_index=True)
    values = values.append(by_name_2, ignore_index=True)
    result = pd.DataFrame(data=values, columns=['value']).drop_duplicates()
    return result

def _partial_filter(input_data, name, partial_element):
    matching = input_data[input_data[name].str.contains(partial_element, flags=re.IGNORECASE)][name]
    return matching

def return_all_identifiers(genes: pd.DataFrame, interactions: pd.DataFrame) -> pd.DataFrame:
    values = genes['ensembl']
    for col in ['gene_name','protein_name','hgnc_symbol']:
        values = pd.concat([values, genes[col]], ignore_index=True)
    for col in ['name_1', 'name_2']:
        values = pd.concat([values, interactions[col]], ignore_index=True)
    result = pd.DataFrame(data=values, columns=['value']).drop_duplicates()
    return result

def search_analysis_results(
        query_cell_types_1: list = None,
        query_cell_types_2: list = None,
        query_genes: list = None,
        query_interactions: list = None,
        significant_means: pd.DataFrame = None,
        deconvoluted: pd.DataFrame = None,
        separator: str = "|",
        long_format: bool = False
) -> pd.DataFrame:
    """
    Searches results of either statistical or DEG analysis for relevant interactions
    matching any of:
        1. A gene in query_genes
        2. A complex containing a gene in query_genes
        3. An interaction name in query_interactions (e.g. 12oxoLeukotrieneB4_byPTGR1)
    where at least one pair of cell types containing one cell type from query_cell_types_1
    and one cell type from query_cell_types_2 has a significant mean.
    NB. If all of query_cell_types_1, query_cell_types_2, query_genes and query_interactions are set to None,
    then all relevant interactions are returned.

    Parameters
    ----------
    query_cell_types_1: list
        A list of cell types
    query_cell_types_2: list
        A list of cell types
    query_genes: list
        A list of gene names
    query_interactions: list
        A list of interactions
    significant_means: pd.DataFrame
    deconvoluted: pd.DataFrame
        Files output by either (by the same) statistical or DEG analysis
    separator: str
        Separator used in cell type pair column names in significant_means dataFrame
    long_format: bool
        Return the search result DataFrame in long format (while dropping rows with NaN in significant_mean column)
    Returns
    -------
    pd.DataFrame
        Relevant interactions from significant_means that match query criteria
    """
    if significant_means is None or deconvoluted is None:
        print("ERROR: Both significant_means and deconvoluted dataframes need to be provided")
        return

    if query_cell_types_1 is None or query_cell_types_1 is None:
        print ("ERROR: Both query_cell_types_1 and query_cell_types_1 need to be provided. " + \
               "If you wish to search for all cell type combinations, set them both to \"All\"")

    # Collect all combinations of cell types (disregarding the order) from query_cell_types_1 and query_cell_types_2
    if query_cell_types_1 == "All" or query_cell_types_2 == "All":
        cols_filter = significant_means.filter(regex="\{}".format(separator)).columns
        all_cts = set([])
        for ct_pair in [i.split(separator) for i in cols_filter.tolist()]:
            all_cts |= set(ct_pair)
        all_cell_types = list(all_cts)
        if query_cell_types_1 == "All":
            query_cell_types_1 = all_cell_types
        if query_cell_types_2 == "All":
            query_cell_types_2 = all_cell_types
    cell_type_pairs = []
    for ct in query_cell_types_1:
        for ct1 in query_cell_types_2:
            cell_type_pairs += ["{}{}{}".format(ct, separator, ct1), "{}{}{}".format(ct1, separator, ct)]
    cols_filter = significant_means.columns[significant_means.columns.isin(cell_type_pairs)]

    # Collect all interactions from query_genes and query_interactions
    interactions = set([])
    if query_genes:
            interactions = interactions.union(frozenset(deconvoluted[deconvoluted['gene_name'].isin(query_genes)]['id_cp_interaction'].tolist()))
    if query_interactions:
        interactions = interactions.union(frozenset(significant_means[significant_means['interacting_pair'].isin(query_interactions)]['id_cp_interaction'].tolist()))

    if interactions:
        result_df = significant_means[significant_means['id_cp_interaction'].isin(interactions)]
    else:
        result_df = significant_means

    # Filter out cell_type_pairs/columns in cols_filter for which no interaction in interactions set is significant
    cols_filter = cols_filter[result_df[cols_filter].notna().any(axis=0)]

    # Filter out interactions which are not significant in any cell_type_pair/column in cols_filter
    result_df = result_df[result_df[cols_filter].notna().any(axis=1)]
    # Select display columns
    result_df = result_df[INTERACTION_COLUMNS + cols_filter.tolist()]

    if long_format:
        # Convert the results DataFrame from (default) wide to long format
        result_df = pd.melt(result_df,
                                 id_vars=result_df.columns[0:len(INTERACTION_COLUMNS)],
                                 value_vars=result_df.columns[len(INTERACTION_COLUMNS):],
                                 value_name='significant_mean',
                                 var_name='interacting_cells') \
            .dropna(subset=['significant_mean'])

    return result_df