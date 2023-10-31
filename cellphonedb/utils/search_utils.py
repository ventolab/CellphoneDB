import re
import time
from cellphonedb.utils import db_utils
from cellphonedb.utils.file_utils import dbg
import pandas as pd

SIMPLE_PFX = "simple:"
COMPLEX_PFX = 'complex:'
ENS_PFX = "ENS"
INTERACTION_COLUMNS = ['interacting_pair', 'partner_a', 'partner_b', 'gene_a', 'gene_b', 'directionality', 'classification']
EXTERNAL_RESOURCE2URI = {'Reactome reaction': 'https://reactome.org/content/detail',
                         'Reactome complex': 'https://reactome.org/content/detail',
                         'ComplexPortal complex': 'https://www.ebi.ac.uk/complexportal/complex',
                         'Rhea reaction': 'https://www.rhea-db.org/rhea'}
SIDENAV_PROPERTY_STYLE = "style=\"padding-left: 60px; font-size: 14px; margin: 20px 0px !important; \""
SIDENAV_A_STYLE = "style=\"padding-left: 60px; font-size: 14px; margin: 20px 0px !important;\""


def populate_proteins_for_complex(complex_name, complex_name2proteins, genes, complex_expanded, complex_composition):
    constituent_proteins = []
    if complex_name not in complex_name2proteins:
        complex_multidata_id = complex_expanded['complex_multidata_id'][
            complex_expanded[['name']].apply(lambda row: row.astype(str).eq(complex_name).any(), axis=1)
        ].to_list()[0]
        for protein_multidata_id in \
                complex_composition['protein_multidata_id'][complex_composition['complex_multidata_id'] ==
                                                            complex_multidata_id].to_list():
            constituent_proteins.append(genes['name'][genes['protein_multidata_id'] == protein_multidata_id].to_list()[0])
        complex_name2proteins[complex_name] = constituent_proteins


def assemble_multidata_ids_for_search(
        query_str: str,
        genes: pd.DataFrame,
        complex_expanded: pd.DataFrame,
        complex_composition: pd.DataFrame,
        gene_synonym2gene_name: dict) -> list:
    multidata_ids = []
    for token in re.split(',\\s*| ', query_str):

        if token in gene_synonym2gene_name:
            # Map any gene synonyms not in gene_input to gene names in gene_input
            token = gene_synonym2gene_name[token]

        complex_multidata_ids = []
        # Attempt to find token in genes (N.B. genes contains protein information also)
        gene_protein_data_list = \
            genes['protein_multidata_id'][
                genes[['ensembl', 'gene_name', 'name', 'protein_name']]
                .apply(lambda row: row.astype(str).eq(token).any(), axis=1)
            ].to_list()
        if (len(gene_protein_data_list) > 0):
            multidata_ids += gene_protein_data_list
            for protein_multidata_id in gene_protein_data_list:
                complex_multidata_ids = \
                    complex_composition['complex_multidata_id'][complex_composition['protein_multidata_id']
                                                                == protein_multidata_id].to_list()
                multidata_ids += complex_multidata_ids
        else:
            # No match in genes - attempt to find token in complex_expanded
            complex_multidata_ids += \
                complex_expanded['complex_multidata_id'][
                    complex_expanded[['name']]
                    .apply(lambda row: row.astype(str).eq(token).any(), axis=1)
                ].to_list()
            multidata_ids += complex_multidata_ids
    return multidata_ids


def search(query_str: str = "",
           cpdb_file_path: str = None) -> (list, map, map, map, map):
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
    interactions, genes, complex_composition, complex_expanded, gene_synonym2gene_name, _ = \
        db_utils.get_interactions_genes_complex(cpdb_file_path)

    interaction_columns = []
    # Cater for DB version-dependent column names
    if 'directionality' in interactions.columns:
        interaction_columns = ['curator', 'source', 'is_ppi', 'directionality', 'classification']

    protein2Info, complex2Info, resource2Complex2Acc, proteinAcc2Name = \
        db_utils.get_protein_and_complex_data_for_web(cpdb_file_path)

    complex_name2proteins = {}
    # Assemble a list of multidata_ids to search interactions DF with
    multidata_ids = assemble_multidata_ids_for_search(
        query_str, genes, complex_expanded, complex_composition, gene_synonym2gene_name)

    # Now search for all multidata_ids in interactions
    duration = time.time() - start
    dbg("Output for query '{}':".format(query_str))
    # Output header
    results.append(['CellphoneDB interaction ', 'Partner A', 'Partner B', 'Gene name A', 'Gene name B', ] +
                   [' Ensembl ID A', 'Ensembl ID B', 'Annotation strategy', 'Curator', 'Source', 'Is PPI'] +
                   ['Directionality', 'Classification'])
    for multidata_id in multidata_ids:
        interactions_data_list = interactions[[
             'id_cp_interaction', 'multidata_1_id', 'multidata_2_id',
             'name_1', 'name_2', 'is_complex_1', 'is_complex_2', 'annotation_strategy'] + interaction_columns][
                interactions[['multidata_1_id', 'multidata_2_id']]
                .apply(lambda row: row.astype(int).eq(multidata_id).any(), axis=1)
        ].values.tolist()
        for interaction in interactions_data_list:
            output_row = [interaction[0]]
            if not interaction[5]:
                data_1 = generate_output(interaction[1], genes)
                prefix_1 = SIMPLE_PFX
            else:
                data_1 = ['', '']
                prefix_1 = COMPLEX_PFX
                populate_proteins_for_complex(interaction[3], complex_name2proteins,
                                              genes, complex_expanded, complex_composition)
            if not interaction[6]:
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
                interaction[7]
                ]
            # Cater for DB version-dependent column names
            if len(interaction) > 8:
                output_row += [
                    interaction[8],
                    interaction[9],
                    interaction[10],
                    interaction[11],
                    interaction[12]
                ]
            results.append(output_row)
    dbg("Total search time: " + str(round(duration, 2)) + "s")
    return results, complex_name2proteins, protein2Info, complex2Info, resource2Complex2Acc, proteinAcc2Name


def generate_output(multidata_id, genes):
    protein_data_list = genes[['gene_name', 'ensembl']][genes['protein_multidata_id'] == multidata_id].values.tolist()
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


def get_sidenav_html(is_complex: bool, bioentity_name: str,
                     complex_name2proteins, protein2Info, complex2Info, resource2Complex2Acc, proteinAcc2Name) -> str:
    if is_complex:
        constituent_proteins = ', '.join(complex_name2proteins[bioentity_name])
        multi_protein_uniprot_url = get_uniprot_url(complex_name2proteins[bioentity_name])
        external_resource_links = ""
        for res in EXTERNAL_RESOURCE2URI:
            if res in resource2Complex2Acc:
                resource_label = res.split(" ")[0]
                if bioentity_name in resource2Complex2Acc[res]:
                    acc = resource2Complex2Acc[res][bioentity_name]
                    res_lookup_id = acc.replace("RHEA:", "")
                    external_resource_links += "<a class=\"teal-text\" href=\"{}/{}\" target=\"blank\" {}>{} {}</a><br>" \
                        .format(EXTERNAL_RESOURCE2URI[res], res_lookup_id, SIDENAV_A_STYLE, resource_label, acc)
        complexInformation = ""
        for item in complex2Info[bioentity_name]:
            complexInformation += "<a {}>{}</a><br> ".format(SIDENAV_PROPERTY_STYLE, item)
        html = ("<li><a class=\"subheader\">Complex Information</a></li>" +
                "<a {}><b>{}</b></a>" +
                "<li><a class=\"subheader\">Members</a></li>" +
                "<a class=\"teal-text\" href=\"{}\" target=\"blank\" {}>{} (see in UniProt)</a>" +
                "<li><a class=\"subheader\">Properties</li>" +
                "{}" +
                "<li><a class=\"subheader\">Cross References</li>" +
                "{}").format(SIDENAV_PROPERTY_STYLE, bioentity_name,
                             multi_protein_uniprot_url, SIDENAV_A_STYLE, constituent_proteins,
                             complexInformation, external_resource_links)
    else:
        proteinInformation = ""
        for item in protein2Info[bioentity_name]:
            proteinInformation += "<a {}>{}</a><br>".format(SIDENAV_PROPERTY_STYLE, item)
        if bioentity_name in proteinAcc2Name:
            proteinName = proteinAcc2Name[bioentity_name]
        else:
            proteinName = ""
        html = ("<li><a class=\"subheader\">Protein Information</a></li>" +
                "<a href=\"https://www.uniprot.org/uniprotkb/{}/entry\" class=\"teal-text\" target=\"blank\" " +
                "{}><b>{}</b> (See in UniProt)</a><br><a {}>{}</a>" +
                "<li><a class=\"subheader\">Properties</li>" +
                "{}").format(bioentity_name, SIDENAV_A_STYLE, bioentity_name, SIDENAV_A_STYLE, proteinName,
                             proteinInformation)
    return html


def get_html_table(data, complex_name2proteins,
                   protein2Info, complex2Info, resource2Complex2Acc, proteinAcc2Name) -> str:
    """
    Parameters
    ----------
    data: list
        A list of sub-lists, where each sub-list represents details of a single interaction
    complex_name2proteins: map
        For all complexes that participate in the interactions found, a mapping to their constituent proteins
    protein2Info: map
        A mapping between a protein and the list of its details (e.g. 'transmembrane') to be displayed in
        the sidebar when the user clicks in a protein accession in the search results table
    protein2Info: map
        A mapping between a protein and the list of its details (e.g. 'transmembrane') to be displayed in
        the sidebar when the user clicks in a protein accession in the search results table
    resource2Complex2Acc: map
        A mapping between a external resource name (a key in EXTERNAL_RESOURCE2URI) and
        a map between a complex name and its corresponding accession in that external resource

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
        for field in [str(x) for x in row]:
            # NB. the two last elements of row contain the mouseover text containing the complex's constituent proteins
            if first_row:
                html += "<th style=\"text-align:left\">{}</th>".format(field)
            else:
                if field.startswith(COMPLEX_PFX):
                    name = ":".join(field.split(":")[1:])
                    constituent_proteins = ', '.join(complex_name2proteins[name])
                    complex_mouseover = "Contains proteins: {}".format(constituent_proteins)
                    inner_html = get_sidenav_html(field.startswith(COMPLEX_PFX), name,
                                                  complex_name2proteins, protein2Info,
                                                  complex2Info, resource2Complex2Acc, proteinAcc2Name)
                    html += ("<td style=\"text-align:left\"><a class=\"teal-text sidenav-trigger\" " +
                             "data-target='sidenav_{}' title=\"{}\" href=\"#\">{}</a>" +
                             "<ul id=\"sidenav_{}\" class=\"sidenav fixed\" style=\"width:410px\">" +
                             inner_html + "</ul></td>").format(name, complex_mouseover, name, name)
                elif field.startswith(SIMPLE_PFX):
                    name = field.split(":")[1]
                    inner_html = get_sidenav_html(False, name,
                                                  complex_name2proteins, protein2Info, complex2Info,
                                                  resource2Complex2Acc, proteinAcc2Name)
                    html += ("<td style=\"text-align:left\"><a class=\"teal-text sidenav-trigger\" data-target='sidenav_{}' " +
                             "title=\"{}\" href=\"#\">{}</a>" +
                             "<ul id=\"sidenav_{}\" class=\"sidenav fixed\" style=\"width:410px\">" +
                             inner_html + "</ul></td>").format(name, name, name, name)
                elif field.startswith(ENS_PFX):
                    html += ("<td style=\"text-align:left\"><a class=\"teal-text\" target=\"_blank\" " +
                             "href=\"https://www.ensembl.org/id/{}\">{}</a></td>") \
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
    for col in ['gene_name', 'protein_name', 'hgnc_symbol']:
        values = pd.concat([values, genes[col]], ignore_index=True)
    for col in ['name_1', 'name_2']:
        values = pd.concat([values, interactions[col]], ignore_index=True)
    result = pd.DataFrame(data=values, columns=['value']).drop_duplicates()
    return result


def collect_celltype_pairs(
        significant_means: pd.DataFrame,
        query_cell_types_1: list,
        query_cell_types_2: list,
        separator: str) -> list:
    if query_cell_types_1 is None or query_cell_types_2 is None:
        cols_filter = significant_means.filter(regex="\\{}".format(separator)).columns
        all_cts = set([])
        for ct_pair in [i.split(separator) for i in cols_filter.tolist()]:
            all_cts |= set(ct_pair)
        all_cell_types = list(all_cts)
        if query_cell_types_1 is None:
            query_cell_types_1 = all_cell_types
        if query_cell_types_2 is None:
            query_cell_types_2 = all_cell_types
    cell_type_pairs = []
    for ct in query_cell_types_1:
        for ct1 in query_cell_types_2:
            cell_type_pairs += ["{}{}{}".format(ct, separator, ct1), "{}{}{}".format(ct1, separator, ct)]
    return cell_type_pairs


def search_analysis_results(
        query_cell_types_1: list = None,
        query_cell_types_2: list = None,
        query_genes: list = None,
        query_interactions: list = None,
        query_classifications: list = None,
        query_minimum_score: int = None,
        significant_means: pd.DataFrame = None,
        deconvoluted: pd.DataFrame = None,
        interaction_scores=None,
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
    query_classifications: list
        A list of query classifications
    query_minimum_score: int
        Find all interactions with at least minimum_score across the selected cell types
    significant_means: pd.DataFrame
    deconvoluted: pd.DataFrame
    interaction_scores: pd.DataFrame
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

    # Collect all combinations of cell types (disregarding the order) from query_cell_types_1 and query_cell_types_2
    cell_type_pairs = collect_celltype_pairs(significant_means, query_cell_types_1, query_cell_types_2, separator)
    cols_filter = significant_means.columns[significant_means.columns.isin(cell_type_pairs)]

    # Collect all interactions from query_genes and query_interactions
    interactions = set([])
    if query_genes:
        interactions = interactions.union(frozenset(deconvoluted[deconvoluted['gene_name']
                                                    .isin(query_genes)]['id_cp_interaction'].tolist()))
    if query_interactions:
        interactions = interactions.union(frozenset(significant_means[significant_means['interacting_pair']
                                                    .isin(query_interactions)]['id_cp_interaction'].tolist()))
    if query_classifications:
        interactions = interactions.union(frozenset(significant_means[significant_means['classification']
                                                    .isin(query_classifications)]['id_cp_interaction'].tolist()))
    # If minimum_score was provided, filter interactions to those with at least minimum_score across the selected cell types
    if query_minimum_score is not None and interaction_scores is not None:
        # Filter out interactions which are below query_minimum_score in any cell_type_pair/column in cols_filter
        interactions_filtered_by_minimum_score = interaction_scores[
            interaction_scores[cols_filter].max(axis=1) >= query_minimum_score
        ]['id_cp_interaction'].tolist()
        if query_genes or query_interactions or query_classifications:
            interactions = interactions.intersection(interactions_filtered_by_minimum_score)
        else:
            # Filter all interactions by query_minimum_score
            interactions = set(interactions_filtered_by_minimum_score)
    result_df = significant_means[significant_means['id_cp_interaction'].isin(interactions)]

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
