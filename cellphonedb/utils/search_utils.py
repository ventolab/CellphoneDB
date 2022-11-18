from utils import db_utils
import re
import time
from utils.utils import dbg, write_to_file
import pandas as pd

SIMPLE_PFX="simple:"
COMPLEX_PFX = 'complex:'

def populate_proteins_for_complex(complex_name, complex_name2proteins, genes, complex_expanded, complex_composition):
    constituent_proteins = []
    if complex_name not in complex_name2proteins:
        complex_multidata_id = complex_expanded['complex_multidata_id'] \
            [complex_expanded[['name']].apply(lambda row: row.astype(str).eq(complex_name).any(), axis=1)].to_list()[0]
        for protein_multidata_id in complex_composition['protein_multidata_id'] \
                [complex_composition['complex_multidata_id'] == complex_multidata_id].to_list():
            constituent_proteins.append(genes['name'][genes['protein_multidata_id'] == protein_multidata_id].to_list()[0])
        complex_name2proteins[complex_name] = constituent_proteins

def search(query_str, user_dir_root, cellophonedb_version):
    start = time.time()
    results = []
    interactions, genes, complex_composition, complex_expanded = \
        db_utils.get_interactions_genes_complex(user_dir_root, cellophonedb_version)

    complex_name2proteins = {}
    # Assemble a list of multidata_ids to search interactions DF with
    multidata_ids = []
    for token in re.split(',\s*| ', query_str):
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
    results.append(['CellphoneDB interaction ', 'Partner A', 'Partner B', 'Gene name A', 'Gene name B', ' Ensembl ID A', 'Ensembl ID B', 'Annotation strategy'])
    for multidata_id in multidata_ids:
         interactions_data_list = interactions[[ \
             'id_cp_interaction','multidata_1_id', 'multidata_2_id', \
             'name_1', 'name_2','is_complex_1', 'is_complex_2', 'annotation_strategy']] \
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
                 interaction[7]
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

def get_html_table(data, complex_name2proteins):
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
                    html += "<td style=\"text-align:left\"><span title=\"{}\">{}</span></td>".format(complex_mouseover, name)
                elif field.startswith(SIMPLE_PFX):
                    name = field.split(":")[1]
                    html += "<td style=\"text-align:left\"><a class=\"teal-text\" href=\"https://www.uniprot.org/uniprotkb/{}/entry\">{}</a></td>" \
                        .format(name, name)
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