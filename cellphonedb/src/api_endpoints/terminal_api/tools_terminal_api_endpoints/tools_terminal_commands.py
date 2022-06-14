import os
import urllib.parse
from typing import Optional

import click
import numpy as np
import pandas as pd

from cellphonedb.src.app.app_logger import app_logger
from cellphonedb.src.app.cellphonedb_app import output_dir, data_dir
from cellphonedb.src.core.generators.complex_generator import complex_generator
from cellphonedb.src.core.generators.gene_generator import gene_generator
from cellphonedb.src.core.generators.protein_generator import protein_generator
from cellphonedb.tools import tools_helper
from cellphonedb.tools.generate_data.filters.non_complex_interactions import only_noncomplex_interactions
from cellphonedb.tools.generate_data.filters.remove_interactions import remove_interactions_in_file
from cellphonedb.tools.generate_data.getters import get_iuphar, get_imex
from cellphonedb.tools.generate_data.mergers.add_curated import add_curated
from cellphonedb.tools.generate_data.mergers.merge_interactions import merge_iuphar_imex_interactions
from cellphonedb.tools.generate_data.parsers import parse_iuphar_guidetopharmacology
from cellphonedb.tools.generate_data.parsers.parse_interactions_imex import parse_interactions_imex
from cellphonedb.utils import utils
from cellphonedb.utils.utils import _get_separator, write_to_file


@click.command("generate_genes")
@click.option('--user-gene', type=click.Path(file_okay=True, exists=True, dir_okay=False))
@click.option('--fetch-uniprot', is_flag=True)
@click.option('--fetch-ensembl', is_flag=True)
@click.option('--result-path', type=str, default=None)
@click.option('--project-name', type=str, default=None)
def generate_genes(user_gene: Optional[str],
                   fetch_uniprot: bool,
                   fetch_ensembl: bool,
                   result_path: str,
                   project_name: str,
                   ) -> None:
    output_path = _set_paths(result_path, project_name)

    # TODO: Add logger
    if fetch_ensembl:
        print('fetching remote ensembl data ... ', end='')
        source_url = 'http://www.ensembl.org/biomart/martservice?query={}'
        query = '<?xml version="1.0" encoding="UTF-8"?><!DOCTYPE Query><Query virtualSchemaName = "default" ' \
                'formatter = "CSV" header = "1" uniqueRows = "1" count = "" datasetConfigVersion = "0.6" >' \
                '<Dataset name = "hsapiens_gene_ensembl" interface = "default" >' \
                '<Attribute name = "ensembl_gene_id" />' \
                '<Attribute name = "ensembl_transcript_id" />' \
                '<Attribute name = "external_gene_name" />' \
                '<Attribute name = "hgnc_symbol" />' \
                '<Attribute name = "uniprotswissprot" />' \
                '</Dataset>' \
                '</Query>'

        url = source_url.format(urllib.parse.quote(query))
        ensembl_db = pd.read_csv(url)
        print('done')
    else:
        ensembl_db = utils.read_data_table_from_file(os.path.join(data_dir, 'sources/ensembl.txt'))
        print('read local ensembl file')

    # additional data comes from given file or uniprot remote url
    if fetch_uniprot:
        try:
            print('fetching remote uniprot file ... ', end='')
            source_url = 'https://www.uniprot.org/uniprot/?query=*&format=tab&force=true' \
                         '&columns=id,entry%20name,reviewed,protein%20names,genes,organism,length' \
                         '&fil=organism:%22Homo%20sapiens%20(Human)%20[9606]%22%20AND%20reviewed:yes' \
                         '&compress=yes'

            uniprot_db = pd.read_csv(source_url, sep='\t', compression='gzip')
            print('done')

        except Exception as e:
            print('Error fetching remote UniProt data, fetching local data')
            uniprot_db = pd.read_csv(os.path.join(data_dir, 'sources/uniprot.tab'), sep='\t')
            print('read local uniprot file')
    else:
        uniprot_db = utils.read_data_table_from_file(os.path.join(data_dir, 'sources/uniprot.tab'))
        print('read local uniprot file')

    ensembl_columns = {
        'Gene name': 'gene_name',
        'Gene stable ID': 'ensembl',
        'HGNC symbol': 'hgnc_symbol',
        'UniProtKB/Swiss-Prot ID': 'uniprot'
    }

    uniprot_columns = {
        'Entry': 'uniprot',
        'Gene names': 'gene_names'
    }

    result_columns = [
        'gene_name',
        'uniprot',
        'hgnc_symbol',
        'ensembl'
    ]

    ensembl_db = ensembl_db[list(ensembl_columns.keys())].rename(columns=ensembl_columns)
    uniprot_db = uniprot_db[list(uniprot_columns.keys())].rename(columns=uniprot_columns)
    hla_genes = utils.read_data_table_from_file(os.path.join(data_dir, 'sources/hla_curated.csv'))
    if user_gene:
        separator = _get_separator(os.path.splitext(user_gene)[-1])
        user_gene = pd.read_csv(user_gene, sep=separator)

    cpdb_genes = gene_generator(ensembl_db, uniprot_db, hla_genes, user_gene, result_columns)

    cpdb_genes[result_columns].to_csv('{}/{}'.format(output_path, 'gene_generated.csv'), index=False)


@click.command("generate_interactions")
@click.argument('proteins', default='protein.csv', type=click.Path(file_okay=True, exists=True, dir_okay=False))
@click.argument('genes', default='gene.csv', type=click.Path(file_okay=True, exists=True, dir_okay=False))
@click.argument('complex', default='complex.csv', type=click.Path(file_okay=True, exists=True, dir_okay=False))
@click.option('--user-interactions', type=click.Path(file_okay=True, exists=True, dir_okay=False))
@click.option('--user-interactions-only', is_flag=True)
@click.option('--result-path', type=str, default=None)
@click.option('--fetch-imex', is_flag=True)
@click.option('--fetch-iuphar', is_flag=True)
@click.option('--project-name', type=str, default=None)
@click.option('--release', is_flag=True)
def generate_interactions(proteins: str,
                          genes: str,
                          complex: str,
                          user_interactions: Optional[str],
                          user_interactions_only: bool,
                          result_path: str,
                          fetch_imex: bool,
                          fetch_iuphar: bool,
                          project_name: str,
                          release: bool
                          ) -> None:
    if user_interactions_only and not user_interactions:
        print("'--user-interactions-only' option detected but no user intearctions file was provided with '--user-interactions'")
        raise Exception('You need to set --user-interactions parameter')

    output_path = utils.set_paths(result_path, project_name)
    downloads_path = utils.set_paths(utils.set_paths(result_path, project_name), 'downloads')

    proteins = utils.read_data_table_from_file(proteins)
    genes = utils.read_data_table_from_file(genes)
    complexes = utils.read_data_table_from_file(complex)

    if not user_interactions_only:
        raw_imex = get_imex.call(genes, downloads_path, fetch_imex)

        interactions_to_remove = utils.read_data_table_from_file(
            os.path.join(data_dir, 'sources/excluded_interaction.csv'))
        interaction_curated = utils.read_data_table_from_file(os.path.join(data_dir, 'sources/interaction_curated.csv'))

    if user_interactions:
        separator = _get_separator(os.path.splitext(user_interactions)[-1])
        user_interactions = pd.read_csv(user_interactions, sep=separator)
        user_interactions['partner_a'] = user_interactions['partner_a'].apply(lambda x: str(x).strip())
        user_interactions['partner_b'] = user_interactions['partner_b'].apply(lambda x: str(x).strip())

        # if --release mark all user interactions as 'curated otherwise keep
        # the ones from CellPhoneDB as 'curated' and the new ones by the user 
        # as 'user_curated'
        if release:
            user_interactions['annotation_strategy'] = 'curated'
        else:
            user_interactions['annotation_strategy'] = 'user_curated'

        # add missing proteing name columns
        if not 'protein_name_a' in user_interactions.columns:
            user_interactions['protein_name_a'] = ''

        if not 'protein_name_b' in user_interactions.columns:
            user_interactions['protein_name_b'] = ''

    result_columns = [
        'partner_a',
        'partner_b',
        'protein_name_a',
        'protein_name_b',
        'annotation_strategy',
        'source'
    ]
    if not user_interactions_only:
        print('Parsing IMEX file')
        imex_interactions = parse_interactions_imex(raw_imex, proteins, genes)

        print('Getting iuphar data')
        raw_iuphar = get_iuphar.call(downloads_path, fetch_iuphar)

        print('Generating iuphar interactions')
        iuphar_interactions = parse_iuphar_guidetopharmacology.call(raw_iuphar, genes, proteins)

        print('Merging iuphar/imex')
        merged_interactions = merge_iuphar_imex_interactions(iuphar_interactions, imex_interactions)

        print('Removing complex interactions')
        no_complex_interactions = only_noncomplex_interactions(merged_interactions, complexes)

        print('Removing selected interactions')
        clean_interactions = remove_interactions_in_file(no_complex_interactions, interactions_to_remove)

        if not release:
            # if not --release we'll add the curated interactions that come with CellPhoneDB
            print('Adding curated interaction')
            interactions_with_curated = add_curated(clean_interactions, interaction_curated)

            result = tools_helper.normalize_interactions(
                interactions_with_curated.append(user_interactions, ignore_index=True, sort=False), 'partner_a',
                'partner_b').drop_duplicates(['partner_a', 'partner_b'], keep='last')
        else:
            # if --release we'll drop previously curated interactions and will keep only 
            # the ones provided by the user and the ones from third party sources
            user_interactions = add_curated(clean_interactions, user_interactions)
            result = tools_helper.normalize_interactions(user_interactions,'partner_a',
                'partner_b').drop_duplicates(['partner_a', 'partner_b'], keep='last')            

    else:
        result = tools_helper.normalize_interactions(user_interactions, 'partner_a', 'partner_b') \
            .drop_duplicates(['partner_a', 'partner_b'], keep='last')

    result[result_columns].sort_values(['partner_a', 'partner_b']).to_csv(
        '{}/interaction_input.csv'.format(output_path), index=False)


@click.command("generate_proteins")
@click.option('--user-protein', type=click.Path(file_okay=True, exists=True, dir_okay=False))
@click.option('--fetch-uniprot', is_flag=True)
@click.option('--result-path', type=str, default=None)
@click.option('--log-file', type=str, default='log.txt')
@click.option('--project-name', type=str, default=None)
def generate_proteins(user_protein: Optional[str],
                      fetch_uniprot: bool,
                      result_path: str,
                      log_file: str,
                      project_name: str):
    uniprot_columns = {
        'Entry': 'uniprot',
        'Entry name': 'protein_name',
    }

    # additional data comes from given file or uniprot remote url
    if fetch_uniprot:
        try:
            source_url = 'https://www.uniprot.org/uniprot/?query=*&format=tab&force=true' \
                         '&columns=id,entry%20name,reviewed,protein%20names,genes,organism,length' \
                         '&fil=organism:%22Homo%20sapiens%20(Human)%20[9606]%22%20AND%20reviewed:yes' \
                         '&compress=yes'

            uniprot_db = pd.read_csv(source_url, sep='\t', compression='gzip')

            print('read remote uniprot file')
        except Exception as e:
            print('Error fetching remote UniProt data, fetching local data')
            uniprot_db = pd.read_csv(os.path.join(data_dir, 'sources/uniprot.tab'), sep='\t')
            print('read local uniprot file')
    else:
        uniprot_db = pd.read_csv(os.path.join(data_dir, 'sources/uniprot.tab'), sep='\t')
        print('read local uniprot file')

    default_values = {
        'transmembrane': False,
        'peripheral': False,
        'secreted': False,
        'secreted_desc': np.nan,
        'secreted_highlight': False,
        'receptor': False,
        'receptor_desc': np.nan,
        'integrin': False,
        'other': False,
        'other_desc': np.nan,
        'tags': 'To_add',
        'tags_reason': np.nan,
        'tags_description': np.nan,
        'pfam': np.nan,
    }

    default_types = {
        'uniprot': str,
        'protein_name': str,
        'transmembrane': bool,
        'peripheral': bool,
        'secreted': bool,
        'secreted_desc': str,
        'secreted_highlight': bool,
        'receptor': bool,
        'receptor_desc': str,
        'integrin': bool,
        'other': bool,
        'other_desc': str,
        'tags': str,
        'tags_reason': str,
        'tags_description': str,
        'pfam': str,
    }

    result_columns = list(default_types.keys())

    output_path = _set_paths(result_path, project_name)
    log_path = '{}/{}'.format(output_path, log_file)
    uniprot_db = uniprot_db[list(uniprot_columns.keys())].rename(columns=uniprot_columns)
    curated_proteins = pd.read_csv(os.path.join(data_dir, 'sources/protein_curated.csv'))
    if user_protein:
        separator = _get_separator(os.path.splitext(user_protein)[-1])
        user_protein = pd.read_csv(user_protein, sep=separator)

    result = protein_generator(uniprot_db, curated_proteins, user_protein, default_values, default_types,
                               result_columns, log_path)

    result[result_columns].to_csv('{}/{}'.format(output_path, 'protein_generated.csv'), index=False)


@click.command("generate_complex")
@click.option('--user-complex', type=click.Path(file_okay=True, exists=True, dir_okay=False))
@click.option('--result-path', type=str, default=None)
@click.option('--log-file', type=str, default='log.txt')
@click.option('--project-name', type=str, default=None)
def generate_complex(user_complex: Optional[str],
                     result_path: str,
                     log_file: str,
                     project_name: str):
    output_path = _set_paths(result_path, project_name)
    log_path = '{}/{}'.format(output_path, log_file)

    curated_complex = pd.read_csv(os.path.join(data_dir, 'sources/complex_curated.csv'))
    if user_complex:
        separator = _get_separator(os.path.splitext(user_complex)[-1])
        user_complex = pd.read_csv(user_complex, sep=separator)

    result = complex_generator(curated_complex, user_complex, log_path)

    result.to_csv('{}/{}'.format(output_path, 'complex_generated.csv'), index=False)


@click.command("filter_all")
@click.option('--input-path', type=str, default=data_dir)
@click.option('--result-path', type=str, default='filtered')
@click.option('--project-name', type=str, default=None)
def filter_all(input_path,
               result_path,
               project_name: str):
    interactions = pd.read_csv(os.path.join(input_path, 'interaction_input.csv'))
    complexes = pd.read_csv(os.path.join(input_path, 'complex_generated.csv'))
    proteins = pd.read_csv(os.path.join(input_path, 'protein_generated.csv'))
    genes = pd.read_csv(os.path.join(input_path, 'gene_generated.csv'))
    output_path = _set_paths(result_path, project_name)

    interacting_partners = pd.concat([interactions['partner_a'], interactions['partner_b']]).drop_duplicates()

    filtered_complexes = _filter_complexes(complexes, interacting_partners)
    write_to_file(filtered_complexes.sort_values('complex_name'), 'complex_input.csv', output_path=output_path)

    filtered_proteins, interacting_proteins = _filter_proteins(proteins, filtered_complexes, interacting_partners)
    write_to_file(filtered_proteins.sort_values('uniprot'), 'protein_input.csv', output_path=output_path)

    filtered_genes = _filter_genes(genes, filtered_proteins['uniprot'])
    write_to_file(filtered_genes.sort_values('gene_name'), 'gene_input.csv', output_path=output_path)

    rejected_members = interacting_partners[~(interacting_partners.isin(filtered_complexes['complex_name']) |
                                              interacting_partners.isin(filtered_proteins['uniprot']))]

    if len(rejected_members):
        app_logger.warning('There are some proteins or complexes not interacting properly: `{}`'.format(
            ', '.join(rejected_members)))


def _filter_genes(genes: pd.DataFrame, interacting_proteins: pd.Series) -> pd.DataFrame:
    filtered_genes = genes[genes['uniprot'].isin(interacting_proteins)]

    return filtered_genes


def _filter_proteins(proteins: pd.DataFrame,
                     filtered_complexes: pd.DataFrame,
                     interacting_partners: pd.DataFrame) -> (pd.DataFrame, pd.DataFrame):
    interacting_proteins = pd.concat(
        [filtered_complexes['uniprot_{}'.format(i)] for i in range(1, 5)]).drop_duplicates()

    filtered_proteins = proteins[
        proteins['uniprot'].isin(interacting_partners) | proteins['uniprot'].isin(interacting_proteins)]

    return filtered_proteins, interacting_proteins


def _filter_complexes(complexes: pd.DataFrame, interacting_partners: pd.DataFrame) -> pd.DataFrame:
    filtered_complexes = complexes[complexes['complex_name'].isin(interacting_partners)]

    return filtered_complexes


def _set_paths(output_path, subfolder):
    if not output_path:
        output_path = output_dir

    if subfolder:
        output_path = os.path.realpath(os.path.expanduser('{}/{}'.format(output_path, subfolder)))

    os.makedirs(output_path, exist_ok=True)

    if _path_is_not_empty(output_path):
        app_logger.warning(
            'Output directory ({}) exist and is not empty. Result can overwrite old results'.format(output_path))

    return output_path


def _path_is_not_empty(path):
    return bool([f for f in os.listdir(path) if not f.startswith('.')])
