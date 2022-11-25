import os
import urllib.parse

import numpy as np
import pandas as pd
import pathlib
import json

from cellphonedb.src.core.generators.complex_generator import complex_generator
from cellphonedb.src.core.generators.gene_generator import gene_generator
from cellphonedb.src.core.generators.protein_generator import protein_generator
from cellphonedb.utils import utils, generate_input_files_helper, db_utils
from cellphonedb.utils.utils import _get_separator, write_to_file

def generate_genes(data_dir,
                   user_gene=None,
                   fetch_uniprot=True,
                   fetch_ensembl=True,
                   result_path=None,
                   project_name=None,
                   ) -> None:
    output_path = utils.set_paths(result_path, project_name)
    if fetch_ensembl:
        print('Fetching remote Ensembl data ... ', end='')
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
        print('Done')
    else:
        ensembl_db = utils.read_data_table_from_file(os.path.join(data_dir, 'sources','ensembl.txt'))
        print('Read local Ensembl file')

    # additional data comes from given file or uniprot remote url
    if fetch_uniprot:
        print('Fetching remote UniProt file ... ', end='')
        try:
            source_url = 'https://www.uniprot.org/uniprot/?query=*&format=tab&force=true' \
                         '&columns=id,entry%20name,reviewed,protein%20names,genes,organism,length' \
                         '&fil=organism:%22Homo%20sapiens%20(Human)%20[9606]%22%20AND%20reviewed:yes' \
                         '&compress=yes'

            uniprot_db = pd.read_csv(source_url, sep='\t', compression='gzip')
            print('done')
        except Exception as e:
            print('Error fetching remote UniProt data, fetching local data')
            uniprot_db = pd.read_csv(os.path.join(data_dir, 'sources/uniprot.tab'), sep='\t')
            print('Read local UniProt file')
    else:
        uniprot_db = utils.read_data_table_from_file(os.path.join(data_dir, 'sources','uniprot.tab'))
        print('Read local UniProt file')

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
    hla_genes = utils.read_data_table_from_file(os.path.join(data_dir, 'sources','hla_curated.csv'))
    if user_gene:
        separator = _get_separator(os.path.splitext(user_gene)[-1])
        user_gene = pd.read_csv(user_gene, sep=separator)

    cpdb_genes = gene_generator(ensembl_db, uniprot_db, hla_genes, user_gene, result_columns)

    cpdb_genes[result_columns].to_csv('{}/{}'.format(output_path, 'gene_generated.csv'), index=False)

def generate_interactions(data_dir,
                          user_interactions=None,
                          user_interactions_only=False,
                          result_path=None,
                          project_name=None,
                          release=True
                          ) -> None:
    if user_interactions_only and not user_interactions:
        print("user-interactions-only' option detected but no user interactions file was provided with '--user-interactions'")
        raise Exception('You need to set user-interactions parameter')

    output_path = utils.set_paths(result_path, project_name)

    if not user_interactions_only:
        interaction_curated = utils.read_data_table_from_file(os.path.join(data_dir, 'sources','interaction_curated.csv'))

    if user_interactions:
        separator = _get_separator(os.path.splitext(user_interactions)[-1])
        user_interactions = pd.read_csv(user_interactions, sep=separator)
        user_interactions['partner_a'] = user_interactions['partner_a'].apply(lambda x: str(x).strip())
        user_interactions['partner_b'] = user_interactions['partner_b'].apply(lambda x: str(x).strip())

        # if release == True mark all user interactions as 'curated'; otherwise keep
        # the ones from CellPhoneDB as 'curated' and the new ones by the user 
        # as 'user_curated'
        if release:
            user_interactions['annotation_strategy'] = 'curated'
        else:
            user_interactions['annotation_strategy'] = 'user_curated'

        # add missing protein_name columns
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
        if not release:
            # if release == False, append user_interactions to interaction_curated
            result = generate_input_files_helper.normalize_interactions(
                pd.concat([interaction_curated, user_interactions], ignore_index=True, sort=False), 'partner_a',
                'partner_b').drop_duplicates(['partner_a', 'partner_b'], keep='last')
        else:
            # if release == True keep user_interactions only
            result = generate_input_files_helper.normalize_interactions(user_interactions,'partner_a',
                'partner_b').drop_duplicates(['partner_a', 'partner_b'], keep='last')            

    else:
        # user_interactions_only == True -> keep user_interactions only
        result = generate_input_files_helper.normalize_interactions(user_interactions, 'partner_a', 'partner_b') \
            .drop_duplicates(['partner_a', 'partner_b'], keep='last')

    result[result_columns].sort_values(['partner_a', 'partner_b']).to_csv(
        '{}/interaction_input.csv'.format(output_path), index=False)

def generate_proteins(data_dir,
                      user_protein=None,
                      fetch_uniprot=True,
                      result_path=None,
                      log_file="log.txt",
                      project_name=None):
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

            print('Downloaded remote UniProt file')
        except Exception as e:
            print('Error fetching remote UniProt data, fetching local data')
            uniprot_db = pd.read_csv(os.path.join(data_dir, 'sources','uniprot.tab'), sep='\t')
            print('Read local UniProt file')
    else:
        uniprot_db = pd.read_csv(os.path.join(data_dir, 'sources','uniprot.tab'), sep='\t')
        print('Read local UniProt file')

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

    output_path = utils.set_paths(result_path, project_name)
    log_path = '{}/{}'.format(output_path, log_file)
    uniprot_db = uniprot_db[list(uniprot_columns.keys())].rename(columns=uniprot_columns)
    curated_proteins = pd.read_csv(os.path.join(data_dir, 'sources','protein_curated.csv'))
    if user_protein:
        separator = _get_separator(os.path.splitext(user_protein)[-1])
        user_protein = pd.read_csv(user_protein, sep=separator)

    result = protein_generator(uniprot_db, curated_proteins, user_protein, default_values, default_types,
                               result_columns, log_path)

    result[result_columns].to_csv('{}/{}'.format(output_path, 'protein_generated.csv'), index=False)

def generate_complex(data_dir,
                     user_complex=None,
                     result_path=None,
                     log_file='log.txt',
                     project_name=None):
    output_path = utils.set_paths(result_path, project_name)
    log_path = '{}/{}'.format(output_path, log_file)

    curated_complex = pd.read_csv(os.path.join(data_dir, 'sources','complex_curated.csv'))
    if user_complex:
        separator = _get_separator(os.path.splitext(user_complex)[-1])
        user_complex = pd.read_csv(user_complex, sep=separator)

    result = complex_generator(curated_complex, user_complex, log_path)

    result.to_csv('{}/{}'.format(output_path, 'complex_generated.csv'), index=False)

def filter_all(data_dir,
               input_path,
               user_complex=None,
               user_interaction=None,
               result_path='filtered',
               project_name=None):
    interactions = pd.read_csv(os.path.join(input_path, 'interaction_input.csv'))
    complexes = pd.read_csv(os.path.join(input_path, 'complex_generated.csv'))
    proteins = pd.read_csv(os.path.join(input_path, 'protein_generated.csv'))
    genes = pd.read_csv(os.path.join(input_path, 'gene_generated.csv'))
    output_path = utils.set_paths(result_path, project_name)

    interacting_partners = pd.concat([interactions['partner_a'], interactions['partner_b']]).drop_duplicates()

    filtered_complexes = _filter_complexes(complexes, interacting_partners)
    write_to_file(filtered_complexes.sort_values('complex_name'), 'complex_input.csv', output_path=output_path)

    filtered_proteins, interacting_proteins = _filter_proteins(proteins, filtered_complexes, interacting_partners)
    write_to_file(filtered_proteins.sort_values('uniprot'), 'protein_input.csv', output_path=output_path)

    filtered_genes = _filter_genes(genes, filtered_proteins['uniprot'])
    write_to_file(filtered_genes.sort_values('gene_name'), 'gene_input.csv', output_path=output_path)

    unknown_interactors = interacting_partners[~(interacting_partners.isin(filtered_complexes['complex_name']) |
                                              interacting_partners.isin(filtered_proteins['uniprot']))]


    if len(unknown_interactors):
        print('WARNING: Some interactions in {} were excluded due to unknown complexes or proteins. Please include any missing complexes in {}, and use proteins present in {} only:\n`{}`'
              .format(user_interaction, user_complex,
                      os.path.join(data_dir, 'sources', 'uniprot.tab'), ', '.join(unknown_interactors)))

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

def download_source_files(user_dir_root, db_version):
    data_dir = os.path.join(db_utils.get_db_path(user_dir_root, db_version), "data")
    sources_path = os.path.join(data_dir, "sources")
    print("Downloading cellphonedb-data/data/sources files into {}:".format(sources_path))
    pathlib.Path(sources_path).mkdir(parents=True, exist_ok=True)
    r = urllib.request.urlopen("https://api.github.com/repos/ventolab/cellphonedb-data/git/trees/master?recursive=1")
    files_data = json.load(r)['tree']
    for rec in files_data:
        if rec['path'].startswith('data/sources/'):
            fname = rec['path'].split('/')[-1]
            url = 'https://raw.githubusercontent.com/ventolab/cellphonedb-data/master/{}'.format(rec['path'])
            print("Downloading: " + fname)
            r = urllib.request.urlopen(url)
            with open(os.path.join(sources_path, fname), 'wb') as f:
                f.write(r.read())

def generate_all(user_dir_root, db_version, user_complex=None, user_interactions=None, user_interactions_only=False ):
    download_source_files(user_dir_root, db_version)
    db_dir = db_utils.get_db_path(user_dir_root, db_version)
    data_dir = os.path.join(db_dir,"data")
    generated_path = os.path.join(data_dir, "generated")
    print("Generating gene_generated.csv file into {}".format(generated_path))
    pathlib.Path(generated_path).mkdir(parents=True, exist_ok=True)
    generate_genes(data_dir,
                   user_gene=None,
                   fetch_uniprot=False,
                   fetch_ensembl=False,
                   result_path=generated_path,
                   project_name=None,
                   )
    print("Generating proteins_generated.csv file into {}".format(generated_path))
    generate_proteins(data_dir,
                      user_protein=None,
                      fetch_uniprot=False,
                      result_path=generated_path,
                      log_file="log.txt",
                      project_name=None)
    print("Generating complex_generated.csv file into {}".format(generated_path))
    generate_complex(data_dir,
                     user_complex=user_complex,
                     result_path=generated_path,
                     log_file='log.txt',
                     project_name=None)
    print("Generating interactions_input.csv file into {}".format(generated_path))
    generate_interactions(data_dir,
                          user_interactions=user_interactions,
                          user_interactions_only=user_interactions_only,
                          result_path=generated_path,
                          project_name=None,
                          release=False)
    print("Generating gene, protein and complex input files file into {}".format(generated_path))
    filter_all(data_dir,
               input_path=generated_path,
               user_complex=user_complex,
               user_interaction=user_interactions,
               result_path=generated_path,
               project_name=None)
