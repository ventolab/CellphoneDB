import os
from io import StringIO

import pandas as pd
import requests
import tqdm

from utils.generate_data_helper import add_to_meta

def call(data_dir, genes: pd.DataFrame, downloads_path: str, fetch: bool, save_backup: bool = True) -> pd.DataFrame:
    proteins = genes['uniprot'].tolist()

    sources = [
        {
            'name': 'InnateDB',
            'base_url': 'https://psicquic.curated.innatedb.com/webservices/current/search/query/species:human',
            'query_parameters': False,
        },
        {
            'name': 'InnateDB-All',
            'base_url': 'https://psicquic.all.innatedb.com/webservices/current/search/query/species:human',
            'query_parameters': False,
        },
        {
            'name': 'IMEx',
            'base_url': 'http://www.ebi.ac.uk/Tools/webservices/psicquic/imex/webservices/current/search/query/{}',
            'query_parameters': True,
        },
        {
            'name': 'IntAct',
            'base_url': 'http://www.ebi.ac.uk/Tools/webservices/psicquic/intact/webservices/current/search/query/{}',
            'query_parameters': True,
        },
        {
            'name': 'bhf-ucl',
            'base_url': 'http://www.ebi.ac.uk/Tools/webservices/psicquic/bhf-ucl/webservices/current/search/query/{}',
            'query_parameters': True,
        },
        {
            'name': 'MatrixDB',
            'base_url': 'http://matrixdb.univ-lyon1.fr:8080/psicquic/webservices/current/search/query/{}',
            'query_parameters': True,
        },
        {
            'name': 'MINT',
            'base_url': 'http://www.ebi.ac.uk/Tools/webservices/psicquic/mint/webservices/current/search/query/{}',
            'query_parameters': True,
        },
        {
            'name': 'I2D',
            'base_url': 'http://ophid.utoronto.ca/psicquic-ws/webservices/current/search/query/{}',
            'query_parameters': True,
        },
        {
            'name': 'UniProt',
            'base_url': 'http://www.ebi.ac.uk/Tools/webservices/psicquic/uniprot/webservices/current/search/query/{}',
            'query_parameters': True,
        },
        {
            'name': 'MBInfo',
            'base_url': 'http://www.ebi.ac.uk/Tools/webservices/psicquic/mbinfo/webservices/current/search/query/{}',
            'query_parameters': True,
        },
        # 'DIP': 'http://imex.mbi.ucla.edu/xpsq-dip-all/service/rest/current/search/query/{}',
    ]

    significant_columns = ['A', 'B', 'altA', 'altB', 'provider']
    carry = pd.DataFrame(columns=significant_columns)

    for source in tqdm.tqdm(sources, desc='Fetching sources'.ljust(20)):
        carry = carry.append(_get_source(data_dir, source, proteins, downloads_path, significant_columns, fetch, save_backup),
                             ignore_index=True,
                             sort=False)
    print()

    carry.drop_duplicates(inplace=True)

    return carry


def _get_source(data_dir, source, proteins, downloads_path, significant_columns, fetch: bool, save_backup: bool = True):
    carry = pd.DataFrame(columns=significant_columns)
    columns_to_save = ['A', 'B', 'altA', 'altB']
    compression = 'xz'
    file_name = '{}_interaction_raw.csv.{}'.format(source['name'], compression)
    download_file_path = os.path.join(downloads_path, file_name)

    def best_path_for(name: str):
        saved_file_path = os.path.join(data_dir, 'sources', name)

        if os.path.exists(download_file_path):
            return download_file_path
        if os.path.exists(saved_file_path):
            return saved_file_path

        print('ERROR: Could not find local source for {}'.format(name))
        exit(1)

    file_path = best_path_for(file_name)

    try:
        if fetch:
            if source['query_parameters']:
                carry = _get_chunked_api_results(carry, columns_to_save, proteins, source)
            else:
                carry = _get_single_api_results(carry, columns_to_save, source)

            carry.drop_duplicates(inplace=True)
            if save_backup:
                carry.to_csv(download_file_path, index=False, compression=compression)
                add_to_meta(file_name, os.path.join(downloads_path, 'meta.json'))

        else:
            tqdm.tqdm.write('Using local version for source {}'.format(source['name']))
            carry = pd.read_csv(file_path, compression=compression)

    except CouldNotFetchFromApiException:
        tqdm.tqdm.write('Could not fetch remote source {}, using available backup'.format(source['name']))
        carry = pd.read_csv(file_path, compression=compression)

    carry['provider'] = source['name']

    return carry


class CouldNotFetchFromApiException(Exception):
    pass


def _get_chunked_api_results(carry, columns_to_save, proteins, source):
    chunk_size = 500
    for idx, chunk in enumerate(tqdm.tqdm(zip(*[iter(proteins)] * chunk_size),
                                          desc=source['name'].ljust(20),
                                          total=int(len(proteins) / chunk_size))):
        url = source['base_url'].format(' or '.join(chunk))
        try:
            response = requests.get(url)

            if response.text:
                s = StringIO(response.text)
                df = pd.read_csv(s, sep='\t', names=columns_to_save, usecols=range(4), na_values='-')
                carry = pd.concat([carry, df], sort=False, axis=0, ignore_index=True)
                carry.drop_duplicates(inplace=True)
            else:
                if response.status_code != 200:
                    raise CouldNotFetchFromApiException()
        except (requests.exceptions.ConnectionError, requests.exceptions.SSLError):
            raise CouldNotFetchFromApiException()

    return carry


def _get_single_api_results(carry, columns_to_save, source):
    url = source['base_url']
    tqdm.tqdm.write('Fetching {}'.format(source['name']))
    try:
        response = requests.get(url)
        if response.text:
            s = StringIO(response.text)
            df = pd.read_csv(s, sep='\t', names=columns_to_save, usecols=range(4), na_values='-')
            carry = pd.concat([carry, df], sort=False, axis=0, ignore_index=True)
            carry.drop_duplicates(inplace=True)
        else:
            if response.status_code != 200:
                raise CouldNotFetchFromApiException()

    except (requests.exceptions.ConnectionError, requests.exceptions.SSLError):
        raise CouldNotFetchFromApiException()

    return carry
