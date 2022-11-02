import os
from io import StringIO

import pandas as pd
import requests

from utils.generate_data_helper import add_to_meta
from utils.generate_data.getters.get_imex import CouldNotFetchFromApiException


def call(data_dir, downloads_path: str, fetch: bool, save_backup: bool = True) -> pd.DataFrame:
    url = 'http://www.guidetopharmacology.org/DATA/interactions.csv'

    compression = 'xz'
    file_name = 'iuphar_interaction_raw.csv.{}'.format(compression)
    download_file_path = os.path.join(downloads_path, file_name)

    def best_path():
        saved_file_path = os.path.join(data_dir, 'sources', file_name)

        if os.path.exists(download_file_path):
            return download_file_path
        if os.path.exists(saved_file_path):
            return saved_file_path

        print('ERROR: Could not find local source for iuphar')
        exit(1)

    try:
        if fetch:
            response = requests.get(url)

            if response.text:
                s = StringIO(response.text)
                df = pd.read_csv(s, dtype=str)

                df.drop_duplicates(inplace=True)
                if save_backup:
                    df.to_csv(download_file_path, index=False, compression=compression)
                    add_to_meta(file_name, os.path.join(downloads_path, 'meta.json'))

                return df
            else:
                if response.status_code != 200:
                    raise CouldNotFetchFromApiException()
        else:
            print('WARNING: Using local version for iuphar')
            df = pd.read_csv(best_path(), compression=compression, dtype=str)
            return df
    except (requests.exceptions.ConnectionError, requests.exceptions.SSLError):
        core_logger.warning('Could not fetch remote source for iuphar, using local backup')
        df = pd.read_csv(best_path(), compression=compression)
        return df
