import io
import os
import zipfile
from distutils.version import LooseVersion
from typing import Union

import requests
from cellphonedb.src.exceptions.NoReleasesException import NoReleasesException

cpdb_releases = os.environ.get('CELLPHONEDB_RELEASE_PATH', '~/.cpdb/releases')
database_file = 'cellphone.db'


def _major(version: LooseVersion) -> int:
    for component in version.version:
        if isinstance(component, int):
            return component

def find_database_for(value: str) -> str:
    file_candidate = os.path.expanduser(value)

    if os.path.exists(file_candidate):
        print('User selected database `{}` is available, using it'.format(file_candidate))
        return file_candidate

    # _ensure_core_version_in_user_dbs()
    user_databases_prefix = os.path.expanduser(cpdb_releases)

    if not os.path.isdir(user_databases_prefix):
        print(f"WARNING: No local databases found. Will download latest database from remote to {user_databases_prefix}.")
        download_database("latest")

    if value == 'latest' or not value:
        available = list_local_versions()
        latest_available = available[0]
        print('Latest local available version is `{}`, using it'.format(latest_available))
        value = latest_available

    downloaded_candidate = os.path.join(user_databases_prefix, value, database_file)
    valid_database = os.path.exists(downloaded_candidate)

    if not valid_database:
        print("WARNING: Database '{value}' is not available. Trying to download it.")
        download_database(value)
        return find_database_for(value)

    return downloaded_candidate

def download_database(version):
    try:
        if not version or version == 'latest':
            latest_release = _latest_release()
            version = latest_release['tag']
            zip_to_download = latest_release['url']
        else:
            releases = _list_releases()
            if version not in releases:
                print(f"ERROR: Unavailable version selected. Available versions are: { ', '.join(releases.keys()) }")
                exit(1)

            zip_to_download = releases[version]['url']

        print(f"Downloading release {version} of CellPhoneDB database")

        output_folder = os.path.expanduser('{}/{}'.format(cpdb_releases, version))
        os.makedirs(output_folder, exist_ok=True)

        zip_response = requests.get(zip_to_download)
        print(f"Download completed!")
        print(f"Copying database to {output_folder}")
        with zipfile.ZipFile(io.BytesIO(zip_response.content)) as thezip:
            root_folder = thezip.namelist()[0]
            for name in thezip.namelist():
                if name.endswith('/'):
                    continue

                file_folder = os.path.dirname(name)
                file_name = os.path.basename(name)

                dest_folder = os.path.realpath(os.path.join(output_folder, os.path.relpath(file_folder, root_folder)))
                dest_file = os.path.join(dest_folder, file_name)
                os.makedirs(dest_folder, exist_ok=True)

                with thezip.open(name) as zf:
                    with open(dest_file, 'wb') as fw:
                        fw.write(zf.read())

    except NoReleasesException:
        print('There are no versions available (or connection could not be made to server to retrieve them)')
        exit(1)


def list_local_versions() -> list:
    try:
        releases_folder = os.path.expanduser(cpdb_releases)
        compatible_versions = [ rf for rf in os.listdir(releases_folder) if os.path.isdir(os.path.join(releases_folder,rf)) ]
        return sorted(compatible_versions, key=LooseVersion, reverse=True)
    except FileNotFoundError:
        return []


def list_remote_database_versions():
    try:
        releases = _list_releases()

        for idx, (_, version) in enumerate(releases.items()):
            note = ' *latest' if idx == 0 else ''
            print('version {}{}: released: {}, url: {}, compatible: {}'.format(version['tag'], note, version['date'],
                                                                               version['link'], version['compatible']))

    except NoReleasesException:
        print('There are no versions available (or connection could not be made to server to retrieve them)')
        exit(1)

def _list_releases() -> dict:
    try:
        releases = _github_query('releases')

        if not releases:
            raise NoReleasesException()

        return _format_releases(*releases)

    except requests.exceptions.ConnectionError:
        raise NoReleasesException()


def _latest_release():
    try:
        compatible_versions = {key: value for key, value in _list_releases().items() if value['compatible']}

        if not compatible_versions:
            raise NoReleasesException()

        latest = sorted(compatible_versions, key=LooseVersion, reverse=True)[0]

        return compatible_versions[latest]

    except requests.exceptions.ConnectionError:
        raise NoReleasesException()


def _github_query(kind) -> Union[dict, list]:
    queries = {
        'releases': 'https://api.github.com/repos/{}/{}/releases'.format('ventolab', 'cellphonedb-data'),
    }

    query = queries[kind]

    if not query:
        raise Exception('unexpected query to github')

    response = requests.get(query, headers={'Accept': 'application/vnd.github.v3+json'})

    if response.status_code != 200:
        raise NoReleasesException()

    return response.json()


def _format_releases(*releases) -> dict:
    return {item['tag_name']: _format_release(item) for item in releases}


def _format_release(item: dict) -> dict:
    tag_name = item['tag_name']

    return {'url': item['zipball_url'],
            'tag': tag_name,
            'date': item['published_at'],
            'link': item['html_url'],
            'compatible': True
            }
