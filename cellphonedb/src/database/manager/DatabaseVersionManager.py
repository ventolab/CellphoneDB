import io
#import json
import os
import zipfile
from distutils.dir_util import copy_tree
from distutils.file_util import copy_file
from distutils.version import LooseVersion
from typing import Union

import requests
from cellphonedb.src.app.app_logger import app_logger
from cellphonedb.src.app.cpdb_app import create_app
from cellphonedb.src.exceptions.NoReleasesException import NoReleasesException
from cellphonedb.src.local_launchers.local_collector_launcher import LocalCollectorLauncher

cpdb_releases = '~/.cpdb/releases'
database_file = 'cellphone.db'


def _major(version: LooseVersion) -> int:
    for component in version.version:
        if isinstance(component, int):
            return component


# def _get_core_version() -> LooseVersion:
#     #with open(os.path.join(core_dir, 'metadata.json')) as metadata_file:
#     #    metadata = json.load(metadata_file)
#     #    core_version = metadata.get('database_version', 'core')
    
#         return LooseVersion(core_version)


# def _ensure_core_version_in_user_dbs():
#     core_version = _get_core_version()

#     user_databases_prefix = os.path.expanduser(cpdb_releases)
#     dest_folder = os.path.join(user_databases_prefix, str(core_version))
#     database_file_location = os.path.join(dest_folder, database_file)

#     if os.path.isfile(database_file_location):
#         return

#     os.makedirs(dest_folder, exist_ok=True)

#     copy_file(os.path.join(core_dir, database_file), dest_folder)
#     copy_tree(os.path.join(core_dir, 'data'), dest_folder)


def find_database_for(value: str) -> str:
    file_candidate = os.path.expanduser(value)

    if os.path.exists(file_candidate):
        app_logger.info('User selected database `{}` is available, using it'.format(file_candidate))
        return file_candidate

    #_ensure_core_version_in_user_dbs()
    user_databases_prefix = os.path.expanduser(cpdb_releases)

    if not os.path.isdir(user_databases_prefix):
        app_logger.warning(f"No local databases found. Will download latest database from remote to {os.path.expanduser(cpdb_releases)}.")
        download_database("latest")

    if value == 'latest' or not value:
        available = list_local_versions()
        latest_available = available[0]
        #app_logger.info('Latest local available version is `{}`, using it'.format(latest_available))
        value = latest_available

    downloaded_candidate = os.path.join(user_databases_prefix, value, database_file)
    valid_database = os.path.exists(downloaded_candidate)

    if not valid_database:
        #app_logger.info(f"Using database '{value}' ({downloaded_candidate})")
    #else:
        app_logger.warning("Database '{value}' is not available. Trying to download it.")
        download_database(value)
        return find_database_for(value)

    return downloaded_candidate


def collect_database(database, output_path, **kwargs):
    database_file_path = os.path.join(output_path, database)

    app = create_app(verbose=True, database_file=database_file_path, collecting=True)
    app.database_manager.database.drop_everything()
    app.database_manager.database.create_all()

    LocalCollectorLauncher(database_file_path).all(**kwargs)


def download_database(version):
    try:
        if not version or version == 'latest':
            latest_release = _latest_release()
            version = latest_release['tag']
            zip_to_download = latest_release['url']
        else:
            releases = _list_releases()
            if version not in releases:
                app_logger.error(f"Unavailable version selected. Available versions are: { ', '.join(releases.keys()) }")
                exit(1)

            zip_to_download = releases[version]['url']

        app_logger.info(f"Downloading release {version} of CellPhoneDB database")

        output_folder = os.path.expanduser('{}/{}'.format(cpdb_releases, version))
        os.makedirs(output_folder, exist_ok=True)

        zip_response = requests.get(zip_to_download)
        app_logger.info(f"Download completed!")
        app_logger.info(f"Copying database to {output_folder}")
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

        # Commenting out older behaviour for multi-versions upport
        # asumes that all versions al compatible and that all folders
        # in ~/.cpdb/releases/ are valid releases
        #core = _get_core_version()
        #local_versions = os.listdir(releases_folder)
        #compatible_versions = [version for version in local_versions if _matching_major(core, version)]

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


def list_local_database_versions():
    releases = list_local_versions()

    if not releases:
        print('There are no versions available')
        exit(1)

    for idx, version in enumerate(releases):
        note = ' *latest' if idx == 0 else ''
        print('version {}{}'.format(version, note))


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
    #core_version = _get_core_version()

    return {item['tag_name']: _format_release(item) for item in releases}


def _format_release(item: dict) -> dict:
    tag_name = item['tag_name']

    return {'url': item['zipball_url'],
            'tag': tag_name,
            'date': item['published_at'],
            'link': item['html_url'],
            'compatible': True
            }


def _matching_major(core, candidate):
    return True
    # make different versions compatible
    #return _major(LooseVersion(candidate)) == _major(core)
