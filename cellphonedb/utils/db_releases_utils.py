from cellphonedb.src.exceptions.NoReleasesException import NoReleasesException
from typing import Union
import requests


def get_remote_database_versions_html(include_file_browsing: bool = False, min_version: float = 4.1):
    """Retrieve a html table containing CellphoneDB database versions and release dates.

        Parameters
        ----------
        include_file_browsing: bool.
            False by default; True when invoked by cellphonedb.org \
            website only - to allow the user to select and browse individual input files from \
            https://github.com/ventolab/cellphonedb-data.
        min_version: float
            If not None, the only versions returned should be >= min_version

        Returns
        -------
        str
            A string containing html table containing database versions and release dates
    """

    result = {'error': None}
    try:
        releases = _github_query('releases')
        if releases:
            first_row = True
            html = "<table class=\"striped\">"
            css_style = "style=\"text-align:center\""
            for rel in releases:
                html += "<tr>"
                if first_row:
                    for header in ['Version', 'Release date']:
                        html += "<th {}>{}</th>".format(css_style, header)
                    if include_file_browsing:
                        html += "<th {}>{}</th>".format(css_style, 'Select file to browse')
                    html += "</tr>"
                    first_row = False
                if min_version:
                    # Skip unless rel_version >= min_version
                    rel_tag = rel['tag_name']
                    rel_version = float('.'.join(rel_tag.replace('v', '').split(".")[0:2]))
                    if rel_version < min_version:
                        continue
                html += get_release_info(rel, rel_version, css_style, include_file_browsing)
                html += "</tr>"
            html += "</table>"
            result['db_releases_html_table'] = html
            return result

    except NoReleasesException:
        err_msg = 'ERROR: There are no versions available (or connection could not be made to server to retrieve them)'
        print(err_msg)
        result['error'] = err_msg
        result['db_releases_html_table'] = None
        return result


def get_release_info(
    rel: dict,
    rel_version: float,
    css_style: str,
    include_file_browsing: bool
) -> str:
    html = "<td {}><a class=\"teal-text\" href=\"{}\">{}</a></td>" \
        .format(css_style, rel['html_url'], rel['tag_name'])
    html += "<td {}>{}</td>".format(css_style, rel['published_at'].split("T")[0])
    if include_file_browsing:
        html += ("<td {}><a class='dropdown-trigger grey lighten-1' href='#' data-target='{}_dropdown'>" +
                 "<i class=\"material-icons teal-text\">pageview</i></a>").format(css_style, rel['tag_name'])
        html += "<ul id='{}_dropdown' class='dropdown-content'>".format(rel['tag_name'])
        input_files = ["gene_input", "protein_input", "complex_input", "interaction_input"]
        if rel_version >= 5.0:
            input_files.append("transcription_factor_input")
        for file_name in input_files:
            html += "<li><a href=\"javascript:get_input_file_as_html_table(\'{}\',\'{}\');\">{}</a></li>" \
                .format(rel['tag_name'], file_name, file_name)
        html += "</ul></td>"
    return html


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
