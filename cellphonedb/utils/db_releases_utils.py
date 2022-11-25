from cellphonedb.src.exceptions.NoReleasesException import NoReleasesException
from typing import Union
import requests

def get_remote_database_versions_html(include_file_browsing=False):
    try:
        releases = _github_query('releases')
        if releases:
            first_row = True
            html = "<table class=\"striped\">"
            css_style = "style=\"text-align:center\""
            for rel in releases:
                html += "<tr>"
                if first_row:
                    for header in ['Version','Release date']:
                        html += "<th {}>{}</th>".format(css_style, header)
                    if include_file_browsing:
                        html += "<th {}>{}</th>".format(css_style, 'Select file to browse')
                    html += "</tr>"
                html += "<td {}><a class=\"teal-text\" href=\"{}\">{}</a></td>".format(css_style, rel['html_url'], rel['tag_name'])
                html += "<td {}>{}</td>".format(css_style, rel['published_at'].split("T")[0])
                if include_file_browsing:
                    html += "<td {}><a class='dropdown-trigger grey lighten-1' href='#' data-target='{}_dropdown'><i class=\"material-icons teal-text\">pageview</i></a>".format(css_style, rel['tag_name'])
                    html += "<ul id='{}_dropdown' class='dropdown-content'>".format(rel['tag_name'])
                    for file_name in ["gene_input", "protein_input", "complex_input", "interaction_input"]:
                        html += "<li><a href=\"javascript:get_input_file_as_html_table(\'{}\',\'{}\');\">{}</a></li>".format(rel['tag_name'], file_name, file_name)
                    html += "</ul></td>"
                html += "</tr>"
                first_row = False
            html += "</table>"
            return html

    except NoReleasesException:
        print('WARNING: There are no versions available (or connection could not be made to server to retrieve them)')

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
