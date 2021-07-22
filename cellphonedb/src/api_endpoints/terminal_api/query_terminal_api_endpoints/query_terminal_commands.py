from typing import Callable

import click

from cellphonedb.src.api_endpoints.terminal_api.util.choose_database import choose_database
from cellphonedb.src.app import cpdb_app
from cellphonedb.src.local_launchers.local_query_launcher import LocalQueryLauncher


def common_options(f: Callable) -> Callable:
    options = [
        click.option('--verbose/--quiet', default=True, help='Print or hide cellphonedb logs [verbose]'),
        click.option('--database', default='latest', callback=choose_database),
    ]

    for option in reversed(options):
        f = option(f)

    return f


@click.command("find_interactions_by_element")
@click.argument('element')
@common_options
@click.option('--output', default=None, help='If provided results will be saved to this path as a csv file [None]')
def find_interactions_by_element(element: str, verbose: bool, database: str, output:str = None):
    LocalQueryLauncher(cpdb_app.create_app(verbose, database)).find_interactions_by_element(element, output)


@click.command("get_interaction_gene")
@click.option('--columns', default=None, help='Columns to set in the result')
@common_options
def get_interaction_gene(columns: str, verbose: bool, database: str):
    LocalQueryLauncher(cpdb_app.create_app(verbose, database)).get_interaction_gene(columns)


@click.command("autocomplete")
@click.argument('partial_element')
@common_options
def autocomplete(partial_element: str, verbose: bool, database: str) -> None:
    LocalQueryLauncher(cpdb_app.create_app(verbose, database)).autocomplete_element(partial_element)
