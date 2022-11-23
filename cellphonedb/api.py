import os
from fastapi import FastAPI
from utils import db_utils, search_utils, db_releases_utils
from fastapi.middleware.cors import CORSMiddleware
import csv, io

RELEASED_VERSION="v5.0.0"
CPDB_ROOT = os.path.join(os.path.expanduser('~'),".cpdb")

interactions, genes, complex_composition, complex_expanded = \
    db_utils.get_interactions_genes_complex(CPDB_ROOT, RELEASED_VERSION)

app = FastAPI()

# See: https://fastapi.tiangolo.com/tutorial/cors/
origins = [
    "http://localhost:8000"
]
app.add_middleware(
    CORSMiddleware,
    allow_origins=origins,
)

# Return cpdb version
@app.get("/dbversion")
def read_root():
    return {"dbversion": RELEASED_VERSION}

# Search for interactions using space- or comma-separated list of terms in {tokens}
@app.get("/search/{tokens}")
def find_interactions(tokens: str):
    results, complex_name2proteins_text = search_utils.search(tokens, CPDB_ROOT, RELEASED_VERSION)
    results_csv = io.StringIO()
    writer = csv.writer(results_csv)
    writer.writerows(results)
    return {"results_html_table" : search_utils.get_html_table(results, complex_name2proteins_text),
            "results_csv" : results_csv.getvalue()}

# Autocomplete {partial_element} using data in genes table or complex names in interactions table
@app.get("/autocomplete/{partial_element}")
def autocomplete(partial_element: str):
    matches_list = search_utils.autocomplete_query(genes, interactions, partial_element)["value"].tolist()
    return { "matches" : matches_list }

# Return all gene, protein and complex ids in cpdb
@app.get("/all_ids")
def get_all_ids():
    all_ids = search_utils.return_all_identifiers(genes, interactions)["value"].tolist()
    return { "all_ids" : all_ids }

@app.get("/db_releases")
def get_db_data_releases():
    db_releases = db_releases_utils.get_remote_database_versions_html(True)
    return { "db_releases_html_table" : db_releases }
