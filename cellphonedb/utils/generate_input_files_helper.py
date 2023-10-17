import json
import os
from datetime import datetime
from json import JSONDecodeError
from cellphonedb.utils import file_utils

import pandas as pd


def interaction_exist(interaction, interactions, interaction_1_key='uniprot_1', interaction_2_key='uniprot_2'):
    """
    Checks if interaction already exists in first dataframe.
    """

    if len(interactions[(interactions[interaction_1_key] == interaction[interaction_1_key]) & (
            interactions[interaction_2_key] == interaction[interaction_2_key])]):
        return True

    if len(interactions[(interactions[interaction_2_key] == interaction[interaction_1_key]) & (
            interactions[interaction_1_key] == interaction[interaction_2_key])]):
        return True

    return False


def sort_interactions_partners_alphabetically(interactions: pd.DataFrame,
                                              names: tuple = ('partner_1', 'partner_2')) -> pd.DataFrame:
    """
    Sorts interactions partners alphabetically:
    ie:
    A->B        A->B
    B->C        B->C
    D->A   =>   A->D
    B->A        A->B
    E->D        D->E
    """
    interactions = interactions.copy()

    def alphabetic_sort(interaction: pd.Series):
        return pd.Series({
            names[0]: min([interaction[names[0]], interaction[names[1]]]),
            names[1]: max([interaction[names[0]], interaction[names[1]]])
        })

    interactions[[names[0], names[1]]] = interactions.apply(alphabetic_sort, axis=1)

    return interactions


def normalize_interactions(interactions, interaction_1_key='protein_1', interaction_2_key='protein_2'):
    """
    Permutes all inversed interactions:
    ie:
    A->B        A->B
    B->A        A->B
    B->A   =>   A->B
    A->B        A->B
    B->A        A->B

    """
    return sort_interactions_partners_alphabetically(interactions, (interaction_1_key, interaction_2_key))


def add_to_meta(file: str, meta_path: str):
    if not os.path.exists(meta_path):
        meta = {}
    else:
        try:
            with open(meta_path, 'r') as metafile:
                meta = json.load(metafile)
        except JSONDecodeError:
            meta = {}

    with open(meta_path, 'w') as metafile:
        meta[file] = {'date': datetime.now().strftime('%Y%m%d')}
        json.dump(meta, metafile, indent=2)


def filter_by_cellphonedb_interactor(uniprots: pd.DataFrame, interactions: pd.DataFrame) -> pd.DataFrame:
    interactions_filtered = pd.merge(interactions, uniprots, left_on=['protein_1'], right_on=['uniprot'], how='inner')

    interactions_filtered = pd.merge(interactions_filtered, uniprots, left_on=['protein_2'],
                                     right_on=['uniprot'], how='inner', suffixes=('_1', '_2'))

    interactions_filtered.rename(columns={'uniprot_1': 'id_multidata_1', 'uniprot_2': 'id_multidata_2'}, inplace=True)
    interactions_filtered.drop_duplicates(inplace=True)

    interactions_filtered = interactions_filtered[
        interactions_filtered.apply(lambda interaction: file_utils.is_cellphonedb_interactor(interaction),
                                    axis=1)]

    return interactions_filtered[interactions.columns.values]
