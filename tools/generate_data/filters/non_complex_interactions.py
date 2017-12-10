import os

import pandas as pd

from tools.app import data_dir, output_dir, current_dir


def only_noncomplex_interactions(interactions, complexes):
    '''

    :type interactions: pd.DataFrame
    :type complexes: pd.DataFrame
    :rtype: pd.DataFrame
    '''

    proteins_in_complex = []

    for i in range(1, 5):
        proteins_in_complex = proteins_in_complex + complexes['protein_%s' % i].dropna().tolist()

    proteins_in_complex = list(set(proteins_in_complex))

    inweb_df_no_complex = interactions[interactions['protein_1'].apply(
        lambda protein: protein not in proteins_in_complex
    )]
    inweb_df_no_complex = inweb_df_no_complex[
        inweb_df_no_complex['protein_2'].apply(
            lambda protein: protein not in proteins_in_complex
        )]

    return inweb_df_no_complex