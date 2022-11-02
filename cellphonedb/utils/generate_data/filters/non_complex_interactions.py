import pandas as pd


def only_noncomplex_interactions(interactions, complexes):
    """
    :type interactions: pd.DataFrame
    :type complexes: pd.DataFrame
    :rtype: pd.DataFrame
    """

    proteins_in_complex = []

    for i in range(1, 5):
        proteins_in_complex = proteins_in_complex + complexes['uniprot_%s' % i].dropna().tolist()

    proteins_in_complex = list(set(proteins_in_complex))

    inweb_df_no_complex = interactions[interactions['uniprot_1'].apply(
        lambda protein: protein not in proteins_in_complex
    )]
    inweb_df_no_complex = inweb_df_no_complex[
        inweb_df_no_complex['uniprot_2'].apply(
            lambda protein: protein not in proteins_in_complex
        )]

    return inweb_df_no_complex
