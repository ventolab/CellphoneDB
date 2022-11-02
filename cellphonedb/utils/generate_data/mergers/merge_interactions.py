import pandas as pd

from utils.generate_data_helper import normalize_interactions, interaction_exist


def merge_interactions(interactions_1, interactions_2, interaction_1_key='protein_1', interaction_2_key='protein_2'):
    """
    Merges two interactions dataframens prioritizing keeping first on duplicates
    """

    interactions_2_not_in_1 = interactions_2[
        interactions_2.apply(lambda row: interaction_exist(row, interactions_1, interaction_1_key, interaction_2_key),
                             axis=1) == False]

    interactions = interactions_1.append(interactions_2_not_in_1, sort=True)

    return interactions


def merge_iuphar_imex_interactions(iuphar_interactions: pd.DataFrame,
                                   imex_interactions: pd.DataFrame) -> pd.DataFrame():
    merge_result = merge_interactions(iuphar_interactions, imex_interactions, 'uniprot_1', 'uniprot_2')
    merge_result.fillna({'iuphar': False}, inplace=True)
    return merge_result


# TODO: Make it easier
def merge_iuphar_other_and_curated_interactions(iuphar_other_interactions: pd.DataFrame,
                                                curated_interactions: pd.DataFrame) -> pd.DataFrame:
    all_interactions = iuphar_other_interactions.append(curated_interactions)

    all_interactions.reset_index(inplace=True, drop=True)
    normalized_interactions = normalize_interactions(all_interactions, 'partner_a', 'partner_b')

    duplicated_interactions = normalized_interactions[
        normalized_interactions.duplicated(['partner_a', 'partner_b'], keep=False)]

    unique_interactions = duplicated_interactions.drop_duplicates(['partner_a', 'partner_b'])

    def merge_values(interaction: pd.Series) -> pd.Series:
        interaction = interaction.copy()
        duplicated = duplicated_interactions[
            (duplicated_interactions['partner_a'] == interaction['partner_a']) & (
                    duplicated_interactions['partner_b'] == interaction['partner_b'])]

        if not duplicated[duplicated['annotation_strategy'] == 'curated'].empty:
            interaction = duplicated[duplicated['annotation_strategy'] == 'curated'].iloc[0]

        if duplicated['iuphar'].any():
            interaction['iuphar'] = True

        return interaction

    merged_duplicated_interactions = unique_interactions.apply(merge_values, axis=1)

    non_repeated_interactions = normalized_interactions.drop_duplicates(['partner_a', 'partner_b'],
                                                                        keep=False)
    interactions_merged = non_repeated_interactions.append(merged_duplicated_interactions, sort=True, ignore_index=True)
    interactions_merged.fillna({'iuphar': False}, inplace=True)
    interactions_merged = interactions_merged.astype(
        {'iuphar': bool})

    return interactions_merged
