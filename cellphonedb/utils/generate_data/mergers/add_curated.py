import pandas as pd

from utils.generate_data.mergers import merge_interactions
from utils.generate_data_helper import normalize_interactions


def add_curated(interactions: pd.DataFrame, interaction_curated: pd.DataFrame) -> pd.DataFrame:
    interactions.rename(index=str, columns={'uniprot_1': 'partner_a', 'uniprot_2': 'partner_b'},
                        inplace=True)

    interactions_curated_normalized = normalize_interactions(interaction_curated, 'partner_a',
                                                             'partner_b')
    curated_duplicated = interactions_curated_normalized[
        interactions_curated_normalized.duplicated(['partner_a', 'partner_b'], keep=False)]

    if not curated_duplicated.empty:
        print(
            'WARNING - Some curated interactions are duplicated.')
        print(interaction_curated.iloc[
            curated_duplicated.sort_values(['partner_a', 'partner_b']).index.values].to_csv(
            index=False))
    interactions.rename(index=str, columns={'protein_1': 'partner_a', 'protein_2': 'partner_b'},
                        inplace=True)

    interaction_curated_result = merge_interactions.merge_iuphar_other_and_curated_interactions(interactions,
                                                                                                interaction_curated)

    return interaction_curated_result
