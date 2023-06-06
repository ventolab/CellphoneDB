import pandas as pd
import numpy as np
from typing import Tuple

def find_active_interactions(
        significant_relevant_df: pd.DataFrame,
        receptor2tfs: dict,
        active_tf2cell_types: dict,
        separator: str) -> Tuple[pd.DataFrame, pd.DataFrame]:
    # Find active interactions by CellSign method:
    # Relevant/significant interaction between cell type A and cell type B is active if one of its partners is or
    # interacts with a transcription factor (TF), and the user tells us that this TF is active in either cell type A or B.
    ct_pair_cols = significant_relevant_df.filter(regex="\{}".format(separator)).columns.tolist()
    active_interactions_deconvoluted = []
    active_interactions = significant_relevant_df.copy()
    active_interactions.set_index('id_cp_interaction', inplace=True)
    # Prepare active_interactions for population
    active_interactions[ct_pair_cols] = 0
    for col_name in ct_pair_cols:
        ct_pair = col_name.split(separator)
        for row in significant_relevant_df[~significant_relevant_df[col_name].isna()] \
            [['id_cp_interaction', 'interacting_pair', 'partner_a', 'partner_b', 'gene_a', 'gene_b', col_name]].values:
            id_cp_interaction = row[0]
            interacting_pair = row[1]
            interaction_active = False
            partner_pair = [row[4], row[5]]
            for idx, partner in enumerate(partner_pair):
                if not partner:
                    # partner is a complex - pick its id from partner_a/b position
                    partner = row[idx+2].split(":")[1]
                if partner in receptor2tfs:
                    for tf in receptor2tfs[partner]:
                        if tf in active_tf2cell_types:
                            cts = active_tf2cell_types[tf]
                            active_cts_for_interaction = list(set(cts) & set(ct_pair))
                            interaction_active = any(active_cts_for_interaction)
                            if interaction_active:
                                for ct in active_cts_for_interaction:
                                    active_interactions_deconvoluted.append([id_cp_interaction, interacting_pair, row[2], row[3], row[4], row[5], tf, col_name, ct])
                if interaction_active:
                    break
            if interaction_active:
                active_interactions.at[id_cp_interaction, col_name] = 1
    # Remove from active_interactions all interactions with 0 across all ct_pair_cols (such interaction are not active)
    active_interactions = active_interactions[active_interactions[ct_pair_cols].apply(lambda row: row.sum() > 0, axis=1)]
    # Drop index - to match other CellphoneDB analysis output DataFrames
    active_interactions.reset_index(drop=False, inplace=True)
    if active_interactions_deconvoluted:
        active_interactions_deconvoluted_df = pd.DataFrame(data=np.array(active_interactions_deconvoluted),
                                                           columns=['id_cp_interaction', 'interacting_pair', 'partner_a', 'partner_b', 'gene_a', 'gene_b', \
                                                                    'active_TF', 'celltype_pairs', 'active_celltype']
                                                           )
    else:
        active_interactions_deconvoluted_df = pd.DataFrame
    return active_interactions, active_interactions_deconvoluted_df

