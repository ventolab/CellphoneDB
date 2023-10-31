import pandas as pd
import numpy as np
from typing import Tuple


def get_active_cts_tf_for_interaction(
    partner: str,
    col_name: str,
    id_cp_interaction: str,
    receptor2tfs: dict,
    separator: str,
    active_tf2cell_types: dict,
    active_interactions: pd.DataFrame
) -> list:
    ct_pair = col_name.split(separator)
    if partner in receptor2tfs:
        for tf in receptor2tfs[partner]:
            if tf in active_tf2cell_types:
                cts = active_tf2cell_types[tf]
                active_cts_for_interaction = list(set(cts) & set(ct_pair))
                partner_active = any(active_cts_for_interaction)
                # Add interaction to active_interactions_deconvoluted only if partner is active
                # and the interaction is relevant
                if partner_active and active_interactions.at[id_cp_interaction, col_name] == 1:
                    return active_cts_for_interaction, tf
    return None, None


def collect_active_interactions_across_cell_types(
    ct_pair_cols: list,
    significant_relevant_df: pd.DataFrame,
    receptor2tfs: dict,
    active_tf2cell_types: dict,
    separator: str,
    active_interactions_deconvoluted: list,
    active_interactions: pd.DataFrame
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    for col_name in ct_pair_cols:
        cols = ['id_cp_interaction', 'interacting_pair', 'partner_a', 'partner_b', 'gene_a', 'gene_b', col_name]
        for row in significant_relevant_df[~significant_relevant_df[col_name].isna()][cols].values:
            id_cp_interaction = row[0]
            interacting_pair = row[1]
            interaction_active = False
            partner_pair = [row[4], row[5]]
            for idx, partner in enumerate(partner_pair):
                if not partner:
                    # partner is a complex - pick its id from partner_a/b position
                    partner = row[idx + 2].split(":")[1]

                active_cts_for_interaction, tf = \
                    get_active_cts_tf_for_interaction(
                        partner, col_name, id_cp_interaction, receptor2tfs, separator,
                        active_tf2cell_types, active_interactions)

                interaction_active = active_cts_for_interaction is not None
                if interaction_active:
                    for ct in active_cts_for_interaction:
                        active_interactions_deconvoluted.append(
                            [id_cp_interaction, interacting_pair, row[2], row[3], row[4], row[5], tf,
                             col_name, ct])
                    break
            if not interaction_active:
                # active_interactions already contains 1 for relevant interactions and 0 for non-relevant, hence
                # the operation below is a no-op if the interactions was already irrelevant
                active_interactions.at[id_cp_interaction, col_name] = 0
    return active_interactions, active_interactions_deconvoluted


def find_active_interactions(
        significant_relevant_df: pd.DataFrame,
        receptor2tfs: dict,
        active_tf2cell_types: dict,
        separator: str) -> Tuple[pd.DataFrame, pd.DataFrame]:
    if not active_tf2cell_types:
        return pd.DataFrame, pd.DataFrame
    # Find active interactions by CellSign method:
    # Relevant/significant interaction between cell type A and cell type B is active if one of its partners is or
    # interacts with a transcription factor (TF), and the user tells us that this TF is active in either cell type A or B.
    ct_pair_cols = significant_relevant_df.filter(regex="\\{}".format(separator)).columns.tolist()
    active_interactions_deconvoluted = []
    active_interactions = significant_relevant_df.copy()
    active_interactions.set_index('id_cp_interaction', inplace=True)
    # If active_interactions was copied from significant_means (in the case of statistical analysis), we want to make it
    # look the same as relevant_interactions (from DEG analysis). This is so that active_interactions can be handled in
    # the same way to find out if an interaction is relevant/significant before recording it as active:
    # Replace all non-nan values with 1, then all nan values with 0
    for col in ct_pair_cols:
        active_interactions[col] = active_interactions[col].apply(lambda x: 1 if not np.isnan(x) else 0)
        active_interactions.astype({col: 'int'})

    active_interactions, active_interactions_deconvoluted = \
        collect_active_interactions_across_cell_types(
            ct_pair_cols,
            significant_relevant_df,
            receptor2tfs,
            active_tf2cell_types,
            separator,
            active_interactions_deconvoluted,
            active_interactions)

    # Remove from active_interactions all interactions with 0 across all ct_pair_cols (such interaction are not active)
    active_interactions = active_interactions[
        active_interactions[ct_pair_cols].apply(lambda row: row.sum() > 0, axis=1)]
    # This is to de-fragment active_interactions - in order to improve performance
    active_interactions = active_interactions.copy()
    # Drop index - to match other CellphoneDB analysis output DataFrames
    active_interactions.reset_index(drop=False, inplace=True)
    if active_interactions_deconvoluted:
        active_interactions_deconvoluted_df = pd.DataFrame(data=np.array(active_interactions_deconvoluted),
                                                           columns=['id_cp_interaction', 'interacting_pair',
                                                                    'partner_a', 'partner_b', 'gene_a', 'gene_b',
                                                                    'active_TF', 'celltype_pairs', 'active_celltype']
                                                           )
    else:
        active_interactions_deconvoluted_df = pd.DataFrame
    return active_interactions, active_interactions_deconvoluted_df
