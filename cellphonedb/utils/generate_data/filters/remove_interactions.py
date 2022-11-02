import pandas as pd

from utils.generate_data_helper import interaction_exist


def remove_interactions_in_file(interactions, interactions_to_remove):
    '''

    :type interactions: pd.DataFrame
    :type interactions_to_remove: pd.DataFrame
    :rtype: pd.DataFrame
    '''

    interactions_filtered = interactions[
        interactions.apply(lambda row: interaction_exist(row, interactions_to_remove), axis=1) == False]

    return interactions_filtered
