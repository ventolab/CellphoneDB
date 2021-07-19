import pandas as pd


def remove_not_defined_columns(data_frame: pd.DataFrame, defined_columns: list) -> pd.DataFrame:
    """
    Remove undefined columns from dataframe

    This method recieves a dataframe and a list of columns and dropps all
    the columns that are not in the given list.

    Parameters
    ----------
    data_farme: pd.DataFrame
        Original DataFrame.
    defined_columns: list
        List of columns to keep.
        
    Returns
    -------
    pd.DataFrame
        DataFrame containing only the columns specified in defined_columns.
    """
    data_frame_keys = list(data_frame.keys())

    for key in data_frame_keys:
        if key not in defined_columns:
            data_frame.drop(key, axis=1, inplace=True)

    return data_frame
