import pandas as pd
import numpy as np


def set_defaults(df: pd.DataFrame, defaults: dict, quiet=False) -> pd.DataFrame:
    df = df.copy()
    for column_name, default_value in defaults.items():
        if column_name not in df:
            if not quiet:
                print('missing column in dataframe: {}, set to default {}'.format(column_name, default_value))
            df[column_name] = default_value
            continue

        df[column_name].replace({np.nan: default_value}, inplace=True)

    return df
