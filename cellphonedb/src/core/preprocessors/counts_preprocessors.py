import pandas as pd
import numpy as np
from cellphonedb.src.core.exceptions.ParseCountsException import ParseCountsException


def counts_preprocessor(counts: pd.DataFrame, meta: pd.DataFrame) -> pd.DataFrame:
    """
    Ensure that counts values are of type float32, and that all cells in meta exist in counts

    Parameters
    ----------
    counts: pd.DataFrame
        Counts data
    meta: pd.DataFrame
        Meta data (a mapping between cells and cell types)

    Returns
    -------
    pd.DataFrame
        counts DataFrame in which counts values are of type float32 and all cells in meta are present

    """
    if not len(counts.columns):
        raise ParseCountsException('Counts values are not decimal values', 'Incorrect file format')
    try:
        if np.any(counts.dtypes.values != np.dtype('float32')):
            counts = counts.astype(np.float32)  # type: pd.DataFrame
    except Exception:
        raise ParseCountsException

    meta.index = meta.index.astype(str)

    if np.any(~meta.index.isin(counts.columns)):
        raise ParseCountsException("Some cells in meta did not exist in counts",
                                   "Maybe incorrect file format")

    if np.any(~counts.columns.isin(meta.index)):
        counts = counts.loc[:, counts.columns.isin(meta.index)]
    return counts
