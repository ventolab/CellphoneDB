import pandas as pd

from cellphonedb.src.core.exceptions.ProcessMetaException import ProcessMetaException


def meta_preprocessor(meta_raw: pd.DataFrame) -> pd.DataFrame:
    meta_raw.columns = map(str.lower, meta_raw.columns)
    try:
        if 'cell' in meta_raw and 'cell_type' in meta_raw:
            meta = meta_raw[['cell', 'cell_type']]
            meta.set_index('cell', inplace=True, drop=True)
            return meta

        if type(meta_raw.index) == pd.core.indexes.multi.MultiIndex:
            raise ProcessMetaException

        elif 'cell_type' in meta_raw:
            meta = meta_raw[['cell_type']]
            if type(meta_raw.index) == pd.core.indexes.range.RangeIndex:
                meta.set_index(meta_raw.iloc[:, 0], inplace=True)
                meta.index.name = 'cell'
                return meta

            if type(meta_raw.index) == pd.core.indexes.base.Index:
                meta.index.name = 'cell'
                return meta

        meta = pd.DataFrame(data={'cell_type': meta_raw.iloc[:, 1]})
        meta.set_index(meta_raw.iloc[:, 0], inplace=True)
        meta.index.name = 'cell'
        meta.index = meta.index.astype(str)
        return meta

    except:
        raise ProcessMetaException
