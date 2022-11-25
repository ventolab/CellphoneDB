import pandas as pd
import numpy as np
from cellphonedb.src.core import core_logger
from cellphonedb.src.core.generators.generator_helper import set_defaults
from cellphonedb.src.exceptions.MissingRequiredColumns import MissingRequiredColumns


def complex_generator(curated_complex: pd.DataFrame, user_complex: pd.DataFrame, log_path: str) -> pd.DataFrame:
    if (isinstance(user_complex, pd.DataFrame) and user_complex.empty) or \
            (not isinstance(user_complex, pd.DataFrame) and not user_complex):
        return curated_complex

    try:
        return _merge_complex(curated_complex, user_complex, log_path)
    except MissingRequiredColumns as e:
        core_logger.error(e)


def _merge_complex(base_complex: pd.DataFrame, additional: pd.DataFrame, log_file: str) -> pd.DataFrame:
    additional = additional.copy()

    defaults = {
        'uniprot_3': np.nan,
        'uniprot_4': np.nan,
        'receptor': False,
        'integrin': False,
        'other': False,
        'other_desc': np.nan,
        'peripheral': False,
        'receptor_desc': np.nan,
        'secreted_desc': np.nan,
        'secreted_highlight': False,
        'secreted': False,
        'transmembrane': False,
        'pdb_structure': False,
        'pdb_id': np.nan,
        'stoichiometry': np.nan,
        'comments_complex': np.nan
    }

    default_types = {
        'complex_name': str,
        'uniprot_1': str,
        'uniprot_2': str,
        'uniprot_3': str,
        'uniprot_4': str,
        'transmembrane': bool,
        'peripheral': bool,
        'secreted': bool,
        'secreted_desc': str,
        'secreted_highlight': bool,
        'receptor': bool,
        'receptor_desc': str,
        'integrin': bool,
        'other': bool,
        'other_desc': str,
        'pdb_id': str,
        'pdb_structure': str,
        'stoichiometry': str,
        'comments_complex': str,
    }

    result_columns = list(default_types.keys())

    required_columns = ['complex_name', 'uniprot_1', 'uniprot_2']

    if not set(required_columns).issubset(additional):
        raise MissingRequiredColumns(list(set(required_columns).difference(additional)))

    # TODO: Fill NA
    # Here we set defaults for additional data
    additional = set_defaults(additional, defaults)

    # we will only use these columns
    additional = additional[result_columns]
    base_complex = base_complex[result_columns]

    # Type casting to ensure they are equal
    base_complex = base_complex.astype(default_types)
    additional = additional.astype(default_types)

    join_key = 'complex_name'

    merged_complex = pd.concat([base_complex, additional], ignore_index=True, sort=False).drop_duplicates()

    if merged_complex.duplicated(join_key).any():
        core_logger.warning('There are differences between merged files: logged to {}'.format(log_file))

        log = merged_complex[merged_complex.duplicated(join_key, keep=False)].sort_values(join_key)
        log.to_csv(log_file, index=False, sep='\t')

    merged_complex.drop_duplicates(join_key, keep='last', inplace=True)
    return merged_complex
