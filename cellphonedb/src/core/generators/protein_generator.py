import pandas as pd

from cellphonedb.src.core import core_logger
from cellphonedb.src.core.generators import generator_helper


def protein_generator(uniprot_db: pd.DataFrame,
                      curated_proteins: pd.DataFrame,
                      user_protein: pd.DataFrame,
                      default_values: dict,
                      default_types: dict,
                      result_columns: list,
                      log_path: str) -> pd.DataFrame:
    merged_proteins = _merge_proteins(uniprot_db, curated_proteins, default_values, default_types, result_columns,
                                      log_path, quiet=True)

    # TODO: Add missing mandatory check
    if isinstance(user_protein, pd.DataFrame) and not user_protein.empty:
        merged_proteins = _merge_proteins(merged_proteins, user_protein, default_values, default_types, result_columns,
                                          log_path, quiet=False)

    return merged_proteins


def _merge_proteins(base_protein: pd.DataFrame,
                    additional: pd.DataFrame,
                    default_values: dict,
                    default_types: dict,
                    result_columns: list,
                    log_file: str,
                    quiet: bool = False) -> pd.DataFrame:
    additional = additional.copy()

    # Here we set defaults for uniprot & curated data
    base_protein = generator_helper.set_defaults(base_protein, default_values, quiet)
    additional = generator_helper.set_defaults(additional, default_values, quiet)

    # we will only use these columns
    additional = additional[result_columns]
    base_protein = base_protein[result_columns]

    # Type casting to ensure they are equal
    additional = additional.astype(default_types)
    base_protein = base_protein.astype(default_types)

    join_key = 'uniprot'

    merged_protein = pd.concat([base_protein, additional], ignore_index=True, sort=False).drop_duplicates()

    if not quiet and merged_protein.duplicated(join_key).any():
        core_logger.warning('There are differences between merged files: logged to {}'.format(log_file))

        log = merged_protein[merged_protein.duplicated(join_key, keep=False)].sort_values(join_key)
        log.to_csv(log_file, index=False, sep='\t')

    # keep='last' -> keep the additional (=curated) duplicate in preference to the one from base_protein (=uniprot)
    merged_protein.drop_duplicates(join_key, keep='last', inplace=True)
    return merged_protein
