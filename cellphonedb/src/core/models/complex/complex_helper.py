import pandas as pd


def get_involved_complex_composition_from_protein(proteins: pd.DataFrame,
                                                  complex_compositions: pd.DataFrame) -> pd.DataFrame:
    """
    Returns a table of complex with all proteins declared in proteins.
    """
    multidata_protein = list(proteins.index)

    complex_counts_filtered = complex_compositions[
        complex_compositions['protein_multidata_id'].apply(lambda protein_multidata: protein_multidata in multidata_protein)]

    if complex_counts_filtered.empty:
        return pd.DataFrame(columns=complex_compositions.columns)

    def all_protein_involved(complex: pd.Series) -> bool:
        number_proteins_in_counts = len(
            complex_counts_filtered[
                complex_counts_filtered['complex_multidata_id'] == complex['complex_multidata_id']])

        if number_proteins_in_counts < complex['total_protein']:
            return False

        return True

    complex_composition_complete = complex_counts_filtered[
        complex_counts_filtered.apply(all_protein_involved, axis=1)]

    return complex_composition_complete
