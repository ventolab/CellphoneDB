import pandas as pd

from utils.generate_data_helper import normalize_interactions


def call(iuphar_guidepharmacology: pd.DataFrame, genes: pd.DataFrame, proteins: pd.DataFrame) -> pd.DataFrame:
    iuphar_filtered = iuphar_guidepharmacology[iuphar_guidepharmacology['target_species'] == 'Human']

    iuphar_filtered = iuphar_filtered[iuphar_filtered['ligand_species'] == 'Human']

    iuphar_filtered = _process_target_uniprot(iuphar_filtered, proteins)

    iuphar_filtered = _process_ligand_gene_symbol(genes, iuphar_filtered)
    iuphar_filtered_cellphone_format = pd.DataFrame()
    iuphar_filtered_cellphone_format[['uniprot_1', 'uniprot_2']] = iuphar_filtered[
        ['uniprot', 'target_uniprot']]

    iuphar_filtered_cellphone_format['annotation_strategy'] = 'guidetopharmacology.org'
    iuphar_filtered_cellphone_format['iuphar'] = True

    iuphar_procesed = _drop_duplicates(iuphar_filtered_cellphone_format)

    return iuphar_procesed


def _drop_duplicates(iuphar_filtered_cellphone_format: pd.DataFrame) -> pd.DataFrame:
    interactions_normalized = normalize_interactions(iuphar_filtered_cellphone_format, 'uniprot_1',
                                                     'uniprot_2')
    interactions_duplicated = interactions_normalized[
        interactions_normalized.duplicated(['uniprot_1', 'uniprot_2'])]
    # if not interactions_duplicated.empty:
    #     print('WARNING: SOME IUPHAR INTERACTIONS ARE DUPLICATED')
    #     print(interactions_duplicated.sort_values(['uniprot_1', 'uniprot_2']).to_csv(index=False))
    iuphar_procesed = interactions_normalized.drop_duplicates(['uniprot_1', 'uniprot_2'])
    return iuphar_procesed


def _process_ligand_gene_symbol(genes: pd.DataFrame, iuphar_filtered: pd.DataFrame) -> pd.DataFrame:
    iuphar_filtered.dropna(subset=['ligand_gene_symbol'], inplace=True)
    iuphar_filtered = iuphar_filtered[iuphar_filtered['ligand_gene_symbol'].apply(lambda uniprot: '|' not in uniprot)]
    missing_genes = iuphar_filtered[~iuphar_filtered['ligand_gene_symbol'].isin(genes['gene_name'])][
        'ligand_gene_symbol'].drop_duplicates()
    if not missing_genes.empty:
        print('WARNING: {} GENE NAMES IN GUIDETOPHARMACOLOGY DATABASE DIDNT EXIST IN CELLPHONEDB DATABASE'.format(
            len(missing_genes)))
        print(missing_genes.tolist())
    iuphar_filtered = pd.merge(iuphar_filtered, genes, left_on='ligand_gene_symbol', right_on='gene_name')
    return iuphar_filtered


def _process_target_uniprot(iuphar_filtered: pd.DataFrame, proteins: pd.DataFrame) -> pd.DataFrame:
    iuphar_filtered.dropna(subset=['target_uniprot'], inplace=True)
    iuphar_filtered = iuphar_filtered[iuphar_filtered['target_uniprot'].apply(lambda uniprot: '|' not in uniprot)]
    missing_uniprots = iuphar_filtered[~iuphar_filtered['target_uniprot'].isin(proteins['uniprot'])][
        'target_uniprot'].drop_duplicates()
    if not missing_uniprots.empty:
        print('WARNING: {} UNIPROTS IN GUIDETOPHARMACOLOGY DATABASE DIDNT EXIST IN CELLPHONEDB DATABASE'.format(
            len(missing_uniprots)))
        print(missing_uniprots.tolist())
    iuphar_filtered = iuphar_filtered[iuphar_filtered['target_uniprot'].isin(proteins['uniprot'])]
    return iuphar_filtered
