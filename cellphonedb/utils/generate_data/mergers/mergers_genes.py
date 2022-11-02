import pandas as pd


def merge_genes_from_uniprot_ensembl_db(ensembls: pd.DataFrame, proteins: pd.DataFrame,
                                        uniprots: pd.DataFrame) -> pd.DataFrame:
    uniprots_filtered = merge_genes_cellphone(proteins, uniprots)
    merged_genes = pd.merge(uniprots_filtered, ensembls, left_on=0, right_on='Gene name')

    result = _merge_ensembl_uniprots(merged_genes)
    result = result[~result['gene_name'].str.contains('HLA')]
    result = result[result['uniprot'].apply(lambda uniprot: uniprot in proteins['uniprot'].tolist())]
    result = result.sort_values(by=list(result.columns.values))

    return result


def _merge_ensembl_uniprots(merged_genes: pd.DataFrame) -> pd.DataFrame:
    merged_from_uniprot = pd.DataFrame()

    merged_from_uniprot[['ensembl', 'uniprot', 'gene_name', 'hgnc_symbol']] = merged_genes[
        ['Gene stable ID', 'Entry', 'Gene name', 'HGNC symbol']]

    merged_from_ensembl = pd.DataFrame()
    merged_from_ensembl[['ensembl', 'uniprot', 'gene_name', 'hgnc_symbol']] = merged_genes[[
        'Gene stable ID', 'UniProtKB/Swiss-Prot ID', 'Gene name', 'HGNC symbol']]

    result = merged_from_uniprot.append(merged_from_ensembl, sort=False).reset_index(drop=True)

    result.dropna(subset=['ensembl', 'uniprot', 'gene_name'], inplace=True)

    check_empty_hgnc(result)

    result['hgnc_symbol'] = result.apply(
        lambda gene: gene['hgnc_symbol'] if gene['hgnc_symbol'] == '' else gene['gene_name'], axis=1)

    result = result.drop_duplicates().reset_index(drop=True)
    return result


def check_empty_hgnc(genes: pd.DataFrame) -> None:
    empty_hgnc = genes[pd.isnull(genes['hgnc_symbol'])].drop_duplicates()
    if not empty_hgnc.empty:
        print('SOME GENES DO NOT HAVE HGNC')
        print(empty_hgnc.to_csv(index=False))


def merge_genes_cellphone(proteins: pd.DataFrame, uniprots: pd.DataFrame) -> pd.DataFrame:
    uniprots_merged = pd.merge(proteins, uniprots, left_on='uniprot', right_on='Entry', indicator=True, how='outer')

    only_in_protein = uniprots_merged[uniprots_merged['_merge'] == 'left_only']

    if not only_in_protein.empty:
        print('SOME PROTEINS DIDNT EXIST IN UNIPROT GENE DATABASE')
        print(only_in_protein['uniprot'].tolist())

    uniprots_filtered = uniprots_merged[uniprots_merged['_merge'] == 'both']
    uniprots_filtered = _deconvolute_genenames(uniprots_filtered)

    return uniprots_filtered


def _deconvolute_genenames(uniprots_filtered: pd.DataFrame) -> pd.DataFrame:
    '''
    Deconvolutes gene names in uniprots database
    ie:
    input:
    Entry   Gene names
    A       a b c
    B       d

    result:

    Entry   Gene name
    A       a
    A       b
    A       c
    B       d
    '''
    uniprots_columns = ['Entry', 'Gene names', 'Ensembl transcript']
    uniprots_filtered = uniprots_filtered[uniprots_columns]

    deconvoluted = uniprots_filtered['Gene names'].str.split(' ').apply(pd.Series)
    deconvoluted.index = uniprots_filtered.set_index(uniprots_columns).index
    deconvoluted = deconvoluted.stack().reset_index(uniprots_columns).reset_index(drop=True)

    return deconvoluted
#
