import pandas as pd


def remove_genes_in_file(base_genes: pd.DataFrame, genes_to_remove: pd.DataFrame) -> pd.DataFrame:
    all_genes = pd.merge(base_genes, genes_to_remove, on=['ensembl', 'gene_name', 'hgnc_symbol', 'uniprot'],
                         how='left')

    genes_filtered = all_genes[all_genes['to_keep'] != False]

    if genes_filtered.duplicated('ensembl').any():
        print('WARNING: There are some ensembl duplicated')
        print(
            genes_filtered[genes_filtered.duplicated('ensembl', keep=False)].sort_values('ensembl').to_csv(index=False))

    columns = ['gene_name', 'uniprot', 'hgnc_symbol', 'ensembl']
    return genes_filtered[columns].sort_values(columns, axis=0)
