import pandas as pd


def gene_generator(ensembl_db: pd.DataFrame,
                   uniprot_db: pd.DataFrame,
                   hla_genes: pd.DataFrame,
                   user_gene: pd.DataFrame,
                   result_columns: list) -> pd.DataFrame:
    def get_first_gene_name(gene_names: str) -> str:
        if type(gene_names) is not str:
            return ''
        gene_names = gene_names.split(' ')
        return gene_names[0]

    uniprot_db = uniprot_db.copy()
    uniprot_db.dropna(inplace=True)
    uniprot_db['gene_name'] = uniprot_db['gene_names'].apply(get_first_gene_name)
    uniprot_db.drop_duplicates(['gene_name', 'uniprot'], inplace=True)

    # Remove hla genes from uniprot
    no_hla_uniprots = uniprot_db[~uniprot_db['gene_name'].str.contains('HLA')]
    # Merge non-hla uniprot entries with ensembl on gene_name only (ignore uniprot column for now)
    cpdb_genes = no_hla_uniprots.merge(ensembl_db, how='inner', on='gene_name',
                                       sort=False, suffixes=('', '_ensembl')) \
        .drop_duplicates(['ensembl', 'uniprot', 'gene_name'])

    # Add additional unique ensembl genes if their corresponding uniprot accession is in no_hla_uniprots['uniprot']
    ensembl_db_filtered = ensembl_db.drop_duplicates()
    ensembl_db_filtered.dropna(inplace=True)
    additional_genes = ensembl_db_filtered[
        ensembl_db_filtered['uniprot'].apply(lambda uniprot: uniprot in no_hla_uniprots['uniprot'].tolist())]
    # Filter out from additional_genes genes that are already in cpdb_genes['ensembl']
    additional_genes = additional_genes[
        additional_genes['ensembl'].apply(lambda ensembl: ensembl not in cpdb_genes['ensembl'].tolist())
    ]
    cpdb_genes = pd.concat([cpdb_genes, additional_genes], ignore_index=True, sort=False)

    # Check for any duplicates in cpdb_genes['ensembl']
    duplicated_ensembl_genes = cpdb_genes[cpdb_genes['ensembl'].duplicated(keep=False)]
    # Remove duplicated ensembl genes if hgnc_symbol != gene_name
    cpdb_genes.drop(duplicated_ensembl_genes[
                        duplicated_ensembl_genes.apply(lambda gene: gene['hgnc_symbol'] != gene['gene_name'],
                                                       axis=1)].index, inplace=True)
    # Append hla_genes to cpdb_genes
    cpdb_genes = pd.concat([cpdb_genes, hla_genes], ignore_index=True, sort=False).drop_duplicates(result_columns)

    # Append user_gene to cpdb_genes, if provided
    cpdb_genes = pd.concat([cpdb_genes, user_gene], ignore_index=True, sort=False).drop_duplicates(result_columns,
                                                                                                   keep='last')
    # Check if exist any duplicated ensembl
    # print('Duplicated ensembl genes')
    # print(len(cpdb_genes[cpdb_genes['ensembl'].duplicated()]))
    # print(cpdb_genes[cpdb_genes['ensembl'].duplicated(keep=False)].to_csv(index=False))

    return cpdb_genes[result_columns]
