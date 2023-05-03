[![Python package](https://img.shields.io/pypi/v/cellphonedb.svg?color=brightgreen&label=python-package)](https://pypi.org/project/cellphonedb)




Documentation
============================================
CellPhoneDB tool provides different methods to assess cellular crosstalk between different cell types by leveraging our CellphoneDB database of interacting molecules with single-cell transcriptome data.

### Novel features in v4
1) New python package that can be easily executed in Jupyter Notebook and Collabs. 
2) A new method to ease the query of CellPhoneDB results.
3) Tutorials to run CellPhoneDB (available [here](https://github.com/ventolab/CellphoneDB/tree/master/notebooks))
4) Improved computational efficiency of method 2 `cpdb_statistical_analysis_method`.
5) A new database ([cellphonedb-data v4.1.0](https://github.com/ventolab/cellphonedb-data)) with more manually curated interactions, making up to a total of 2,923 interactions.


# Installation

We highly recommend using an isolated python environment (as described in steps 1 and 2) using [conda](https://docs.conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html) or [virtualenv](https://docs.python.org/3/library/venv.html) but you could of course omit these steps and install via `pip` immediately.

1. Create python=>3.8 environment
   - Using conda: `conda create -n cpdb python=3.8`
   - Using virtualenv: `python -m venv cpdb`

2. Activate environment
   - Using conda: `source activate cpdb`
   - Using virtualenv: `source cpdb/bin/activate`

3. Install CellPhoneDB `pip install cellphonedb`

4. Set up the kernel for the Jupyter notebooks.
   - Install the ipython kernel: `pip install -U ipykernel`.
   - Add the environment as a jupyter kernel: `python -m ipykernel install --user --name 'cpdb'`.
   - Open/Start Jupyter and select the created kernel.

5. Download the database.
   - Follow this [tutorial](https://github.com/ventolab/CellphoneDB/blob/master/notebooks/T00_DownloadDB.ipynb).

> NOTE: Works with Python v3.8 or greater. If your default Python interpreter is for `v2.x` (you can check it with `python --version`), calls to `python`/`pip` should be substituted by `python3`/`pip3`.


# Analysis & Methods
There are three ways of running CellphoneDB, each producing a specific output:


![CellphoneDB methods](./cellphoneDB_overview.png)



- METHOD 1 simple **analysis** (>= v1): Here, no statistical analysis is performed. CellphoneDB will output the mean for all the interactions for each cell type pair combination. Note that CellphoneDB will report the means only if all the gene members of the interactions are expressed by at least a fraction of cells in a cell type (`threshold`). If the condition `threshold` is not met, the interaction will be ignored in the corresponding cell type pairs. 


- METHOD 2 **statistical_analysis** (>= v1): This is a statistical analysis that evaluates for significance all the interactions that can potentially occur in your dataset: i.e. between ALL the potential cell type pairs. Here, CellphoneDB uses empirical shuffling to calculate which ligand–receptor pairs display significant cell-type specificity. Specifically, it estimates a null distribution of the mean of the average ligand and receptor expression in the interacting clusters by randomly permuting the cluster labels of all cells. The P value for the likelihood of cell-type specificity of a given receptor–ligand complex is calculated on the basis of the proportion of the means that are as high as or higher than the actual mean. 


- METHOD 3 **degs_analysis** (>= v3): This method is proposed as an alternative to the statistical inference approach. This approach allows the user to design more complex comparisons to retrieve interactions specific to a cell type of interest. This is particularly relevant when your research question goes beyond comparing "one" cell type vs "the rest". Examples of alternative contrasts are hierarchical comparisons (e.g. you are interested in a specific lineage, such epithelial cells, and wish to identify the genes changing their expression within this lineage) or comparing disease vs control (e.g. you wish to identify upregulated genes in disease T cells by comparing them against control T cells).  For this CellphoneDB method (`cpdb_degs_analysis_method`), the user provides an input file (`test_DEGs.txt` in the command below) indicating which genes are relevant for a cell type (for example, marker genes or significantly upregulated genes resulting from a differential expression analysis (DEG)). CellphoneDB will select interactions where: 
  1) all the genes in the interaction are expressed in the corresponding cell type by more than 10% of cells (`threshold = 0.1`) and 
  2) at least one gene-cell type pair is in the provided `DEG.tsv` file. 

   The user can identify marker genes or DEGs using their preferred tool (we provide [notebooks](https://github.com/ventolab/CellphoneDB/tree/master/notebooks) for both Seurat and Scanpy users) and input the information to CellphoneDB via a [text file](https://github.com/ventolab/CellphoneDB/blob/master/README.md#preparing-your-degs-file-optional-if-method-degs_analysis). 


## METHOD 1. Retrieval of interaction expression means

With this simple `analysis` method, no analysis of significance is performed. This option will output the mean of each interaction in each cell type pair. The mean expression of a **simple interaction** is computed by averaging the expression of all the gene participants in the corresponding producing cells. To compute the mean of an interaction involving **multi-subunit heteromeric complexes** we use the member of the complex with the minimum expression.

![means](./cellphoneDB_computing_means.png)

Only interactions involving receptors and ligands expressed by more than a fraction of the cells (`threshold` default is 0.1, which is 10%) in the specific cluster are included. We generally do not consider that an interaction is feasible if one of their gene participants is expressed by less than 10% of cells (users can modify this fraction `threshold`).

   - Example command: 
   ```shell
   from cellphonedb.src.core.methods import cpdb_analysis_method

   means, deconvoluted = cpdb_analysis_method.call(
            cpdb_file_path = cellphonedb.zip,
            meta_file_path = test_meta.txt,
            counts_file_path = test_counts.h5ad,
            counts_data = 'hgnc_symbol',
            output_path = out_path)
   ```
   -  Output: Without running statistical inference of receptor-ligand interactions only `means.csv` and `deconvoluted.csv` are generated.


## METHOD 2. Statistical inference of interaction specificity

With this `statistical_analysis` method, we predict enriched receptor–ligand interactions between two cell types based on expression of a receptor by one cell type and a ligand by another cell type, using scRNA-seq data. To identify the most relevant interactions between cell types, we look for the cell-type specific interactions between ligands and receptors. 

Importantly:
1. Only receptors and ligands expressed in more than a user-specified threshold percentage of the cells in the specific cluster (`threshold` default is 0.1) are tested and will get a mean value in the significant.txt output.  
2. For the multi-subunit heteromeric complexes, we require that: 
   1. all subunits of the complex are expressed by a proportion of cells (`threshold`), and then 
   2. We use the member of the complex with the minimum expression to compute the interaction means and perform the random shuffling.

We then perform pairwise comparisons between all cell types. First, we randomly permute the cluster labels of all cells (1,000 default) and determine the mean of the average receptor expression level in a cluster and the average ligand expression level in the interacting cluster. For each receptor–ligand pair in each pairwise comparison between two cell types, this generates a null distribution. By calculating the proportion of the means which are equal or higher than the actual mean, we obtain a p-value for the likelihood of cell-type specificity of a given receptor–ligand complex. We then prioritise interactions that are highly enriched between cell types based on the number of significant pairs, so that the user can manually select biologically relevant ones.

    - Example command:
    ```shell
    from cellphonedb.src.core.methods import cpdb_statistical_analysis_method

    deconvoluted, means, pvalues, significant_means = cpdb_statistical_analysis_method.call(
            cpdb_file_path = cellphonedb.zip,
            meta_file_path = test_meta.txt,
            counts_file_path = test_counts.h5ad,
            counts_data = 'hgnc_symbol',
            output_path = out_path)
    ```
    -  Output: Apart from the outputs in method 1, additional `pvalues.csv` and `significant_means.csv` files are generated with the values for the significant interactions. In this last file, ligand–receptor pairs are ranked on the basis of their total number of significant P values across the cell populations. 

#### Cell subsampling for accelerating analyses (Optional METHOD 2)
Sc-RNA-seq datasets are growing in size exponentially as technological developments and protocol improvements enable the sequencing of more and more cells. Large-scale datasets can profile hundreds of thousands cells, which presents a challenge for the existing analysis methods in terms of both memory usage and runtime. In order to improve the speed and efficiency of our protocol and facilitate its broad accessibility, we integrated subsampling as described in Hie *et al.* 2019 (PMID: 31176620). This "geometric sketching" approach aims to maintain the transcriptomic heterogeneity within a dataset with a smaller subset of cells. The subsampling step is optional, enabling users to perform the analysis either on all cells, or with other subsampling methods of their choice.

Alternatively, the user can downsample the number of cells using their preferred method. We recommend the users use downsample their dataset to even out the contribution of each cell type (i.e. the number of cells in each cell type). This will ensure that the null distribution is representing all the cell types evenly (i.e. not biased towards cell types with larger numbers of cells). 


## METHOD 3. Retrieval of differentially expressed interactions

With this `degs_analysis` method introduced in version 3 the user can retrieve interactions where all their gene participants are expressed (in the corresponding cell type pair) and at least one gene participant is differentially expressed (list provided by the user). More specifically, this method will retrieve as **relevant** those interactions meeting both of these criteria: 
   1. all the genes in the interaction are expressed in the corresponding cell type by more than 10% of cells (`threshold = 0.1`)
   2. at least one gene-cell type pair is in the provided `DEG.tsv` file. 

The relevant/selected interactions will be labelled as 1 in the `relevant_interactions.txt` file and will get a mean assigned in the "significant_means.csv" file.  

   - Example command: 
   ```shell
   from cellphonedb.src.core.methods import cpdb_degs_analysis_method

   deconvoluted, means, relevant_interactions, significant_means = cpdb_degs_analysis_method.call(
            cpdb_file_path = cellphonedb.zip,
            meta_file_path = test_meta.txt,
            counts_file_path = test_counts.h5ad,
            degs_file_path = degs_file.txt,
            counts_data = 'hgnc_symbol',
            threshold = 0.1,
            output_path = out_path)
   ```
   -  Output: This approach will output `relevant_interactions.txt` (instead of "pvalues.csv") and the `significant_means.csv` files.  


This method gives the user the freedom to design their gene expression comparison in a way that better matches their research question. With method 2, our null hypothesis (and background distribution) considers all the cell types in the dataset and performs a "one" cell type vs "the rest" comparison. However the user may wish to use a different approach to better reflect their research scenario. Find below a list of example cases:
- The analysis needs to account for technical batch or biological covariates. Here is better to rely on differential expression approaches that can include such confounders and provide CellphoneDB the results directly.
- The user is interested on the specificities within specific lineages and wish to perform a hierarchical differential expression analysis (e.g. the user is interested in a specific lineage, such epithelial cells, and wishes to identify the genes changing their expression within this epithelial lineage; RESEARCH QUESTION: What are the interactions upregulated in epithelial-A compared to epithelial-B?).
- The user wishes to compare specific populations in a disease vs control fashion (e.g. identify upregulated genes in disease T cells by comparing them against control T cells; RESEARCH QUESTION: What are the interactions upregulated by disease T-cells?).
 
The user should perform their differential expression analysis using their preferred tool and strategy. We provide [example notebooks](https://github.com/ventolab/CellphoneDB/tree/master/notebooks) to compute DEGs for both Seurat and Scanpy users. **The user is responsible for designing a DEG analysis appropriated to the experimental design** / **research question.**

See below for how to prepare the DEGs file.


## Inclusion of Spatial Information:  Microenvironments

CellPhoneDB can prioritise interactions occurring between neighbouring cell types. The tool will restrict the cell type interacting pairs to those sharing a microenviroment (i.e. only test a combination of clusters if these coexist in a microenviroment). 

Spatial information of the cells is provided via the `microenvironments` file.  This is a two columns file indicating which cell type is in which spatial microenvironment (see [example](https://github.com/ventolab/CellphoneDB/blob/master/in/endometrium_atlas_example/endometrium_example_microenviroments.tsv) ). CellphoneDB will use this information to define possible pairs of interacting cells (i.e. pairs of clusters sharing/coexisting in a microenvironment). 

To consider microenvironments in any of the methods, add:

   ```shell
   microenvs = test_microenvs.txt
   ```
   

You can define microenvironments with prior knowledge, imaging or Visium analysis with [cell2location](https://cell2location.readthedocs.io/en/latest/notebooks/cell2location_short_demo_downstream.html#4.-Identify-groups-of-co-located-cell-types-using-matrix-factorisation).


# Input files

## Counts file

For large datasets, do not use .txt files to input counts `test_counts.txt`. Please, input counts as h5ad (recommended), h5 or a path to a folder containing a 10x output with mtx/barcode/features files. 

> NOTE that your gene/protein ids must be **HUMAN**. If you are working with another species such as mouse, we recommend you to convert the gene ids to their corresponding HUMAN orthologues.

> NOTE that by default, CellphoneDB will assume that you are using ensembl gene ids (`counts-data` ensembl as default). If you are using gene symbols, please indicate it by adding this in your run `counts-data = 'hgnc_symbol'`.

## Meta file

This is the file linking bacodes/cells to clusters/cell types. This file is generated by the users after they have annotated each cluster identified by scRNA-seq data (e.g., by using packages such as Seurat and SCANPY). The file contains two columns: ‘Cell’, indicating the name of the cell; and ‘cell_type’, indicating the name of the cluster considered. Formats accepted are .csv, .txt, .tsv, .tab and .pickle.

## DEGs file

This file is only used by is METHOD 3 `degs_analysis`. It is a .txt with two columns: the first column should be the cell type name and the second column the associated significant gene id. The remaining columns are ignored. See example [here](https://github.com/ventolab/CellphoneDB/blob/master/in/endometrium_atlas_example/endometrium_example_DEGs.tsv). 

> NOTE that CellphoneDB does not perform any filtering and all the genes in the file will be considered significant. Please ensure you filter the genes using your preferred cut-offs.  

> NOTE that the cell type/cluster name should match those in your `meta.txt`.  


## Microenvironment file

This is a .txt with two columns indicating which cell type (1st column) is in which spatial microenvironment (end column). See an example [here](https://github.com/ventolab/CellphoneDB/blob/master/in/endometrium_atlas_example/endometrium_example_microenviroments.tsv). 

CellphoneDB will use this information to restrict the pairs of interacting cell types (i.e. pairs of clusters sharing/coexisting in a microenvironment). 

You can define microenvironments with prior knowledge, imaging or Visium analysis with [cell2location](https://cell2location.readthedocs.io/en/latest/notebooks/cell2location_short_demo_downstream.html#4.-Identify-groups-of-co-located-cell-types-using-matrix-factorisation).

> NOTE that the cell type/cluster name should match those in your `meta.txt`.  

# Output files

All files (except "deconvoluted.txt") follow the same structure: rows depict interacting proteins while columns represent interacting cell type pairs. 

- The "means.txt" file contains mean values for each ligand-receptor interaction (rows) for each cell-cell interaction pair (columns). 
- The "pvalues.txt" contains the P values for the likelihood of cell-type specificity of a given receptor–ligand complex (rows) in each cell-cell interaction pair (columns), resulting from the `statistical_analysis`. 
- The "significant_means.txt" contains the mean expression (same as "means.txt") of the significant receptor–ligand complex only. This is the result of crossing "means.csv" and "pvalues.txt".
- The "relevant_interactions.txt" contains a binary matrix indicating if the interaction is relevant (1) or not (0). An interaction is classified as relevant if a gene is a DEG in a cluster/cell type (information provided by the user in the DEG.tsv file) and all the participant genes are expressed. Alternatively, the value is set to 0. This file is specific to `degs_analysis`. Each row corresponds to a ligand-receptor interaction, while each column corresponds to a cell-cell interacting pair.
- The "deconvoluted.txt" file gives additional information for each of the interacting partners. This is important as some of the interacting partners are heteromers. In other words, multiple molecules have to be expressed in the same cluster in order for the interacting partner to be functional. 


See below the meaning of each column in the outputs:

### P-value (pvalues.txt), Mean (means.txt), Significant mean (significant_means.txt) and Relevant interactions (relevant_interactions.txt)

* id_cp_interaction: Unique CellphoneDB identifier for each interaction stored in the database.
* interacting_pair: Name of the interacting pairs separated by “|”.
* partner A or B: Identifier for the first interacting partner (A) or the second (B). It could be: UniProt (prefix `simple:`) or complex (prefix `complex:`)
* gene A or B: Gene identifier for the first interacting partner (A) or the second (B). The identifier will depend on the input user list.
* secreted: True if one of the partners is secreted.
* Receptor A or B: True if the first interacting partner (A) or the second (B) is annotated as a receptor in our database.
* annotation_strategy: Curated if the interaction was annotated by the CellphoneDB developers. Otherwise, the name of the database where the interaction has been downloaded from.
* is_integrin: True if one of the partners is integrin.
* rank: Total number of significant p-values for each interaction divided by the number of cell type-cell type comparisons. (Only in significant_means.txt)
* means: Mean values for all the interacting partners: mean value refers to the total mean of the individual partner average expression values in the corresponding interacting pairs of cell types. If one of the mean values is 0, then the total mean is set to 0. (Only in means.txt)
* p.values: p-values for all the interacting partners: p.value refers to the enrichment of the interacting ligand-receptor pair in each of the interacting pairs of cell types. (Only in pvalues.txt)
* significant_mean: Significant mean calculation for all the interacting partners. If p.value < 0.05, the value will be the mean. Alternatively, the value is set to 0. (Only in significant_means.txt)
* relevant_interactions: Indicates if the interaction is relevant (1) or not (0). If a gene in the interaction is a DEG (i.e. a gene in the DEG.tsv file), and all the participant genes are expressed, the interaction will be classified as relevant. Alternatively, the value is set to 0. ( Only in relevant_interactions.txt)

Again, remember that the interactions are not symmetric. It is not the same `IL12`-`IL12 receptor` for clusterA clusterB (i.e. receptor is in clusterB) that `IL12`-`IL12 receptor` for clusterB clusterA (i.e. receptor is in clusterA).


### Deconvoluted (deconvoluted.txt)

* gene_name: Gene identifier for one of the subunits that are participating in the interaction defined in the “means.csv” file. The identifier will depend on the input of the user list.
* uniprot: UniProt identifier for one of the subunits that are participating in the interaction defined in the “means.csv” file.
* is_complex: True if the subunit is part of a complex. Single if it is not, complex if it is.
* protein_name: Protein name for one of the subunits that are participating in the interaction defined in the “means.csv” file.
* complex_name: Complex name if the subunit is part of a complex. Empty if not.
* id_cp_interaction: Unique CellphoneDB identifier for each of the interactions stored in the database.
* mean: Mean expression of the corresponding gene in each cluster.



# Interpreting the outputs


## How to read and interpret the results?

The key files are `significant_means.txt` (for statistical_analysis) or `relevant_interactions.txt` (for degs_analysis), see below. When interpreting the results, we recommend you **first define your questions of interest**. Next, focus on specific cell type pairs and manually review the interactions prioritising those with lower p-value and/or higher mean expression. For graphical representation we recommend [@zktuong](https://github.com/zktuong) repository: [ktplots](https://www.github.com/zktuong/ktplots/) in R and [ktplotspy](https://www.github.com/zktuong/ktplotspy/) in python.



CellphoneDB output is high-throughput. CellphoneDB provides all cell-cell interactions that may potentially occur in your dataset, given the expression of the cells. The size of the output may be overwhelming, but if you apply some rationale (which will depend on the design of your experiment and your biological question), you will be able to narrow it down to a few candidate interactions. The new method `degs_analysis` will allow you to perform a more tailored analysis towards specific cell-types or conditions, while the option `microenvs` will allow you to restrict the combinations of cell-type pairs to test.


It may be that not all of the cell-types of your input dataset co-appear in time and space. Cell types that do not co-appear in time and space will not interact. For example, you might have cells coming from different in vitro systems, different developmental stages or disease and control conditions. Use this prior information to restrict and ignore infeasible cell-type combinations from the outputs (i.e., columns) as well as their associated interactions (i.e. rows). You can restrict the analysis to feasible cell-type combinations using the option `microenvs`. Here you can input a two columns file indicating which cell type is in which spatiotemporal microenvironment (see [example](https://github.com/ventolab/CellphoneDB/blob/master/README.md#preparing-your-microenviroments-file-optional-if---microenvs) ). CellphoneDB will use this information to define possible pairs of interacting cells (i.e. pairs of clusters co-existing in a microenvironment) ignoring the rest of combinations.


## Why values of clusterA-clusterB are different to the values of clusterB-clusterA?

When __reading the outputs__, is IMPORTANT to note that the interactions are not symmetric. Partner A expression is considered for the first cluster/cell type (clusterA), and partner B expression is considered on the second cluster/cell type (clusterB). Thus, `IL12`-`IL12 receptor` for clusterA-clusterB (i.e. the receptor is in clusterB) is not the same that `IL12`-`IL12 receptor` for clusterB-clusterA (i.e. the receptor is in clusterA), and will have different values.

In other words:
* clusterA_clusterB = clusterA expressing partner A and clusterB expressing partner B.
* clusterA_clusterB and clusterB_clusterA  values will be different.


![CellphoneDB methods](./interpreting_results.png)


## How can I query my CellPhoneDB results?

CellPhoneDB results can be queried by making use of the `search_analysis_results` method. This method requires two of the files generated by CellPhoneDB: `significant_means` and `deconvoluted`.

Through this method, users can specify the cell pairs of interest and both; the genes `query_genes` participating in the interaction and/or the name of the interaction itself `query_interactions`. This method will search for significant/relevant interactions in which any cell specified in `query_cell_types_1` is found to any cell specified in `query_cell_types_2`. Cell pairs within any of these two lists will not be queried, that is to say, no interaction between cells A and B or C and D will be queried.

```shell
from cellphonedb.utils import search_utils

search_results = search_utils.search_analysis_results(
    query_cell_types_1 = ['cell_A', 'cell B'],  
    query_cell_types_2 = ['cell_C', 'cell D'],
    query_genes = ['Gene A', 'GeneB'], 
    query_interactions = ['interaction_name_1', 'interaction_name_2'],  
    significant_means = significant_means,
    deconvoluted = deconvoluted,
    long_format = True 
)
```



DATABASE of interactions
============================================

CellphoneDB has its own database of interactions called **CellphoneDB-data**, which can be found at https://github.com/ventolab/cellphonedb-data. CellphoneDB database is a **manually curated** repository of receptors, ligands and their interactions.  


#### Key features of CellphoneDB
- Subunit architecture is included for both ligands and receptors, representing **heteromeric complexes** accurately. This is crucial, as cell-cell communication relies on multi-subunit protein complexes that go beyond the binary representation used in most databases and studies. 
- Includes interactions involving **non-peptidic molecules** (i.e., not encoded by a gene) acting as ligands. Examples of these include steroid hormones (e.g., estrogen). To do so, we have reconstructed the biosynthetic pathways and used the last representative enzyme as a proxy of ligand abundance. We retrieve this information by manually reviewing and curating relevant literature and peer-reviewed pathway resources such as REACTOME. We include more than 200 interactions involving non-peptidic ligands!
- Only includes HUMAN interactions. 


## Database design: input files
CellphoneDB stores ligand-receptor and other types of interactions as well as other properties of the interacting partners, including their subunit architecture and gene and protein identifiers. In order to create the content of the database, four main .csv data files are required: "gene_input.csv", "protein_input.csv", " complex_input.csv" and "interaction_input.csv" (See [Figure 4](https://www.nature.com/articles/s41596-020-0292-x/figures/4) Efremova et al 2018 ).

### 1. "gene_input"
Mandatory fields: "gene_name"; "uniprot"; "hgnc_symbol" and "ensembl"
This file is crucial for establishing the link between the scRNAseq data and the interaction pairs stored at the protein level. It includes the following gene and protein identifiers: 1) gene name ("gene_name"); 2) UniProt identifier ("uniprot"); 3) HUGO nomenclature committee symbol (HGNC) ("hgnc_symbol") and iv) gene ensembl identifier (ENSG) ("ensembl"). In order to create this file, lists of linked protein and gene identifiers are downloaded from UniProt and merged using gene names. Several rules need to be considered when merging the files:

- UniProt annotation prevails over the gene Ensembl annotation when the same gene Ensembl points towards different UniProt identifiers.
- UniProt and Ensembl lists are also merged by their UniProt identifier but this information is only used when the UniProt or Ensembl identifier are missing in the original merged list by gene name.
- If the same gene name points towards different HGNC symbols, only the HGNC symbol matching the gene name annotation is considered.
- Only one HLA isoform is considered in our interaction analysis and it is stored in a manually HLA-curated list of genes, named "HLA_curated".

### 2. "protein_input"
Mandatory fields: "uniprot"; "protein_name"
Optional fields: "transmembrane"; "peripheral"; "secreted"; "secreted_desc"; "secreted_highlight"; "receptor"; "receptor_desc" ; "integrin"; "other"; "other_desc"; "tags"; "tags_ description"; "tags_reason"
Two types of input are needed to create this file: i) systematic input using UniProt annotation, and ii) manual input using curated annotation both from developers of CellphoneDB ("protein_curated") and users. For the systematic input, the UniProt identifier ("uniprot") and the name of the protein ("protein_name") are downloaded from UniProt. For the curated input, developers and users can introduce additional fields relevant to the future systematic assignment of ligand-receptor interactions (see below the "Systematic input from other databases" section for interaction_list). Importantly, the curated information always has priority over the systematic information. The optional input is organised using the fields described below:

#### a. Location of the protein in the cell
There are four non-exclusive options: transmembrane ("transmembrane"), peripheral ("peripheral") and secreted ("secreted", "secreted_desc" and "secreted_highlight").

Plasma membrane proteins are downloaded from UniProt using the keyword KW-1003 (cell membrane). Peripheral proteins from the plasma membrane are annotated using the UniProt keyword SL-9903, and the remaining proteins are annotated as transmembrane proteins. We complete our lists of plasma transmembrane proteins by doing an extensive manual curation using literature mining and UniProt description of proteins with transmembrane and immunoglobulin-like domains.

Secreted proteins are downloaded from UniProt using the keyword KW-0964 (secreted), and are further annotated as cytokines (KW-0202), hormones (KW-0372), growth factors (KW-0339) and immune-related using UniProt keywords and manual annotation. Cytokines, hormones, growth factors and other immune-related proteins are indicated as "secreted_highlight" in the protein_input lists. "secreted_desc" indicates a quality for the secreted protein.

All the manually annotated information is carefully tagged and can be identified. Please see the "curation tags" section below.

#### b. Receptors and integrins
Three fields are allocated to annotate receptors or integrins: "receptor", "receptor_desc" and "integrin".

Receptors are defined by the UniProt keyword KW-0675. The receptors list is extensively reviewed and new receptors are added based on UniProt description and bibliography revision. Receptors involved in immune-cell communication are carefully curated. For some of the receptors, a short description is included in "receptor_desc".

Integrin is a manual curation field that indicates the protein is part of the integrin family. All the annotated information is carefully tagged and can be identified. For details, see the "curation tags" section below.

#### c. Membrane and secreted proteins not considered for the cell-cell communication analysis ("others")
We created another column named "others" that consists of proteins that are excluded from our analysis. We also added "others_desc" to add a brief description of the excluded protein. Those proteins include: (i) co-receptors; (ii) nerve-specific receptors such as those related to ear-binding, olfactory receptors, taste receptors and salivary receptors; (iii) small molecule receptors; (iv) immunoglobulin chains; (v) viral and retroviral proteins, pseudogenes, cancer antigens and photoreceptors.

Curation "tags"
Three fields indicate whether the protein has been manually curated: "tags", "tags_ description" and "tags_reason".

The "tags" field is related to the manual curation of a protein and contains three fields: (i) ‘N/A’, which indicates that the protein matched with UniProt description in all fields; (ii) ‘To_add’, which indicates that secreted and/or plasma membrane protein annotation has been added; and (iii) ‘To_comment’, which indicates that the protein is either secreted (KW-0964) or membrane-associated (KW-1003), but that we manually added a specific property of the protein (that is, the protein is annotated as a receptor).

The "tags_reason" field is related to the protein properties and has five possible values: (i) ‘extracellular_add’, which indicates that the protein is manually annotated as plasma membrane; (ii) ‘peripheral_add’, which indicates that the protein is manually annotated as a peripheral protein instead of plasma membrane; (iii) ‘secreted_add’, which indicates that the protein is manually annotated as secreted; (iv) ‘secreted_high’, which indicates that the protein is manually annotated as secreted highlight. For cytokines, hormones, growth factors and other immune-related proteins, the option (v) ‘receptor_add’ indicates that the protein is manually annotated as a receptor.

Finally, the "tags_description" field is a brief description of the protein, function or property related to the manually curated protein.

### 3. "complex_input"
Mandatory fields: "complex_name"; "uniprot1, 2,..."
Optional fields: "transmembrane"; "peripheral"; "secreted"; "secreted_desc"; "secreted_highlight"; "receptor"; "receptor_desc" ; "integrin"; "other"; "other_desc"; "pdb_id"; "pdb_structure" ; "stoichiometry"; "comments_complex"
Heteromeric receptors and ligands - that is, proteins that are complexes of multiple gene products - are annotated by reviewing the literature and UniProt descriptions. Cytokine complexes, TGF family complexes and integrin complexes are carefully annotated.

These lists contain the UniProt identifiers for each of the heteromeric ligands and receptors ("uniprot1", "uniprot2",...) and a name given to the complex ("complex_name"). Common fields with "protein_input" include: "transmembrane", "peripheral", "secreted", "secreted_desc", "secreted_highlight", "receptor", "receptor_desc", "integrin", "other", "other_desc" (see description in the above "protein_input" section for clarification). In addition, we include other optional information that may be relevant for the stoichiometry of the heterodimers. If heteromers are defined in the RCSB Protein Data Bank (http://www.rcsb.org), structural information is included in our CellphoneDB annotation in the "pdb_structure", "pdb_id" and "stoichiometry". An additional field "comments_complex" was created to add a brief description of the heteromer.

### 4. "interaction_input"
Mandatory fields: "partner_a"; "partner_b"; "annotation_strategy"; "source"
Optional fields: "protein_name_a"; "protein_name_b"
Interactions stored in CellphoneDB are annotated using their UniProt identifier (binary interactions) or the name of the complex (interactions involving heteromers) ("partner_a" and "partner_b"). The name of the protein is also included, yet not mandatory ("protein_name_a" and "protein_name_b"). Protein names are not stored in the database.

There are two main inputs of interactions: i) a systematic input querying other databases, and ii) a manual input using curated information from CellphoneDB developers ("interaction_curated") and users. The method used to assign the interaction is indicated in the "annotation_strategy" column.

Each interaction stored has a CellphoneDB unique identifier ("id_cp_interaction") generated automatically by the internal pipeline.


## User-defined database
Our system allows users to create their own database of interactions and complexes. In order to do so, the format of the users’ lists must be compatible with the input files. 

Do you want to contribute  our curation effort?
Users can submit their lists using the Python package version of CellphoneDB, and then send them via <contact@cellphonedb.org> or a pull request to the CellphoneDB data [repository](https://github.com/ventolab/CellphoneDB-data).


Plotting results
============================================

Currently CellPhoneDB relies on external plotting implementations to represent the results. Examples are provided in the [tutorials](notebooks).

We recommend using tools such as the ktplots:
[@zktuong](https://github.com/zktuong):
- [ktplots](https://www.github.com/zktuong/ktplots/) (R; preferred)
- [ktplotspy](https://www.github.com/zktuong/ktplotspy/) (python; under development)


Release notes
============================================

### cellphonedb v4.1.0
1) New python package that can be easily executed in Jupyter Notebook and Collabs. 
2) A new method to ease the query of CellPhoneDB results.
3) Tutorials to run CellPhoneDB (available [here](https://github.com/ventolab/CellphoneDB/tree/master/notebooks))
4) Improved computational efficiency of method 2 `cpdb_statistical_analysis_method`.

### cellphonedb-data v4.1.0
1) New database (v4.1.0) with more manually curated interactions, making up to a total of 2,923 interactions.
2) Non-curated exernal databases are discarded.

### cellphonedb-data v4.0.0
1) More manually curated interactions added, with special focus on protein acting as heteromeric complexes. This version fo the database includes almost 2,000 high-confidence interactions, including heteromeric complexes! We believe modelling complexes is key to minimise false positives in the predictions.
2) Includes interactions involving non-peptidic molecules (i.e., not encoded by a gene) acting as ligands. Examples of these include steroid hormones (e.g., estrogen). To do so, we have reconstructed the biosynthetic pathways and used the last representative enzyme as a proxy of ligand abundance. We retrieve this information by manually reviewing and curating relevant literature and peer-reviewed pathway resources such as REACTOME. We include more than 200 interactions involving non-peptidic ligands!

### cellphonedb v3.0.0
1) New method to incorporate spatial information. CellPhoneDB now allows the incorporation of spatial information of the cells via the `microenvironments` file. CellphoneDB will use this information to define possible pairs of interacting cells (i.e. pairs of clusters sharing/coexisting in a microenvironment). You can define microenvironments with prior knowledge, imaging or Visium analysis with [cell2location](https://cell2location.readthedocs.io/en/latest/notebooks/cell2location_short_demo_downstream.html#4.-Identify-groups-of-co-located-cell-types-using-matrix-factorisation).
2) New DEG analysis method. This method relies on Differentially Expressed Genes (`cpdb_degs_analysis_method` or `cellphonedb method degs_analysis`) as an alternative to the permutation-based approach. The user identifies the DEGs using their preferred tool and provides the information to CellphoneDB via text file. 

### cellphonedb-data v3.0.0
1) Updated interactions involving WNT pathway.


Tutorials
============================================
Tutorials to run CellPhoneDB are available [here](https://github.com/ventolab/CellphoneDB/tree/master/notebooks)



FAQs
============================================
### 1. What are the counts input files accepted? 
CellphoneDB accepts counts files in the following formats: as a text file (with columns indicating individual cells and rows indicating genes), as a h5ad (recommended), a h5 or a path to a folder containing a 10x output with mtx/barcode/features files.

### 2. How to extract the CellphoneDB input files from a Seurat object? 
We recommend using normalised count data. This can be obtained by taking the normalised slot from the Seurat object or by taking the raw data slot and applying the normalisation manually. The user can also normalise the data using their preferred method.


```R
# R
# take raw data and normalise it
count_raw <- seurat_obj@assays$RNA@counts[,seurat_obj@cell.names]
count_norm <- apply(count_raw, 2, function(x) (x/sum(x))*10000)
write.table(count_norm, ‘CellphoneDB_count.txt’, sep=’\t’, quote=F)
```                        

Note that you can also export your Seurat object as a 10x .mtx output format containing the matrix.mtx.gz, features.tsv.gz and barcodes.tsv.gz. To generate these you can use:

```R
# R
library(Matrix)
writeMM(obj = seurat_obj@assays$RNA@counts, file = 'outdir/matrix.mtx') # or the normalised counts stored in seurat_obj@assays$RNA@data
features <- as.data.frame(rowData(seurat_obj))
features$gene_names<- rownames(features)
features <- features[c('gene_ids','gene_names','feature_types')]    
write.table(features, file = 'outdir/features.tsv', sep = '\t', quote=F, row.names = F, col.names = F)  
write.table(colData(seurat_obj), file = 'outdir/barcodes.tsv', sep = '\t', quote=F, row.names = F, col.names = F)
# and then compress the files to get .gz
```
                    
### 3. How to extract the CellphoneDB input files from a scanpy anndata? 

You can provide an anndata as .h5ad file.
                   
                    
### 4. Should the input file with the count data be with HGNC symbols (gene names) or Ensembl IDs? 

CellphoneDB.2 allows the use of both HGNC symbols and Ensembl IDs.

Please, specify HGNC symbols with  `counts-data = hgnc_symbol`.


### 5. What is the purpose of subsampling? 
The datasets that are generated are increasing in the number of sequenced cells exponentially. In order to increase the speed of CellphoneDB, we included an optional step in the analysis the method described in (Hie B, Cho H, DeMeo B, Bryson B and Berger B, Geometric Sketching Compactly Summarizes the Single-Cell Transcriptomic Landscape, Cell Systems 2019). The user can also choose another method to subsample and simply input the subsampled data into CellphoneDB in the same way as described above. We recommend subsampling for very big datasets; the minimum number of cells to use the subsampling option is 1000.

### 6. What is the meaning of “Rank” in the “significant_means.txt” output file? 
The rank is calculated by counting the significant p-values per interaction pair (per row) and dividing with the total number of cluster-cluster comparisons. The idea is to prioritise interactions that are highly specific, that is they have only one or few significant p-values and to have on the bottom of the list the interactions that are present everywhere or not present anywhere at all.


Citing
============================================

The first version of CellphoneDB was originally developed at the [Teichmann Lab](http://www.teichlab.org/) in the Wellcome Sanger Institute (Cambridge, UK) by Roser Vento-Tormo and Mirjana Efremova. Currently, it is being further developed and supported by the [Vento-Tormo Lab](https://ventolab.org/) (CellphoneDB ≥v3).

If you use CellphoneDB or CellphoneDB-data, please cite our papers:

- **CellPhoneDB v1 (original)**: Single-cell reconstruction of the early maternal-fetal interface in humans. Vento-Tormo R, Efremova M, et al., Nature. 2018 [link](https://www.nature.com/articles/s41586-018-0698-6)

- **CellPhoneDB v2**: Inferring cell-cell communication from combined expression of multi-subunit receptor-ligand complexes. Efremova M, Vento-Tormo M, Teichmann S, Vento-Tormo R. Nat Protoc. 2020 [link](https://www.nature.com/articles/s41596-020-0292-x)

- **CellphoneDB v3**: Mapping the temporal and spatial dynamics of the human endometrium in vivo and in vitro. L Garcia-Alonso, L-François Handfield, K Roberts, K Nikolakopoulou et al. Nature Genetics 2021 [link](https://www.nature.com/articles/s41588-021-00972-2)

- **CellphoneDB v4 (latest)**: Single-cell roadmap of human gonadal development. L Garcia-Alonso, V Lorenzi et al. 2022 Nature [link](https://www.nature.com/articles/s41586-022-04918-4)
