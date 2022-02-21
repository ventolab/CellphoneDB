CELLPHONEDB GUIDE
============================================

## Three analysis types in CellPhoneDB
There are three ways of running cellphoneDB, each producing a specific output:

- **statistical_analysis** (>= v1): This is a statistical analsyis to retrieve all the interactions that can potentially occur in your dataset between ALL cell type pairs. Here, CellphoneDB uses empirical shuffling to calculate which ligand–receptor pairs display significant cell-type specificity. Specifically, it estimates a null distribution of the mean of the average ligand and receptor expression in the interacting clusters by randomly permuting the cluster labels of all cells. The P value for the likelihood of cell-type specificity of a given receptor–ligand complex is calculated on the basis of the proportion of the means that are as high as or higher than the actual mean. 
    - Example command: 
    ```shell
    cellphonedb method statistical_analysis test_meta.txt test_counts.txt
    ```
    -  Output: If the user uses the statistical inference approach (`method statistical_analysis`), additional "pvalues.csv" and "significant_means.csv" file are generated with the values for the significant interactions. Finally, ligand–receptor pairs are ranked on the basis of their total number of significant P values across the cell populations. 

- **degs_analysis** (>= v3): We recently introduced a novel __method to query the database__ alternative to the statistical inference approach. This approach allows the user to design more complex comparisons to define genes specific to a cell type. This is particularly relevant when comparison "one" vs "rest" does not apply. Examples are hierarchical comparisons (comparing cell states within a lineage) or comparing disease vs control.  Here, the user provides an input file (`test_DEGs.txt` in the command below) indicating which genes are relevant for a cell type (for example, marker genes or upregulated genes resulting from a differential expression analysis (DEG)).  The new CellphoneDB method (`method degs_analysis`) will select interactions where: (i) all the genes are expressed by a fraction of cells above a threshold and (ii) at least one gene-cell type pair is in the provided DEG.tsv file. The user can identify marker genes or DEGs using their preferred tool (we provide [notebooks](https://github.com/ventolab/CellphoneDB/tree/master/notebooks) for both Seurat and Scanpy users) and provide the information to CellphoneDB via a text file. 
    - Example command: 
    ```shell
    cellphonedb method degs_analysis test_meta.txt test_counts.txt test_DEGs.txt
    ```
    -  Output: If the user uses the `degs_analysis` approach, "relevant_interactions.txt" (instead of  "pvalues.csv") and "significant_means.csv" files are generated.  

- **analysis** (>= v1): Here, no statistical analysis is performed. CellphoneDB will output all the interactions where all the gene members are expressed in above a fraction of cells (`--threshold`). OUTPUT: Without running statistical inference of receptor-ligand interactions, only "means.csv" and "desconvoluted.csv" are generated. 
    - Example command: 
    ```shell
    cellphonedb method analysis test_meta.txt test_counts.h5ad
    ```
    -  Output: If the user uses the statistical inference approach (`method statistical_analysis`), additional "pvalues.csv" and "significant_means.csv" file are generated with the values for the significant interactions. Finally, ligand–receptor pairs are ranked on the basis of their total number of significant P values across the cell populations. 

For large datasets, do not use .txt files for counts `test_counts.txt`. Input counts as h5ad (recommended), h5 or a path to a folder containing a 10x output with mtx/barcode/features files. NOTE that your gene/protein ids must be HUMAN. If you are working with another specie such as mouse, we recommend you to convert the gene ids to their corresponding orthologous.

Please, Check https://www.cellphonedb.org/documentation for more info


## Interpreting the outputs

### Why values of clusterA-clusterB are different to the values of clusterB-clusterA?

When __reading the outputs__, is IMPORTANT to note that the interactions are not symmetric. Partner A expression is considered for the first cluster/cell type, and partner B expression is considered on the second cluster/cell type. In other words:
* clusterA_clusterB = clusterA expressing partner A and clusterB expressing partner B.
* clusterA_clusterB and clusterB_clusterA  values will be different.


### Output files

All files (except "deconvoluted.txt") follow the same structure: rows depict interacting proteins while columns represent interacting cell type pairs. 

- The "means.txt" file contains mean values for each ligand-receptor interaction (rows) for each cell-cell interacting pair (columns). 
- The "pvalues.txt" contains the P values for the likelihood of cell-type specificity of a given receptor–ligand complex (rows) in each cell-cell interacting pair (columns), resulting from the `statistical_analysis`. 
- The "significant_means.txt" contains the mean expression (same as "means.txt") of the significant receptor–ligand complex, only. This is the result of crossing "means.csv" and "pvalues.txt".
- The "relevant_interactions.txt" contains a binary matrix indicating if the interaction is relevant (1) or not (0). An interaction is relevant if a gene is a DEG in a cluster/cell type (information from the user provided in the DEG.tsv file) and all the participant genes are expressed, the interaction will be classified as relevant. Alternatively, the value is set to 0. This file is specific to `degs_analysis`. Each row corresponds to a ligand-receptor interaction, while each column corresponds to a cell-cell interacting pair.
- The "deconvoluted.txt" file gives additional information for each of the interacting partners. This is important as some of the interacting partners are heteromers. In other words, multiple molecules have to be expressed in the same cluster in order for the interacting partner to be functional. 


See below the meaning of each column in the outputs:

#### P-value (pvalues.txt), Mean (means.txt), Significant mean (significant_means.txt) and Relevant interactions (relevant_interactions.txt)

* id_cp_interaction: Unique CellPhoneDB identifier for each interaction stored in the database.
* interacting_pair: Name of the interacting pairs separated by “|”.
* partner A or B: Identifier for the first interacting partner (A) or the second (B). It could be: UniProt (prefix `simple:`) or complex (prefix `complex:`)
* gene A or B: Gene identifier for the first interacting partner (A) or the second (B). The identifier will depend on the input user list.
* secreted: True if one of the partners is secreted.
* Receptor A or B: True if the first interacting partner (A) or the second (B) is annotated as a receptor in our database.
* annotation_strategy: Curated if the interaction was annotated by the CellPhoneDB developers. Otherwise, the name of the database where the interaction has been downloaded from.
* is_integrin: True if one of the partners is integrin.
* rank: Total number of significant p-values for each interaction divided by the number of cell type-cell type comparisons. (Only in significant_means.txt)
* means: Mean values for all the interacting partners: mean value refers to the total mean of the individual partner average expression values in the corresponding interacting pairs of cell types. If one of the mean values is 0, then the total mean is set to 0. (Only in mean.txt)
* p.values: p-values for all the interacting partners: p.value refers to the enrichment of the interacting ligand-receptor pair in each of the interacting pairs of cell types. (Only in pvalues.txt)
* significant_mean: Significant mean calculation for all the interacting partners. If p.value < 0.05, the value will be the mean. Alternatively, the value is set to 0. (Only in significant_means.txt)
* relevant_interactions: Indicates if the interaction is relevant (1) or not (0). If a gene in the interaction is a DEG (i.e. a gene in the DEG.tsv file), and all the participant genes are expressed, the interaction will be classified as relevant. Alternatively, the value is set to 0. ( Only in relevant_interactions.txt)

Importantly, the interactions are not symmetric. Partner A expression is considered on the first cluster, and partner B expression is considered on the second cluster. In other words:
* clusterA_clusterB = clusterA expressing partner A and clusterB expressing partner B.
* clusterA_clusterB and clusterB_clusterA  values will be different.


#### Deconvoluted (deconvoluted.txt)

* gene_name: Gene identifier for one of the subunits that are participating in the interaction defined in the “means.csv” file. The identifier will depend on the input of the user list.
* uniprot: UniProt identifier for one of the subunits that are participating in the interaction defined in the “means.csv” file.
* is_complex: True if the subunit is part of a complex. Single if it is not, complex if it is.
* protein_name: Protein name for one of the subunits that are participating in the interaction defined in the “means.csv” file.
* complex_name: Complex name if the subunit is part of a complex. Empty if not.
* id_cp_interaction: Unique CellPhoneDB identifier for each of the interactions stored in the database.
* mean: Mean expression of the corresponding gene in each cluster.


#### How to read and interpret the results?

CellphoneDB output is high-throughput. CellphoneDB provides all cell-cell interactions that may potentially occur in your dataset, given the expression of the cells in your dataset. The size of the output may be overwhelming, but if you apply some rationale, you will be able to narrow it down to a few candidate interactions. 

When interpreting the results, we recommend you **first define your questions of interest**. Next, focus on specific cell type pairs and manually review the interactions prioritising those with lower p-value and/or higher mean expression. Then, select the cell type pairs and proteins of interest to generate the heatmap plots for a visual representation.

It may be that not all of the cell types of your input dataset co-appear in time and space. Cell types that do not co-appear in time and space in your individuals will not interact. Examples are cells specific to different in vitro systems, different developmental stages or disease vs control comparisons. You might wish to ignore unfeasible cell type combinations from the outputs (i.e., columns) as well as their associated interactions (i.e. rows).





