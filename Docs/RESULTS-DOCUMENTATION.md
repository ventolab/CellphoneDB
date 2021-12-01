CELLPHONEDB TABLE GUIDE
============================================
Please, Check https://www.cellphonedb.org/documentation for more info

Without running statistical inference of receptor-ligand interactions, only "means.csv" and "desconvoluted.csv" are generated. The "means.csv" file contains mean values for each ligand-receptor interaction. The "deconvoluted.csv" file gives additional information for each of the interacting partners. This is important as some of the interacting partners are heteromers. In other words, multiple molecules have to be expressed in the same cluster in order for the interacting partner to be functional. If the user uses the statistical inference approach (`method statistical_analysis`), additional "pvalues.csv" and "significant_means.csv" file are generated with the values for the significant interactions.

We recently introduced a novel method to query the database alternative to the statistical inference approach. The new method (`method degs_analysis`) will select interactions where all the genes are expressed by a fraction of cells above a threshold and at least one gene is a DEG. The user can identify DEGs using their preferred tool and provide the information to CellphoneDB via text file.  If the user uses the degs_analysis approach, "relevant_interactions.txt" (instead of  "pvalues.csv") and "significant_means.csv" file are generated with the values for the significant interactions. 


P-value (pvalues.txt), Mean (means.txt), Significant mean (significant_means.txt) and Relevant interactions (relevant_interactions.txt)
---------------------
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
* p.values: p-values for the all the interacting partners: p.value refers to the enrichment of the interacting ligand-receptor pair in each of the interacting pairs of cell types. (Only in pvalues.txt)
* significant_mean: Significant mean calculation for all the interacting partners. If p.value < 0.05, the value will be the mean. Alternatively, the value is set to 0. (Only in significant_means.txt)
* relevant_interactions: Indicates if the interaction is relevant (1) or not (0). If a gene in the interaction is a DEG (i.e. a gene in the DEG,tsv file), and all the participant genes are expressed, the interaction will be classified as relevant. Alternatively, the value is set to 0. ( Only in relevant_interactions.txt)

Importantly, the interactions are not symmetric. Partner A expression is considered on the first cluster, and partner B expression is considered on the second cluster. In other words:
* clusterA_clusterB = clusterA expressing partner A and clusterB expressing partner B.
* clusterA_clusterB and clusterB_clusterA  values will be different.


Deconvoluted (deconvoluted.txt)
-------------------------------
* gene_name: Gene identifier for one of the subunits that are participating in the interaction defined in “means.csv” file. The identifier will depend on the input of the user list.
* uniprot: UniProt identifier for one of the subunits that are participating in the interaction defined in “means.csv” file.
* is_complex: True if the subunit is part of a complex. Single if it is not, complex if it is.
* protein_name: Protein name for one of the subunits that are participating in the interaction defined in “means.csv” file.
* complex_name: Complex name if the subunit is part of a complex. Empty if not.
* id_cp_interaction: Unique CellPhoneDB identifier for each of the interactions stored in the database.
* mean: Mean expression of the corresponding gene in each cluster.
