CELLPHONEDB GUIDE
============================================

CellPhoneDB tool provides different methods to assess cellular crosstalk between different cell types by levereging our CellPhoneDB database of interacting molecules with single-cell transcriptome data. 


- [Analysis types in CellPhoneDB](https://github.com/ventolab/CellphoneDB/blob/master/Docs/RESULTS-DOCUMENTATION.md#analysis-types-in-cellphonedb)
- [Interpreting the outputs](https://github.com/ventolab/CellphoneDB/blob/master/Docs/RESULTS-DOCUMENTATION.md#interpreting-the-outputs)
   - [How to read and interpret the results?](https://github.com/ventolab/CellphoneDB/blob/master/Docs/RESULTS-DOCUMENTATION.md#how-to-read-and-interpret-the-results)
   - [Why values of clusterA-clusterB are different to the values of clusterB-clusterA?](https://github.com/ventolab/CellphoneDB/blob/master/Docs/RESULTS-DOCUMENTATION.md#why-values-of-clustera-clusterb-are-different-to-the-values-of-clusterb-clustera)
- [Output files](https://github.com/ventolab/CellphoneDB/blob/master/Docs/RESULTS-DOCUMENTATION.md#output-files)
- [Database design and generation](https://github.com/ventolab/CellphoneDB/blob/master/Docs/RESULTS-DOCUMENTATION.md#Database-design-and-generation)



# Analysis types in CellPhoneDB
There are three ways of running cellphoneDB, each producing a specific output:


![cellphoneDB methods](https://github.com/ventolab/CellphoneDB/blob/master/Docs/cellphoneDB_overview.png)


- METHOD 1 simple **analysis** (>= v1): Here, no statistical analysis is performed. CellphoneDB will output the mean for all the interactions for each cell type pari combination. Note that CellphoneDB will report the means only if all the gene members of the interactions are expressed by at least a fraction of cells in a celltype (`--threshold`). If the condition `--threshold` is not meet, the interaction will be ignores in the corresponding celltype pairs.
- 
   - Example command: 
   ```shell
   cellphonedb method analysis test_meta.txt test_counts.h5ad
   ```
   -  Output: Without running statistical inference of receptor-ligand interactions, only "means.csv" and "desconvoluted.csv" are generated. 

- METHOD 2 **statistical_analysis** (>= v1): This is a statistical analysis that evaluates for significance all the interactions that can potentially occur in your dataset: i.e. between ALL the potential celltype pairs. Here, CellphoneDB uses empirical shuffling to calculate which ligand–receptor pairs display significant cell-type specificity. Specifically, it estimates a null distribution of the mean of the average ligand and receptor expression in the interacting clusters by randomly permuting the cluster labels of all cells. The P value for the likelihood of cell-type specificity of a given receptor–ligand complex is calculated on the basis of the proportion of the means that are as high as or higher than the actual mean. 
    - Example command:
    ```shell
    cellphonedb method statistical_analysis test_meta.txt test_counts.txt
    ```
    -  Output: If the user uses the statistical inference approach (`method statistical_analysis`), additional "pvalues.csv" and "significant_means.csv" file are generated with the values for the significant interactions. Finally, ligand–receptor pairs are ranked on the basis of their total number of significant P values across the cell populations. 

- METHOD 3 **degs_analysis** (>= v3): We recently introduced a novel __method to query the database__ alternative to the statistical inference approach. This approach allows the user to design more complex comparisons to define genes specific to a cell type. This is particularly relevant when your research question goes beyond comparing "one" cell type vs "rest". Examples of alternative contrasts are hierarchical comparisons (e.g. you are interested in a specific lineage, such epithelial cells, and wish to identify the genes changing their expression within this lineage) or comparing disease vs control (e.g. you wish to identify upregulated genes in disease T cells by comparing them against control T cells).  For this CellphoneDB method (`method degs_analysis`), the user provides an input file (`test_DEGs.txt` in the command below) indicating which genes are relevant for a cell type (for example, marker genes or significantly upregulated genes resulting from a differential expression analysis (DEG)). CellphoneDB will select interactions where: 
   - (i) all the genes in the interaction are expressed in the corresponding celltype by more than 10% of cells (`--threshold 0.1`) and 
   - (ii) at least one gene-cell type pair is in the provided DEG.tsv file. 

   The user can identify marker genes or DEGs using their preferred tool (we provide [notebooks](https://github.com/ventolab/CellphoneDB/tree/master/notebooks) for both Seurat and Scanpy users) and input the information to CellphoneDB via a [text file](https://github.com/ventolab/CellphoneDB/blob/master/README.md#preparing-your-degs-file-optional-if-method-degs_analysis). 
   - Example command: 
   ```shell
   cellphonedb method degs_analysis test_meta.txt test_counts.txt test_DEGs.txt --threshold 0.1
   ```
   -  Output: If the user uses the `degs_analysis` approach, "relevant_interactions.txt" (instead of  "pvalues.csv") and "significant_means.csv" files are generated.  

For large datasets, do not use .txt files for counts `test_counts.txt`. Input counts as h5ad (recommended), h5 or a path to a folder containing a 10x output with mtx/barcode/features files. NOTE that your gene/protein ids must be HUMAN. If you are working with another specie such as mouse, we recommend you to convert the gene ids to their corresponding orthologous.


## METHOD 2 Statistical inference of receptor-ligand specificity

With this `statistical_analysis` method, we predict enriched receptor–ligand interactions between two cell types based on expression of a receptor by one cell type and a ligand by another cell type, using scRNA-seq data. To identify the most relevant interactions between cell types, we look for the cell-type specific interactions between ligands and receptors. Only receptors and ligands expressed in more than a user-specified threshold percentage of the cells in the specific cluster are considered significant (default is 0.1).

We then perform pairwise comparisons between all cell types. First, we randomly permute the cluster labels of all cells (1,000 times as a default) and determine the mean of the average receptor expression level in a cluster and the average ligand expression level in the interacting cluster. For each receptor–ligand pair in each pairwise comparison between two cell types, this generates a null distribution. By calculating the proportion of the means which are as or higher than the actual mean, we obtain a p-value for the likelihood of cell-type specificity of a given receptor–ligand complex. We then prioritize interactions that are highly enriched between cell types based on the number of significant pairs, so that the user can manually select biologically relevant ones. For the multi-subunit heteromeric complexes, we require that all subunits of the complex are expressed (using a user-specified threshold), and therefore we use the member of the complex with the minimum average expression to perform the random shuffling.

### Optional: Cell subsampling for accelerating analyses
Technological developments and protocol improvements have enabled an exponential growth of the number of cells obtained from scRNA-seq experiments27. Large-scale datasets can profile hundreds of thousands cells, which presents a challenge for the existing analysis methods in terms of both memory usage and runtime. In order to improve the speed and efficiency of our protocol and facilitate its broad accessibility, we integrated subsampling as described in Hie et al.28. This "geometric sketching" approach aims to maintain the transcriptomic heterogeneity within a dataset with a smaller subset of cells. The subsampling step is optional, enabling users to perform the analysis either on all cells, or with other subsampling methods of their choice.




# Interpreting the outputs


## How to read and interpret the results?

The key files are `significant_means.txt` (for statistical_analysis) or `relevant_interactions.txt` (for degs_analysis), see below. When interpreting the results, we recommend you **first define your questions of interest**. Next, focus on specific cell type pairs and manually review the interactions prioritising those with lower p-value and/or higher mean expression. Then, select the cell type pairs and proteins of interest to generate the dotplots for a visual representation. See the options `--columns` and `--rows` to tweak the `dotplot` [here](https://github.com/ventolab/CellphoneDB/blob/master/README.md#dot_plot).


CellphoneDB output is high-throughput. CellphoneDB provides all cell-cell interactions that may potentially occur in your dataset, given the expression of the cells. The size of the output may be overwhelming, but if you apply some rationale (which will depend on the design of your experiment and your biological question), you will be able to narrow it down to a few candidate interactions. The new method `degs_analysis` will allow you perform a more tailored analysis towards specific cell-types or conditions, while the option `---microenvs` will allow you to restrict the combinations of cell-type pairs to test.


It may be that not all of the cell-types of your input dataset co-appear in time and space. Cell types that do not co-appear in time and space will not interact. For example, you might have cells coming from different in vitro systems, different developmental stages or disease and control conditions. Use this prior information to restrict and ignore unfeasible cell-type combinations from the outputs (i.e., columns) as well as their associated interactions (i.e. rows). You can restrict the analysis to feasible cell-type combinations using the option `---microenvs`. Here you can input a two columns file indicating which cell type is in which spatiotemporal microenvironment (see [example](https://github.com/ventolab/CellphoneDB/blob/master/README.md#preparing-your-microenviroments-file-optional-if---microenvs) ). CellphoneDB will use this information to define possible pairs of interacting cells (i.e. pairs of clusters co-existing in a microenvironment) ignoring the rest of combinations. 


## Why values of clusterA-clusterB are different to the values of clusterB-clusterA?

When __reading the outputs__, is IMPORTANT to note that the interactions are not symmetric. Partner A expression is considered for the first cluster/cell type (clusterA), and partner B expression is considered on the second cluster/cell type (clusterB). Thus, `IL12`-`IL12 receptor` for clusterA-clusterB (i.e. the receptor is in clusterB) is not the same that `IL12`-`IL12 receptor` for clusterB-clusterA (i.e. the receptor is in clusterA), and will have different values.

In other words:
* clusterA_clusterB = clusterA expressing partner A and clusterB expressing partner B.
* clusterA_clusterB and clusterB_clusterA  values will be different.



![cellphoneDB methods](https://github.com/ventolab/CellphoneDB/blob/master/Docs/interpreting_results.png)


# Output files

All files (except "deconvoluted.txt") follow the same structure: rows depict interacting proteins while columns represent interacting cell type pairs. 

- The "means.txt" file contains mean values for each ligand-receptor interaction (rows) for each cell-cell interacting pair (columns). 
- The "pvalues.txt" contains the P values for the likelihood of cell-type specificity of a given receptor–ligand complex (rows) in each cell-cell interacting pair (columns), resulting from the `statistical_analysis`. 
- The "significant_means.txt" contains the mean expression (same as "means.txt") of the significant receptor–ligand complex, only. This is the result of crossing "means.csv" and "pvalues.txt".
- The "relevant_interactions.txt" contains a binary matrix indicating if the interaction is relevant (1) or not (0). An interaction is relevant if a gene is a DEG in a cluster/cell type (information from the user provided in the DEG.tsv file) and all the participant genes are expressed, the interaction will be classified as relevant. Alternatively, the value is set to 0. This file is specific to `degs_analysis`. Each row corresponds to a ligand-receptor interaction, while each column corresponds to a cell-cell interacting pair.
- The "deconvoluted.txt" file gives additional information for each of the interacting partners. This is important as some of the interacting partners are heteromers. In other words, multiple molecules have to be expressed in the same cluster in order for the interacting partner to be functional. 


See below the meaning of each column in the outputs:

### P-value (pvalues.txt), Mean (means.txt), Significant mean (significant_means.txt) and Relevant interactions (relevant_interactions.txt)

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

Again, remember that the interactions are not symmetric. It is not the same `IL12`-`IL12 receptor` for clusterA clusterB (i.e. receptor is in clusterB) that `IL12`-`IL12 receptor` for clusterB clusterA (i.e. receptor is in clusterA).


### Deconvoluted (deconvoluted.txt)

* gene_name: Gene identifier for one of the subunits that are participating in the interaction defined in the “means.csv” file. The identifier will depend on the input of the user list.
* uniprot: UniProt identifier for one of the subunits that are participating in the interaction defined in the “means.csv” file.
* is_complex: True if the subunit is part of a complex. Single if it is not, complex if it is.
* protein_name: Protein name for one of the subunits that are participating in the interaction defined in the “means.csv” file.
* complex_name: Complex name if the subunit is part of a complex. Empty if not.
* id_cp_interaction: Unique CellPhoneDB identifier for each of the interactions stored in the database.
* mean: Mean expression of the corresponding gene in each cluster.


# Database design and generation
## Database input files
CellPhoneDB stores ligand-receptor interactions as well as other properties of the interacting partners, including their subunit architecture and gene and protein identifiers. In order to create the content of the database, four main .csv data files are required: "genes_input.csv", "proteins_input.csv", " complexes_input.csv" and "interactions_input.csv" (Figure 4).

### 1. "gene_input"
Mandatory fields: "gene_name"; "uniprot"; "hgnc_symbol" and "ensembl"
This file is crucial for establishing the link between the scRNAseq data and the interaction pairs stored at the protein level. It includes the following gene and protein identifiers: i) gene name ("gene_name"); ii) UniProt identifier ("uniprot"); iiii) HUGO nomenclature committee symbol (HGNC) ("hgnc_symbol") and iv) gene ensembl identifier (ENSG) ("ensembl"). In order to create this file, lists of linked proteins and genes identifiers are downloaded from UniProt and merged using gene names. Several rules need to be considered when merging the files:

UniProt annotation prevails over the gene Ensembl annotation when the same gene Ensembl points towards different UniProt identifiers.
UniProt and Ensembl lists are also merged by their UniProt identifier but this information is only used when the UniProt or Ensembl identifier are missing in the original merged list by gene name.
If the same gene name points towards different HGNC symbols, only the HGNC symbol matching the gene name annotation is considered.
Only one HLA isoform is considered in our interaction analysis and it is stored in a manually HLA-curated list of genes, named "HLA_curated".

### 2. "protein_input"
Mandatory fields: "uniprot"; "protein_name"
Optional fields: "transmembrane"; "peripheral"; "secreted"; "secreted_desc"; "secreted_highlight"; "receptor"; "receptor_desc" ; "integrin"; "other"; "other_desc"; "tags"; "tags_ description"; "tags_reason"
Two types of input are needed to create this file: i) systematic input using UniProt annotation, and ii) manual input using curated annotation both from developers of CellPhoneDB ("proteins_curated") and users. For the systematic input, the UniProt identifier ("uniprot") and the name of the protein ("protein_name") are downloaded from UniProt. For the curated input, developers and users can introduce additional fields relevant to the future systematic assignment of ligand-receptor interactions (see below the "Systematic input from other databases" section for interaction_list). Importantly, the curated information always has priority over the systematic information. The optional input is organised using the fields described below:

#### a. Location of the protein in the cell
There are four non-exclusive options: transmembrane ("transmembrane"), peripheral ("peripheral") and secreted ("secreted", "secreted_desc" and "secreted_highlight").

Plasma membrane proteins are downloaded from UniProt using the keyword KW-1003 (cell membrane). Peripheral proteins from the plasma membrane are annotated using the UniProt keyword SL-9903, and the remaining proteins are annotated as transmembrane proteins. We complete our lists of plasma transmembrane proteins by doing an extensive manual curation using literature mining and UniProt description of proteins with transmembrane and immunoglobulin-like domains.

Secreted proteins are downloaded from UniProt using the keyword KW-0964 (secreted), and are further annotated as cytokines (KW-0202), hormones (KW-0372), growth factors (KW-0339) and immune-related using UniProt keywords and manual annotation. Cytokines, hormones, growth factors and other immune-related proteins are indicated as "secreted_highlight" in the protein_input lists. "secreted_desc" indicates a quality for the secreted protein.

All the manually annotated information is carefully tagged and can be identified. Please see the "curation tags" section below.

#### b. Receptors and integrins
Three fields are allocated to annotate receptors or integrins: "receptor", "receptor_desc" and "integrin".

Receptors are defined by the UniProt keyword KW-0675. The receptors list is extensively reviewed and new receptors are added based on UniProt description and bibliography revision. Receptors involved in immune-cell communication are carefully curated. For some of the receptors, a short description is included in "receptor_desc".

Integrin is a manual curation field that indicates the protein is part of the integrin family. All the annotated information is carefully tagged and can be identified. For details, see "curation tags" section below.

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

These lists contain the UniProt identifiers for each of the heteromeric ligands and receptors ("uniprot1", "uniprot2",...) and a name given to the complex ("complex_name"). Common fields with "proteins_input" include: "transmembrane", "peripheral", "secreted", "secreted_desc", "secreted_highlight", "receptor", "receptor_desc", "integrin", "other", "other_desc" (see description in the above "protein_input" section for clarification). In addition, we include other optional information that may be relevant for the stoichiometry of the heterodimers. If heteromers are defined in the RCSB Protein Data Bank (http://www.rcsb.org), structural information is included in our CellPhoneDB annotation in the "pdb_structure", "pdb_id" and "stoichiometry". An additional field "comments_complex" was created to add a brief description of the heteromer.

### 4. "interaction_input"
Mandatory fields: "partner_a"; "partner_b"; "annotation_strategy"; "source"
Optional fields: "protein_name_a"; "protein_name_b"
Interactions stored in CellPhoneDB are annotated using their UniProt identifier (binary interactions) or the name of the complex (interactions involving heteromers) ("partner_a" and "partner_b"). The name of the protein is also included, yet not mandatory ("protein_name_a" and "protein_name_b"). Protein names are not stored in the database.

There are two main inputs of interactions: i) a systematic input querying other databases, and ii) a manual input using curated information from CellPhoneDB developers ("interactions_curated") and users. The method used to assign the interaction is indicated in the "annotation_strategy" column.

Each interaction stored has a CellPhoneDB unique identifier ("id_cp_interaction") generated automatically by the internal pipeline.

#### a. Systematic input from other databases
We consider interacting partners as: (i) binary interactions annotated by IUPHAR (http://www.guidetopharmacology.org) and (ii) cytokines, hormones and growth factors interactions annotated by innateDB (https://www.innatedb.com) and the iMEX consortium (https://www.imexconsortium.org).

Binary interactions from IUPHAR are directly downloaded from "http://www.guidetopharmacology.org/DATA/interactions.csv" and "guidetopharmachology.org"  is indicated in the "annotation_strategy" field. For the iMEX consortium all protein-protein interactions are downloaded using the PSICQUIC REST APIs19. The IMEx20, IntAct21, InnateDB22, UCL-BHF, MatrixDB23, MINT24, I2D25, UniProt, MBInfo registries are used. Interacting partners are defined as follows: 

- Interacting partner A has to be transmembrane, receptor and cannot be classified as "others" (see the previous "protein_list" section for more information)
- Interacting partner B has to be "secreted_highlight". This group of proteins includes cytokines, hormones, growth factors and other immune-related proteins (see the previous "protein_list" section for more information).
In both cases, some interactions are excluded: i) interactions where one of the components is part of a complex (see "complex_curated" list in the above section); ii) interactions which are not involved in cell-cell communication or are wrongly annotated by our systematic method. These are stored in a curated list of proteins named "excluded_interaction". The "excluded_interaction" file contains five fields: a) uniprot_1: name of the interacting partner A that is going to be excluded; b) uniprot_2: name of the interacting partner B that is going to be excluded; c) name.1: name of the protein to be excluded corresponding to uniprot_1; d) name.2: name of the protein to be excluded corresponding to uniprot_2; e) comments: information about the exclusion of the protein.

Homomeric complexes - proteins interacting with themselves - are excluded from the systematic analysis. Importantly, in cases where both the systematic and the curated input detect the interactions, the curated input always prevails over the systematic information.

#### b. Curated approach
The majority of ligand–receptor interactions are manually curated by reviewing UniProt descriptions and PubMed information on membrane receptors. Cytokine and chemokine interactions were annotated following the International Union of Pharmacology annotation 26. The interactions of other groups of cell-surface proteins were manually reviewed, including the TGF family, integrins, lymphocyte receptors, semaphorins, ephrins, Notch and TNF receptors. The bibliography used to annotate the interaction is stored in "source". ‘Uniprot’ indicates that the interaction has been annotated using UniProt descriptions.

## User-defined receptor-ligand datasets
Our system allows users to create their own lists of curated proteins and complexes. In order to do so, the format of the users’ lists must be compatible with the input files. Users can submit their lists using the Python package version of CellPhoneDB, and then send them via email, the cellphonedb.org form, or a pull request to the CellPhoneDB data repository (https://github.com/Teichlab/cellphonedb-data).

### Database structure
Information is stored in an SQLite relational database (https://www.sqlite.org). SQLAlchemy (www.sqlalchemy.org) and Python 3 were used to build the database structure and the query logic. The application is designed to allow analysis on potentially large count matrices to be performed in parallel. This requires an efficient database design, including optimisation for query times, indices and related strategies. All application code is open source and uploaded both to github and the web server.  An explanation of the content of the tables and the database schema is available in Supplementary methods and Supplementary Figure 1 and Figure 2.



