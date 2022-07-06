# CellPhoneDB

## What is CellPhoneDB?

CellPhoneDB is a publicly available repository of curated receptors, ligands and their interactions in **HUMAN**. Subunit architecture is included for both ligands and receptors, representing heteromeric complexes accurately. This is crucial, as cell-cell communication relies on multi-subunit protein complexes that go beyond the binary representation used in most databases and studies. CellPhoneDB integrates [existing datasets](Docs/ppi-resources.md) that pertain to cellular communication and new manually reviewed information. Databases from which CellPhoneDB gets information are: UniProt, Ensembl, PDB, the IMEx consortium, IUPHAR.

CellPhoneDB can be used to search for a particular ligand/receptor, or interrogate your own single-cell transcriptomics data (or even bulk transcriptomics data if your samples represent pure populations!). 

For more details on the analysis check the [documentation here](Docs/RESULTS-DOCUMENTATION.md), our protocols paper [Efremova et al 2020](https://www.nature.com/articles/s41596-020-0292-x) or [Garcia-Alonso et al](https://www.nature.com/articles/s41588-021-00972-2) (for CellphoneDB v3).

## New in CellPhoneDB v4


This release involves a major **Database Update**. We have invested quite some time curating more cell-cell communication interactions validated experimentally. Specifically, we have:

 1. Manually curated more protein-protein interactions involved in cell-cell communication, with special focus on protein acting as heteromeric complexes. The new database includes almost **2,000 high-confidence interactions**! We believe modelling complexes is key to minimise false positive predictions.
 2. Annotated non-peptidic molecules (i.e., not encoded by a gene) acting as ligands. Examples of these include steroid hormones (e.g., estrogen). To do so, we have reconstructed the biosynthetic pathways and use the last representative enzyme as a proxy of ligand abundance. We retrieve this information by manually reviewing and curating relevant literature and peer-reviewed pathway resources such as REACTOME. We include more than **200 interactions involving non-peptidic ligands**!


Check [Garcia-Alonso & Lorenzi et al](https://www.nature.com/articles/s41586-022-04918-4) for an example applying CellphoneDB v4.

To use the lastest version, re-install cellphonedb & download the database:

```
pip install -U cellphonedb
cellphonedb database download
```



## New in CellPhoneDB v3

1. **Incorporate spatial information** CellPhoneDB now allows the incorporation of spatial information of the cells via the `microenvironments` file. This is a two columns file indicating which cell type is in which spatial microenvironment (see [example](https://github.com/ventolab/CellphoneDB/blob/master/in/endometrium_atlas_example/endometrium_example_microenviroments.tsv) ). CellphoneDB will use this information to define possible pairs of interacting cells (i.e. pairs of clusters sharing/coexisting in a microenvironment). You can define microenvironments with prior knowledge, imaging or Visium analysis with [cell2location](https://cell2location.readthedocs.io/en/latest/notebooks/cell2location_short_demo_downstream.html#4.-Identify-groups-of-co-located-cell-types-using-matrix-factorisation).
2. **New analysis method added**,  using differentially expressed genes (DEGs) instead of random shuffling (`cellphonedb method degs_analysis`). This approach will select interactions where all the genes are expressed by a fraction of cells above a `--threshold` and at least one gene is a DEG. The user can identify DEGs using their preferred tool and provide the information to CellphoneDB via text file. The first column should be the cell type/cluster and the second column the associated gene id. The remaining columns are ignored (see [example](https://github.com/ventolab/CellphoneDB/blob/master/in/endometrium_atlas_example/endometrium_example_DEGs.tsv) ). We provide [notebooks](https://github.com/ventolab/CellphoneDB/tree/master/notebooks) for both Seurat and Scanpy users.  
3. **Database update** WNT pathway has been further curated.

Check [Garcia-Alonso et al](https://www.nature.com/articles/s41588-021-00972-2) for an example applying CellphoneDB v3.

## Installing CellPhoneDB
NOTE: Works with Python v3.6 or greater. If your default Python interpreter is for `v2.x` (you can check it with `python --version`), calls to `python`/`pip` should be substituted by `python3`/`pip3`.

We highly recommend using an isolated python environment (as described in steps 1 and 2) using [conda](https://docs.conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html) or [virtualenv](https://docs.python.org/3/library/venv.html) but you could of course omit these steps and install via `pip` immediately.

1. Create python=>3.6 environment
- Using conda: `conda create -n cpdb python=3.7`
- Using virtualenv: `python -m venv cpdb`

2. Activate environment
- Using conda: `source activate cpdb`
- Using virtualenv: `source cpdb/bin/activate`

3. Install CellPhoneDB `pip install cellphonedb`


## Running CellPhoneDB Methods

Please, activate your environment if you didn't previously
- Using conda: `source activate cpdb`
- Using virtualenv: `source cpdb/bin/activate`

To use the example data, please [download meta/counts test data](https://github.com/Teichlab/cellphonedb/blob/master/in/example_data/cellphonedb_example_data.zip?raw=true). i.e.
```shell
curl https://raw.githubusercontent.com/Teichlab/cellphonedb/master/in/example_data/test_counts.txt --output test_counts.txt
curl https://raw.githubusercontent.com/Teichlab/cellphonedb/master/in/example_data/test_meta.txt --output test_meta.txt
```

### Prepatring INPUTS
#### Preparing your counts input file (mandatory)
Counts file can be a text file or a `h5ad` (recommended), `h5` or a path to a folder containing a 10x output with `mtx/barcode/features` files. NOTE: Your gene/protein **ids must be HUMAN**. If you are working with another specie such as mouse, we recommend you to convert the gene ids to their corresponding orthologous. 

#### Preparing your DEGs file (optional, if `method degs_analysis`)
This is a two columns file indicanting which gene is specific or upregulated in a cell type (see [example](https://github.com/ventolab/CellphoneDB/blob/master/in/endometrium_atlas_example/endometrium_example_DEGs.tsv) ). The first column should be the cell type/cluster name (matching those in `meta.txt`) and the second column the associated gene id. The remaining columns are ignored. We provide [notebooks](https://github.com/ventolab/CellphoneDB/tree/master/notebooks) for both Seurat and Scanpy users. It is on you to design a DEG analysis appropiated for your research question. 

#### Preparing your microenviroments file (optional, if `--microenvs`)
This is a two columns file indicating which cell type is in which spatial microenvironment (see [example](https://github.com/ventolab/CellphoneDB/blob/master/in/endometrium_atlas_example/endometrium_example_microenviroments.tsv) ). CellphoneDB will use this information to define possible pairs of interacting cells (i.e. pairs of clusters co-appearing in a microenvironment). 

### RUN examples
####  Example with running the DEG-based method
```shell
cellphonedb method degs_analysis test_meta.txt test_counts.txt test_DEGs.txt
```

####  Example with running the statistical method
```shell
cellphonedb method statistical_analysis test_meta.txt test_counts.txt
```

#### Example without using the statistical method
 - **Using text files**
```shell
cellphonedb method analysis test_meta.txt test_counts.txt 
```

 - **Using h5ad count file**
```shell
cellphonedb method analysis test_meta.txt test_counts.h5ad
```

####  Example running a microenviroments file
```shell
cellphonedb method statistical_analysis test_meta.txt test_counts.txt --microenvs test_microenvs.txt
```

####  Example running the DEG-based method with microenviroments file
```shell
cellphonedb method degs_analysis test_meta.txt test_counts.txt test_DEGs.txt --microenvs test_microenvs.txt
```

To understand the different analysis and results, please check the [results documentation](Docs/RESULTS-DOCUMENTATION.md).


### Optional Parameters

~ **Optional Method parameters**:
- `--counts-data`: [ensembl \| gene_name \| hgnc_symbol] Type of gene identifiers in the counts data
- `--project-name`: Name of the project. A subfolder with this name is created in the output folder
- `--iterations`: Number of iterations for the statistical analysis [1000]
- `--threshold`: % of cells expressing the specific ligand/receptor
- `--result-precision`: Number of decimal digits in results [3]
- `--output-path`: Directory where the results will be allocated (the directory must exist) [out]
- `--output-format`: Output format of the results files (extension will be added to filename if not present) [txt]
- `--means-result-name`: Means result filename [means]
- `--significant-means-result-name`: Significant mean result filename [significant_means]
- `--deconvoluted-result-name`: Deconvoluted result filename [deconvoluted]
- `--verbose/--quiet`: Print or hide CellPhoneDB logs [verbose]
- `--subsampling`: Enable subsampling
- `--subsampling-log`: Enable subsampling log1p for non log-transformed data inputs !!mandatory!!
- `--subsampling-num-pc`: Subsampling NumPC argument (number of PCs to use) [100]
- `--subsampling-num-cells`: Number of cells to subsample to [1/3 of cells]


~ **Optional Method Statistical parameters**
- `--microenvs`: Spatial microenviroments input file. Restricts the cluster/cell_type interacting pairs to the cluster/cell_type sharing a microenviroment (i.e. only test a combination of clusters if these coexist in a microenviroment). This file should contain two columns: 1st column indicates the cluster/cell_type, 2nd column indicates the microenviroment name.  See example [here](https://github.com/ventolab/CellphoneDB/tree/master/in). 
- `--pvalues-result-name`: P-values result filename [pvalues]
- `--pvalue`: P-value threshold [0.05]
- `--debug-seed`: Debug random seed -1. To disable it please use a value >=0 [-1]
- `--threads`: Number of threads to use. >=1 [4]

### Usage Examples

Set number of iterations and threads
```shell
cellphonedb method statistical_analysis yourmetafile.txt yourcountsfile.txt --iterations=10 --threads=2
```

Set project subfolder
```shell
cellphonedb method analysis yourmetafile.txt yourcountsfile.txt --project-name=new_project
```

Set output path
```shell
mkdir custom_folder
cellphonedb method statistical_analysis yourmetafile.txt yourcountsfile.txt --output-path=custom_folder
```

Subsampling
```shell
cellphonedb method analysis yourmetafile.txt yourcountsfile.txt --subsampling --subsampling-log false --subsampling-num-cells 3000
```

## Plotting statistical method results

In order to plot results from the statistical methods, you need to run it first.

Currently there are two plot types available: `dot_plot` & `heatmap_plot`

Once you have the needed files (`means` & `pvalues` from method statistical_analysis or `means` & `relevant_interactions` from method degs_analysis) you can proceed as follows:
```shell
cellphonedb plot dot_plot
```

```shell
cellphonedb plot heatmap_plot yourmeta.txt
```

### `dot_plot`
This plot type requires `ggplot2` R package installed and working

You can tweak the options for the plot with these arguments:
- `--means-path`: The means output file [./out/means.txt]
- `--pvalues-path`: The pvalues output file [./out/pvalues.txt]
- `--output-path`: Output folder [./out]
- `--output-name`: Filename of the output plot [plot.pdf]
- `--rows`: File with a list of rows to plot, one per line [all available]
- `--columns`: File with a list of columns to plot, one per line [all available]
- `--verbose / --quiet`: Print or hide CellPhoneDB logs [verbose]

Available output formats are those supported by `R's ggplot2` package, among others they are:
- `pdf`
- `png`
- `jpeg`

This format will be inferred from the `--output-name` argument

To plot only desired rows/columns (samples for [rows](in/example_data/rows.txt) and [columns](in/example_data/columns.txt) based in example data files):
```shell
cellphonedb plot dot_plot --rows in/rows.txt --columns in/columns.txt
``` 

### `heatmap_plot`
This plot type requires `pheatmap` R package installed and working 
This plot type includes two features `count` & `log_count` 

You can tweak the options for the plot with these arguments:
- `--pvalues-path`: The pvalues output file [./out/pvalues.txt]
- `--output-path`: Output folder [./out]
- `--count-name`: Filename of the output plot [heatmap_count.pdf]
- `--log-name`: Filename of the output plot using log-count of interactions [heatmap_log_count.pdf]
- `--count-network-name`: Filename of the output network file [count_network.txt]
- `--interaction-count-name`: Filename of the output interactions-count file [interactions_count.txt] 
- `--pvalue`: pvalue threshold to consider when plotting [0.05]
- `--verbose / --quiet`: Print or hide cellphonedb logs [verbose]

Available output formats are those supported by `R's pheatmap` package, among others they are:
- `pdf`
- `png`
- `jpeg`

This format will be inferred from the `--count-name` & `--log-name` arguments.


## Using different database versions
CellPhoneDB databases can be updated from the remote repository through our tool. Furthermore, available versions can be listed and downloaded for use. 

To use one of those versions, a user must provide the argument `--database <version_or_file>` to the method to be executed.

If the given parameter is a readable database file, it will be used as is. Otherwise it will use some of the versions matching the selected version.

If the selected version does not exist in the local environment it will be downloaded from the remote repository. (See below.) If no `--database` argument is given in methods execution, it will use the latest local version available.

Downloaded versions will be stored in a user folder under `~/.cpdb/releases`

## Listing remote available versions
The command to list available versions from the remote repository is:
```shell
cellphonedb database list_remote
``` 

## Listing local available versions
The command to list available versions from the local repository is:
```shell
cellphonedb database list_local
``` 

## Download version
The command to download a version from the remote repository is:
```shell
cellphonedb database download
``` 
or 

```shell 
cellphonedb database download --version <version_spec|latest> 
``` 

whereby `version_spec` must be one of the listed in the `database list_remote` command.
If no version is specified or `latest` is used as a `version_spec`, the newest available version will be downloaded


## Generating user-specific custom database
A user can generate custom databases and use them. In order to generate a new database, a user can provide his/her own lists.

These lists can be: genes, proteins, complexes and/or interactions. In the generation process they will get merged with the ones from the CellPhoneDB release sources. The user lists have higher precedence than the ones included in CellPhoneDB package.

To generate such a database the user has to issue this command:
```shell
cellphonedb database generate  
```
Generate specific parameters:

- `--user-protein`: Protein input file
- `--user-gene`: Gene input file
- `--user-complex`: Complex input file
- `--user-interactions`: Interactions input file
- `--fetch`: Some lists can be downloaded from original sources while creating the database, eg: uniprot, ensembl. By default, the snapshots included in the CellPhoneDB package will be used; to enable a fresh copy `--fetch` must be appended to the command
- `--result-path`: Output folder
- `--log-file`: Log file
- `--user-interactions-only`: Use only interactions provided.

Result database file is generated in the folder `out` with `cellphonedb_user_{datetime}.db`. The user defined input tables will be merged with the current CellPhoneDB input tables. To use this database, please use `--database` parameter in methods.
E.g:
```
 cellphonedb method statistical_analysis in/example_data/test_meta.txt in/example_data/test_counts.txt --database out/cellphonedb_user_2019-05-10-11_10.db
```

### Examples for user-specific custom database

1. To add or correct some interactions

    **Input**:
    
    - `your_custom_interaction_file.csv`: Comma separated file (use mandatory columns!) with interactions to add/correct.
     
    **Command**:
    
    ```shell
    cellphonedb database generate --user-interactions your_custom_interaction_file.csv 
    ```
    **Result**:
    
    New database file with CellPhoneDB interactions + user custom interactions. For duplicated interactions, user lists overwrite the CellPhoneDB original data.

2. To use only user-specific interactions
    
    **Input**:
    
    - `your_custom_interaction_file.csv`: Comma separated file (use mandatory columns!) with interactions to use.
    
    **Command**:
    
    ```shell
    cellphonedb database generate --user-interactions your_custom_interaction_file.csv --user-interactions-only
    ```
    **Result**:
    
    New database file with **only** user custom interactions. 

3. To correct any protein data
    
    **Input**:
    - `your_custom_protein_file.csv`: Comma separated file (use mandatory columns!) with proteins to overwrite. 
     
    **Command**:
    ```shell
    cellphonedb database generate --user-protein your_custom_protein_file.csv 
    ```
    **Result**:
    
    New database file with CellPhoneDB interactions + user custom interactions. For duplicated interactions or proteins, user lists overwrite the CellPhoneDB original data.

4. To add some interactions and correct any protein data
    
    **Input**:
    
    - `your_custom_interaction_file.csv`: Comma separated file (use mandatory columns!) with interactions to add/correct.
    - `your_custom_protein_file.csv`: Comma separated file (use mandatory columns!) with proteins to overwrite. 
    
    **Command**:
    
    ```shell
    cellphonedb database generate --user-interactions your_custom_interaction_file.csv --user-protein your_custom_protein_file.csv 
    ```
    **Result**:
    
    New database file with CellPhoneDB interactions + user custom interactions. On duplicated interactions or proteins, user list overwrites CellPhoneDB original data.    
     
5. To update remote sources (UniProt, IMEx, ensembl, etc.)     
    
    **IMPORTANT**
    
    This command uses external resources allocated in external servers. The command may not end correctly if external servers are not available. The timing of this step depends on external servers and the user's internet connection and can take a lot of time. 
    
    **Input**:
    
    - `your_custom_interaction_file.csv`: Comma separated file (use mandatory columns!) with interactions to add/correct.
    - `your_custom_protein_file.csv`: Comma separated file (use mandatory columns!) with proteins to overwrite. 
    
    **Command**:
    
    ```shell
    cellphonedb database generate --fetch 
    ```
    **Result**:
    
    New database file with CellPhoneDB interactions + user custom interactions. For duplicated interactions or proteins, user lists overwrite the CellPhoneDB original data.
    
    
Some lists can be downloaded from original sources while creating the database, e.g. uniprot, ensembl. By default, the snapshots included in the CellPhoneDB package will be used---to enable a fresh copy `--fetch` must be appended to the command.

In order to use specific lists those can be specified like this `--user-protein`, `--user-gene`, `--user-complex`, `--user-interactions`, `--user-interactions-only`
followed by the corresponding file path.

The database file can be then used as explained below. The intermediate lists used for the generation will be saved along the database itself.

As the lists are processed, then filtered, and lastly collected, two versions may exist: `_generated` is the unfiltered one whereas `_input` is the final state prior being inserted in the database.


## Contributing to CellPhoneDB

CellPhoneDB is an open-source project. If you are interested in contributing to this project, please let us know.

You can check all project documentation in the [Docs](Docs) section
