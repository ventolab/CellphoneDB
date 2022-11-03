# CellPhoneDB 

[![Database](https://img.shields.io/github/v/release/ventolab/cellphonedb-data.svg?color=blue&label=database)](https://github.com/ventolab/CellphoneDB-data)
[![Python package](https://img.shields.io/pypi/v/cellphonedb.svg?color=brightgreen&label=python-package)](https://pypi.org/project/cellphonedb)

## What is CellPhoneDB?

CellPhoneDB is a publicly available repository of curated receptors, ligands and their interactions in **HUMAN**. Subunit architecture is included for both ligands and receptors, representing heteromeric complexes accurately. This is crucial, as cell-cell communication relies on multi-subunit protein complexes that go beyond the binary representation used in most databases and studies. <font color="red">CellPhoneDB integrates [existing datasets](Docs/ppi-resources.md) that pertain to cellular communication and new manually reviewed information. Databases from which CellPhoneDB gets information are: UniProt, Ensembl, PDB, the IMEx consortium, IUPHAR.</font>

CellPhoneDB can be used to search for a particular ligand/receptor, or interrogate your own single-cell transcriptomics data (or even bulk transcriptomics data if your samples represent pure populations!). 

<font color="red">For more details on the analysis check the [documentation here](Docs/RESULTS-DOCUMENTATION.md), our protocols paper [Efremova et al 2020](https://www.nature.com/articles/s41596-020-0292-x) or [Garcia-Alonso et al](https://www.nature.com/articles/s41588-021-00972-2) (for CellphoneDB v3).</font>

## New in CellPhoneDB v4

This release involves a major **Database Update**. We have invested quite some time curating more cell-cell communication interactions validated experimentally. Specifically, we have:

 1. Manually curated more protein-protein interactions involved in cell-cell communication, with special focus on protein acting as heteromeric complexes. The new database includes almost **2,000 high-confidence interactions**, including heteromeric complexes! We believe modelling complexes is key to minimise false positives in the predictions.
 2. Annotated non-peptidic molecules (i.e., not encoded by a gene) acting as ligands. Examples of these include steroid hormones (e.g., estrogen). To do so, we have reconstructed the biosynthetic pathways and used the last representative enzyme as a proxy of ligand abundance. We retrieve this information by manually reviewing and curating relevant literature and peer-reviewed pathway resources such as REACTOME. We include more than **200 interactions involving non-peptidic ligands**!


Check [Garcia-Alonso & Lorenzi et al](https://www.nature.com/articles/s41586-022-04918-4) for an example applying CellphoneDB v4.

## New in CellPhoneDB v3

1. **Incorporate spatial information** CellPhoneDB now allows the incorporation of spatial information of the cells via the `microenvironments` file. This is a two columns file indicating which cell type is in which spatial microenvironment (see [example](https://github.com/ventolab/CellphoneDB/blob/master/in/endometrium_atlas_example/endometrium_example_microenviroments.tsv) ). CellphoneDB will use this information to define possible pairs of interacting cells (i.e. pairs of clusters sharing/coexisting in a microenvironment). You can define microenvironments with prior knowledge, imaging or Visium analysis with [cell2location](https://cell2location.readthedocs.io/en/latest/notebooks/cell2location_short_demo_downstream.html#4.-Identify-groups-of-co-located-cell-types-using-matrix-factorisation).
2. **New analysis method added**,  using differentially expressed genes (DEGs) instead of random shuffling (`cellphonedb method degs_analysis`). This approach will select interactions where all the genes are expressed by a fraction of cells above a `--threshold` and at least one gene is a DEG. The user can identify DEGs using their preferred tool and provide the information to CellphoneDB via text file. The first column should be the cell type/cluster and the second column the associated gene id. The remaining columns are ignored (see [example](https://github.com/ventolab/CellphoneDB/blob/master/in/endometrium_atlas_example/endometrium_example_DEGs.tsv) ). We provide [notebooks](https://github.com/ventolab/CellphoneDB/tree/master/notebooks) for both Seurat and Scanpy users.  
3. **Database update** WNT pathway has been further curated.

Check [Garcia-Alonso et al](https://www.nature.com/articles/s41588-021-00972-2) for an example applying CellphoneDB v3.

## Installing CellPhoneDB
NOTE: Works with Python v3.9 or greater. If your default Python interpreter is for `v2.x` (you can check it with `python --version`), calls to `python`/`pip` should be substituted by `python3`/`pip3`.

We highly recommend using an isolated python environment (as described in steps below) using [conda](https://docs.conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html) or [virtualenv](https://docs.python.org/3/library/venv.html) but you could of course omit these steps and install via `pip` immediately.

Note that the instructions below refer to running CellphoneDB in development environment
- `cd <your_workspace>`
- `git clone git@github.com:ventolab/CellphoneDB.git`
- `cd CellphoneDB`
- `git checkout bare-essentials`
- `conda create -n cpdb python=3.9`
- `source activate cpdb`
- `pip install -r requirements.txt`

## Running CellPhoneDB Methods
Please, activate your environment if you didn't previously
- Using conda: `source activate cpdb`
- Using virtualenv: `source cpdb/bin/activate`
- `mkdir -p ~/.cpdb/user_files`
- `cp <your_workspace>/CellphoneDB/example_data/* ~/.cpdb/user_files`
- `jupyter notebook &`
- Follow instructions on http://localhost:8888/notebooks/cellphonedb/cellphonedb.ipynb

### Prepatring INPUTS
#### Preparing your counts input file (mandatory)
Counts file can be a text file or a `h5ad` (recommended), `h5` or a path to a folder containing a 10x output with `mtx/barcode/features` files. NOTE: Your gene/protein **ids must be HUMAN**. If you are working with another specie such as mouse, we recommend you to convert the gene ids to their corresponding orthologous. 

#### Preparing your DEGs file (optional, if `method degs_analysis`)
This is a two columns file indicanting which gene is specific or upregulated in a cell type (see [example](https://github.com/ventolab/CellphoneDB/blob/master/in/endometrium_atlas_example/endometrium_example_DEGs.tsv) ). The first column should be the cell type/cluster name (matching those in `meta.txt`) and the second column the associated gene id. The remaining columns are ignored. We provide [notebooks](https://github.com/ventolab/CellphoneDB/tree/master/notebooks) for both Seurat and Scanpy users. It is on you to design a DEG analysis appropiated for your research question. 

#### Preparing your microenviroments file (optional, if `--microenvs`)
This is a two columns file indicating which cell type is in which spatial microenvironment (see [example](https://github.com/ventolab/CellphoneDB/blob/master/in/endometrium_atlas_example/endometrium_example_microenviroments.tsv) ). CellphoneDB will use this information to define possible pairs of interacting cells (i.e. pairs of clusters co-appearing in a microenvironment). 

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

<font color="red">
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
</font>

## Generating user-specific custom database
A user can generate custom databases and use them. In order to generate a new database, a user can provide his/her own lists.

These lists can be: genes, proteins, complexes and/or interactions. In the generation process they will get merged with the ones from the CellPhoneDB release sources. The user lists have higher precedence than the ones included in CellPhoneDB package.

To generate such a database please do the following (taking v5.0.0 as an example)
- cd `~/.cpdb/releases/'
- `mkdir -p v5.0.0_own; cd v5.0.0_own`
- `cp v5.0.0/*_input.csv .`
- Modify *_input.csv as appropriate
- In http://localhost:8888/notebooks/cellphonedb/cellphonedb.ipynb:
    - `cellophonedb_version = "v5.0.0_own"`
    - `use_local_files=True`
    - `user_dir_root = os.path.join(os.path.expanduser('~'),".cpdb")`
    - `controller.create_db(user_dir_root, cellophonedb_version, use_local_files)`
  The above command will create a ~/.cpdb/releases/v5.0.0_own/cellphonedb.zip file that will be used by all the analysis methods in the above notebook.

## Contributing to CellPhoneDB

CellPhoneDB is an open-source project. If you are interested in contributing to this project, please let us know.

You can check all project documentation in the [Docs](Docs) section
