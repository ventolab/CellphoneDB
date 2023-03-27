[![Database](https://img.shields.io/github/v/release/ventolab/cellphonedb-data.svg?color=blue&label=database)](https://github.com/ventolab/CellphoneDB-data) [![Python package](https://img.shields.io/pypi/v/cellphonedb.svg?color=brightgreen&label=python-package)](https://pypi.org/project/cellphonedb)
# CellPhoneDB 

## What is CellPhoneDB?

CellPhoneDB is a publicly available repository of curated receptors, ligands and their interactions in **HUMAN**. CellPhoneDB can be used to search for a particular ligand/receptor, or interrogate your own single-cell transcriptomics data (or even bulk transcriptomics data if your samples represent pure populations!). 

Subunit architecture is included for both ligands and receptors, representing heteromeric complexes accurately. This is crucial, as cell communication relies on multi-subunit protein complexes that go beyond the binary representation used in most databases and studies. CellPhoneDB also incorporates biosynthetic pathways in which we use the last representative enzyme as a proxy of ligand abundance, by doing so, we include interactions involving non-peptidic CellPhoneDB includes only manually curated & reviewd molecular interactions with evidenced role in cellular communication.

For more details on the analysis check the [DOCUMENTATION](https://cellphonedb.readthedocs.io/en/latest/#). 

Please cite our papers [Vento-Tormo R, Efremova M, et al., 2018](https://www.nature.com/articles/s41596-020-0292-x) (original CellphoneDB) or [Garcia-Alonso et al., 2021](https://www.nature.com/articles/s41586-018-0698-6) (for CellphoneDB method 3).


## New in CellPhoneDB-data v4.1.0

This release of CellphoneDB database integrates new manually reviewed interactions with evidenced roles in cell-cell communication together with existing datasets that pertain to cellular communication (such as Shilts *et al.* 2022 and Kanemura *et al.* 2023). Recently, the database expanded to include non-protein molecules acting as ligands.

1. CellPhoneDB has been implemented as a python package, improving its efficiency and adding new methods, such as the CellPhoneDB results query function.
2. Manually curated more protein-protein interactions involved in cell-cell communication, with a special focus on proteins acting as heteromeric complexes [cellphonedb-data v4.1.0](https://github.com/ventolab/cellphonedb-data). The new database includes more than [2,900 high-confidence interactions](https://www.cellphonedb.org/database.html), including heteromeric complexes. In this version we haved added new G-protein-coupled receptors interactions from Kanemura *et al.* 2023 and  Shilts *et al.* 2022.
3. Interactions retrieved from external resources have been removed from this release to include high-confidence interactions only.
4. [Tutorials](notebooks) for the new CellPhoneDB implementation.

See updates from [previous releases here](https://github.com/ventolab/CellphoneDB/blob/master/docs/RESULTS-DOCUMENTATION.md#release-notes).


## Installing CellPhoneDB 
NOTE: Works with Python v3.8 or greater. If your default Python interpreter is for `v2.x` (you can check it with `python --version`), calls to `python`/`pip` should be substituted by `python3`/`pip3`.

We highly recommend using an isolated python environment (as described in steps 1 and 2) using [conda](https://docs.conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html) or [virtualenv](https://docs.python.org/3/library/venv.html) but you could of course omit these steps and install via `pip` immediately.

1. Create python=>3.8 environment
- Using conda: `conda create -n cpdb python=3.8`
- Using virtualenv: `python -m venv cpdb`

2. Activate environment
- Using conda: `source activate cpdb`
- Using virtualenv: `source cpdb/bin/activate`

3. Install CellPhoneDB `pip install cellphonedb`


## Running CellPhoneDB Methods

Please, activate your environment if you didn't previously
- Using conda: `source activate cpdb`
- Using virtualenv: `source cpdb/bin/activate`

We have created a set of tutorials that can be accessed for each To use the example data, please [tutorials and data](notebooks).

### Prepatring INPUTS
#### Preparing your counts input file (mandatory)
Counts file can be a text file or a `h5ad` (recommended), `h5` or a path to a folder containing a 10x output with `mtx/barcode/features` files. NOTE: Your gene/protein **ids must be HUMAN**. If you are working with another specie such as mouse, we recommend you to convert the gene ids to their corresponding orthologous. 

#### Preparing your DEGs file (optional, if `method degs_analysis`)
This is a two columns file indicanting which gene is specific or upregulated in a cell type (see [example](in/endometrium_atlas_example/endometrium_example_DEGs.tsv) ). The first column should be the cell type/cluster name (matching those in `meta.txt`) and the second column the associated gene id. The remaining columns are ignored. We provide [notebooks](notebooks) for both Seurat and Scanpy users. It is on you to design a DEG analysis appropiated for your research question. 

#### Preparing your microenviroments file (optional, if `microenvs_file_path`)
This is a two columns file indicating which cell type is in which spatial microenvironment (see [example](in/endometrium_atlas_example/endometrium_example_microenviroments.tsv) ). CellphoneDB will use this information to define possible pairs of interacting cells (i.e. pairs of clusters co-appearing in a microenvironment). 

### RUN examples

For more detailed examples refer to out tutorials [here](notebooks).
####  Example with running the DEG-based method
```shell
from cellphonedb.src.core.methods import cpdb_degs_analysis_method

deconvoluted, means, relevant_interactions, significant_means = cpdb_degs_analysis_method.call(
        cpdb_file_path = cellphonedb.zip,
        meta_file_path = test_meta.txt,
        counts_file_path = test_counts.h5ad,
        degs_file_path = degs_file_path,
        counts_data = 'hgnc_symbol',
        threshold = 0.1,
        output_path = out_path)
```

####  Example with running the statistical method
```shell
from cellphonedb.src.core.methods import cpdb_statistical_analysis_method

deconvoluted, means, pvalues, significant_means = cpdb_statistical_analysis_method.call(
        cpdb_file_path = cellphonedb.zip,
        meta_file_path = test_meta.txt,
        counts_file_path = test_counts.h5ad,
        counts_data = 'hgnc_symbol',
        output_path = out_path)
```

#### Example without using the statistical method
 - **Using text files**
```shell
from cellphonedb.src.core.methods import cpdb_analysis_method

means, deconvoluted = cpdb_analysis_method.call(
        cpdb_file_path = cellphonedb.zip,
        meta_file_path = test_meta.txt,
        counts_file_path = test_counts.txt,
        counts_data = 'hgnc_symbol',
        output_path = out_path)
```

 - **Using h5ad count file**
```shell
from cellphonedb.src.core.methods import cpdb_analysis_method

means, deconvoluted = cpdb_analysis_method.call(
        cpdb_file_path = cellphonedb.zip,
        meta_file_path = test_meta.txt,
        counts_file_path = test_counts.h5ad,
        counts_data = 'hgnc_symbol',
        output_path = out_path)
```

####  Example running a microenviroments file
```shell
from cellphonedb.src.core.methods import cpdb_analysis_method

means, deconvoluted = cpdb_analysis_method.call(
        cpdb_file_path = cellphonedb.zip,
        meta_file_path = test_meta.txt,
        counts_file_path = test_counts.h5ad,
        counts_data = 'hgnc_symbol',
        microenvs_file_path = microenvs_file_path,
        output_path = out_path)
```

####  Example running the DEG-based method with microenviroments file
```shell
from cellphonedb.src.core.methods import cpdb_degs_analysis_method

deconvoluted, means, relevant_interactions, significant_means = cpdb_degs_analysis_method.call(
        cpdb_file_path = cellphonedb.zip,
        meta_file_path = test_meta.txt,
        counts_file_path = test_counts.h5ad,
        counts_data = 'hgnc_symbol',
        microenvs_file_path = microenvs_file_path,
        output_path = out_path)
```

To understand the different analysis and results, please check the [results documentation](docs/RESULTS-DOCUMENTATION.md).


### Optional Parameters

~ **Optional Method parameters**:
- `counts_data`: [ensembl \| gene_name \| hgnc_symbol] Type of gene identifiers in the counts data
- `iterations`: Number of iterations for the statistical analysis [1000]
- `threshold`: % of cells expressing the specific ligand/receptor
- `result_precision`: Number of decimal digits in results [3]
- `output_path`: Directory where the results will be allocated (the directory must exist) [out]
- `output_suffix`: Output format of the results files (time stamp will be added to filename if not present) [txt]
- `subsampling`: Enable subsampling
- `subsampling_log`: Enable subsampling log1p for non log-transformed data inputs !!mandatory!!
- `subsampling_num_pc`: Subsampling NumPC argument (number of PCs to use) [100]
- `subsampling_num_cells`: Number of cells to subsample the dataset [1/3 of cells]


~ **Optional Method Statistical parameters**
- `microenvs_file_path`: Spatial microenviroments input file. Restricts the cluster/cell_type interacting pairs to the cluster/cell_type sharing a microenviroment (i.e. only test a combination of clusters if these coexist in a microenviroment). This file should contain two columns: 1st column indicates the cluster/cell_type, 2nd column indicates the microenviroment name.  See example [here](https://github.com/ventolab/CellphoneDB/tree/master/in). 
- `pvalue`: P-value threshold [0.05]
- `debug_seed`: Debug random seed -1. To disable it please use a value >=0 [-1]
- `threads`: Number of threads to use. >=1 [4]


### Query results

CellPhoneDB results can be queried by making use of the `search_analysis_results` method. This method requires two of the files generated by CellPhoneDB; `significant_means` and  `deconvoluted`. 

Through this method, users can specify the cell pairs of interest and both; the genes `query_genes` participating in the interaction and/or the name of the interaction itself `query_interactions`. This method will search for significant/relevant interactions in which any cell specified in `query_cell_types_1` is found to any cell specified in `query_cell_types_2`. Cell pairs within any of these two lists will not be queried, that is to say, no interaction between cells A and B or C and D will be queried.

```shell
from cellphonedb.utils import search_utils

search_results = search_utils.search_analysis_results(
    query_cell_types_1 = list_of_cells_1,
    query_cell_types_2 = list_of_cells_2,
    query_genes = list_of_genes,
    query_interactions = list_of_interaction_names,
    significant_means = significant_means,
    deconvoluted = cpdb_deconvoluted,
    separator = '|',
    long_format = True
)

Examples of this are provided in the [tutorials](notebooks).
```
## Plotting results

Currently CellPhoneDB relies on external plotting implementations to represent the results. Some examples are provided in the [tutorials](notebooks).

Currently we recommend using tools such as: seaborn, ggplot or a more specific and tailored implementation as the ktplots:
[@zktuong](https://github.com/zktuong):
- [ktplots](https://www.github.com/zktuong/ktplots/) (R)
- [ktplotspy](https://www.github.com/zktuong/ktplotspy/) (python)


## Using different database versions
CellPhoneDB databases can be updated from the remote repository through our tool. Furthermore, available versions can be listed and downloaded for use. Please, refer to our tutorials for a comprehensive [example](notebooks).

First, the user must download the database to its preferred directory, once this is done, the user must provide the argument `cpdb_file_path` to the CellPhoneDB method to be executed with the provided version of the database.

The database is downloaded in a `zip` format along with the input files employed to generate it. These input files can be modified to update the database with new interactions.

## Listing remote available versions
The command to list available versions from the remote repository is:
```shell
from IPython.display import HTML, display
from cellphonedb.utils import db_releases_utils

display(HTML(db_releases_utils.get_remote_database_versions_html()['db_releases_html_table']))
``` 
See [examples](notebooks).
## Download version
The command to download a version from the remote repository is:
```shell
from cellphonedb.utils import db_utils

db_utils.download_database(cpdb_target_dir, cpdb_version)
``` 
See [examples](notebooks).

## Generating user-specific custom database
A user can generate custom databases and use them. In order to generate a new database, a user can provide his/her own lists.

We recommend first to download CellPhoneDB database, move these downloaded files into a new folder and then modify its content to add new interactions. Once this process is completed, the `created_db` method will create a new database in `zip` format in the same folder where the inputs are located. Examples of how to download and create the database can be found here [example](notebooks).

To generate such a database the user has to issue this command:
```shell
from cellphonedb.utils import db_utils

db_utils.create_db(cpdb_input_dir) 
```

Result database file is generated in the folder `out` with `cellphonedb_{datetime}.zip`.
Do not change the name of the input files, otherwise CellPhoneDB will not recognize them and and error will be thrown.


## Contributing to CellPhoneDB

CellPhoneDB is an open-source project. If you are interested in contributing to this project, please let us know.

You can check all project documentation in the [docs](https://cellphonedb.readthedocs.io/en/latest/#) section


## Citing CellphoneDB

**CellphoneDB v4**: Single-cell roadmap of human gonadal development. L Garcia-Alonso, V Lorenzi et al. 2022 Nature [link](https://www.nature.com/articles/s41586-022-04918-4)

**CellphoneDB v3**: Mapping the temporal and spatial dynamics of the human endometrium in vivo and in vitro. L Garcia-Alonso, L-François Handfield, K Roberts, K Nikolakopoulou et al. Nature Genetics 2021 [link](https://www.nature.com/articles/s41588-021-00972-2)

**CellPhoneDB v2**: Inferring cell-cell communication from combined expression of multi-subunit receptor-ligand complexes. Efremova M, Vento-Tormo M, Teichmann S, Vento-Tormo R. Nat Protoc. 2020 [link](https://www.nature.com/articles/s41596-020-0292-x)

**CellPhoneDB v1 original**: Single-cell reconstruction of the early maternal-fetal interface in humans. Vento-Tormo R, Efremova M, et al., Nature. 2018 [link](https://www.nature.com/articles/s41586-018-0698-6)

The first version of CellphoneDB was developed at the [Teichmann Lab](http://www.teichlab.org/) by Roser Vento-Tormo and Mirjana Efremova (CellphoneDB v1 and v2). Currently, it is being further developed and supported by the [Vento-Tormo Lab](https://ventolab.org/) (Wellcome Sanger Institute, Cambridge, UK; CellphoneDB ≥v3).
