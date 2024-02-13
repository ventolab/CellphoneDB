[![Database](https://img.shields.io/github/v/release/ventolab/cellphonedb-data.svg?color=blue&label=database)](https://github.com/ventolab/CellphoneDB-data) [![Python package](https://img.shields.io/pypi/v/cellphonedb.svg?color=brightgreen&label=python-package)](https://pypi.org/project/cellphonedb)
# CellphoneDB 

## What is CellphoneDB?

CellphoneDB is a publicly available repository of **HUMAN** curated receptors, ligands and their interactions paired with a tool to interrogate your own single-cell transcriptomics data (or even bulk transcriptomics data if your samples represent pure populations!). 

> A distinctive feature of CellphoneDB is that the subunit architecture of either ligands and receptors is taken into account, representing heteromeric complexes accurately. This is crucial, as cell communication relies on multi-subunit protein complexes that go beyond the binary representation used in most databases and studies. CellphoneDB also incorporates biosynthetic pathways in which we use the last representative enzyme as a proxy of ligand abundance, by doing so, we include interactions involving non-peptidic molecules. CellphoneDB includes only manually curated & reviewed molecular interactions with evidenced role in cellular communication.

## Documentation
> For more details on using CellphoneDB and a more detailed description of the methods, visit the [DOCUMENTATION](https://cellphonedb.readthedocs.io/en/latest/#). 


## Novel features in v5
1) New python package that can be easily executed in Jupyter Notebook and Collabs. 
2) A scoring methodology to rank interaction based on the expression specificity of the interacting partners.
3) A CellSign module to leverage interactions based on the activity of the transcription factor downstream the receptor. This module is accompanied by a collection of 211 well described receptor-transcription factor direct relationships.
4) A new method of querying of CellphoneDB results `search_utils.search_analysis_results`.
5) Tutorials to run CellphoneDB (available [here](https://github.com/ventolab/CellphoneDB/tree/master/notebooks))
6) Improved computational efficiency of method 2 `cpdb_statistical_analysis_method`.
7) A new database ([cellphonedb-data v5.0](https://github.com/ventolab/cellphonedb-data)) with more manually curated interactions, making up to a total of ~3,000 interactions. This release of CellphoneDB database has three main changes:
    - Integrates new manually reviewed interactions with evidenced roles in cell-cell communication. 
    - Includes non-protein molecules acting as ligands.
    - For interactions with a demonstrated signalling directionality, partners have been ordered according (ligand is partner A, receptor partner B).
    - Interactions have been classified within signaling pathways.
    - CellphoneDB no longer imports interactions from external resources. This is to avoid the inclusion of low-confidence interactions.

See updates from [previous releases here](https://github.com/ventolab/CellphoneDB/blob/master/docs/RESULTS-DOCUMENTATION.md#release-notes).


## Installing CellphoneDB 

We highly recommend using an isolated python environment (as described in steps 1 and 2) using [conda](https://docs.conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html) or [virtualenv](https://docs.python.org/3/library/venv.html) but you could of course omit these steps and install via `pip` immediately.

1. Create python=>3.8 environment
   - Using conda: `conda create -n cpdb python=3.8`
   - Using virtualenv: `python -m venv cpdb`

2. Activate environment
   - Using conda: `source activate cpdb`
   - Using virtualenv: `source cpdb/bin/activate`

3. Install CellphoneDB `pip install cellphonedb`

4. Set up the kernel for the Jupyter notebooks.
   - Install the ipython kernel: `pip install -U ipykernel`.
   - Add the environment as a jupyter kernel: `python -m ipykernel install --user --name 'cpdb'`.
   - Open/Start Jupyter and select the created kernel.

5. Download the database.
   - Follow this [tutorial](https://github.com/ventolab/CellphoneDB/blob/master/notebooks/T0_DownloadDB.ipynb).

NOTE: Works with Python v3.8 or greater. If your default Python interpreter is for `v2.x` (you can check it with `python --version`), calls to `python`/`pip` should be substituted by `python3`/`pip3`.

#### Running CellphoneDB Methods

Please, activate your environment if you didn't previously
- Using conda: `source activate cpdb`
- Using virtualenv: `source cpdb/bin/activate`

We have created a set of tutorials that can be accessed for each To use the example data, please [tutorials and data](notebooks).

## Preparing INPUTS
#### Preparing your counts input file (mandatory)
Counts file can be a text file or a `h5ad` (recommended), `h5` or a path to a folder containing a 10x output with `mtx/barcode/features` files. NOTE: Your gene/protein **ids must be HUMAN**. If you are working with another specie such as mouse, we recommend you to convert the gene ids to their corresponding orthologous. 

#### Preparing your DEGs file (optional, if `method degs_analysis`)
This is a two columns file indicating which gene is specific or up-regulated in a cell type (see [example](notebooks/data_tutorial.zip) ). The first column should be the cell type/cluster name (matching those in `meta.txt`) and the second column the associated gene id. The remaining columns are ignored. We provide [notebooks](notebooks) for both Seurat and Scanpy users. It is on you to design a DEG analysis that is appropriate to your research question. 

#### Preparing your microenvironments file (optional, if `microenvs_file_path`)
This is a two columns file indicating which cell type is in which spatial microenvironment (see [example](notebooks/data_tutorial.zip) ). CellphoneDB will use this information to define possible pairs of interacting cells (i.e. pairs of clusters co-appearing in a microenvironment). 

#### Preparing your active transcription factor file (optional, if `active_tfs_file_path`)
This is a two columns file indicating which cell type and which TFs are active (see [example](notebooks/data_tutorial.zip)). CellphoneDB will use this information to highlight relevant/significant interactions whose downstream TF is active. The information defined in this file (which TFs are active per cell) must be provided by the user.

### RUN examples

For more detailed examples refer to our tutorials [here](notebooks).

####  Example with running the DEG-based method
```shell
from cellphonedb.src.core.methods import cpdb_degs_analysis_method

cpdb_results = cpdb_degs_analysis_method.call(
        cpdb_file_path = cellphonedb.zip,
        meta_file_path = test_meta.txt,
        counts_file_path = test_counts.h5ad,
        degs_file_path = degs_file_path,
        counts_data = 'hgnc_symbol',
        active_tfs_file_path = active_tf.txt,
        score_interactions = True,
        microenvs_file_path = microenvs_file_path,
        threshold = 0.1,
        output_path = out_path)
```

####  Example with running the statistical method
```shell
from cellphonedb.src.core.methods import cpdb_statistical_analysis_method

cpdb_results = cpdb_statistical_analysis_method.call(
        cpdb_file_path = cellphonedb.zip,
        meta_file_path = test_meta.txt,
        counts_file_path = test_counts.h5ad,
        counts_data = 'hgnc_symbol',
        active_tfs_file_path = active_tf.txt,
        microenvs_file_path = microenvs_file_path
        score_interactions = True,
        threshold = 0.1,
        output_path = out_path)
```

#### Example without using the statistical method
 - **Using text files**
```shell
from cellphonedb.src.core.methods import cpdb_analysis_method

cpdb_results = cpdb_analysis_method.call(
        cpdb_file_path = cellphonedb.zip,
        meta_file_path = test_meta.txt,
        counts_file_path = test_counts.txt,
        counts_data = 'hgnc_symbol',
        score_interactions = True,
        output_path = out_path)
```

 - **Using h5ad count file**
```shell
from cellphonedb.src.core.methods import cpdb_analysis_method

cpdb_results = cpdb_analysis_method.call(
        cpdb_file_path = cellphonedb.zip,
        meta_file_path = test_meta.txt,
        counts_file_path = test_counts.h5ad,
        counts_data = 'hgnc_symbol',
        output_path = out_path)
```

Results are saved as files in `output_path` and as a dictionary of dataframes in the output variable `cpdb_results`.

To understand the different analysis and results, please check the [results documentation](docs/RESULTS-DOCUMENTATION.md).


### Optional Parameters

 **Optional Method parameters**:
- `counts_data`: Type of gene identifiers in the counts data [ensembl \| gene_name \| hgnc_symbol]
- `iterations`: Number of iterations for the statistical analysis [1000]
- `threshold`: % of cells expressing the specific ligand/receptor
- `result_precision`: Number of decimal digits in results [3]
- `output_path`: Directory where the results will be saved (the directory must exist) [out]
- `output_suffix`: Output format of the results files (time stamp will be added to filename if not present) [txt]
- `subsampling`: Enable subsampling based on geometric sketching
- `subsampling_log`: Enable subsampling log1p for non log-transformed data inputs !!mandatory!!
- `subsampling_num_pc`: Subsampling NumPC argument (number of PCs to use) [100]
- `subsampling_num_cells`: Number of cells to subsample the dataset [1/3 of cells by default]


~ **Optional Method Statistical parameters**
- `microenvs_file_path`: Spatial microenvironments input file. Restricts the cluster/cell_type interacting pairs to the cluster/cell_type sharing a microenviroment (i.e. only test a combination of clusters if these coexist in a microenviroment). This file should contain two columns: 1st column indicates the cluster/cell_type, 2nd column indicates the microenviroment name.  See example [here](https://github.com/ventolab/CellphoneDB/tree/master/in). 
- `pvalue`: P-value threshold [0.05]
- `debug_seed`: Debug random seed -1. To disable it please use a value >=0 [-1]
- `threads`: Number of threads to use. >=1 [4]


### Query results

CellphoneDB results can be queried by making use of the `search_analysis_results` method. This method requires two of the files generated by CellphoneDB `significant_means` and  `deconvoluted`, optionally `interaction_scores` can be used to subset interactions by score.

Through this method, users can specify the cell pairs of interest and both; the genes `query_genes` participating in the interaction and/or the name of the interaction itself `query_interactions`. This method will search for significant/relevant interactions in which any cell specified in `query_cell_types_1` is found to any cell specified in `query_cell_types_2`. Cell pairs within any of these two lists will not be queried, that is to say, no interaction between cells A and B or C and D will be queried.

```shell
from cellphonedb.utils import search_utils

search_results = search_utils.search_analysis_results(
        query_cell_types_1 = ['cell_A', 'cell B'],  
        query_cell_types_2 = ['cell_C', 'cell D'],
        query_genes = ['Gene_1', 'Gene_2', 'Gene_3'], 
        query_interactions = ['interaction_name_1', 'interaction_name_2'],  
        significant_means = cpdb_results['significant_means'],
        deconvoluted = cpdb_results['cpdb_deconvoluted'],
        interaction_scores = cpdb_results['interaction_scores'],
        query_minimum_score = 50,
        query_classifications = ['pathway_A', 'pathway_B'],
        separator = '|',
        long_format = True
)
```
Examples provided in [tutorials](notebooks).
## Plotting results

Currently, CellphoneDB relies on external plotting implementations to represent the results. Examples of the use are provided in the [tutorials](notebooks).

We recommend using tools such as: seaborn, ggplot, or a more specific and tailored implementation such as the [@ktplots](https://github.com/zktuong):
- [ktplots](https://www.github.com/zktuong/ktplots/) (R)
- [ktplotspy](https://www.github.com/zktuong/ktplotspy/) (python)


## Using different database versions
CellphoneDB databases can be updated from the remote repository through our tool. Furthermore, available versions can be listed and downloaded for use. Please, refer to our tutorials for a comprehensive [example](notebooks).

First, the user must download the database to its preferred directory, once this is done, the user must provide the argument `cpdb_file_path` to the CellphoneDB method to be executed with the provided version of the database.

The database is downloaded in a `zip` format along with the input files employed to generate it. These input files can be modified to update the database with new interactions.

> CellphoneDB v5 is compatible with database version 4.1.0 or newer.

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

We recommend first to download CellphoneDB database, move these downloaded files into a new folder and then modify its content to add new interactions. Once this process is completed, the `created_db` method will create a new database in `zip` format in the same folder where the inputs are located. Examples of how to download and create the database can be found here [example](notebooks).

To generate such a database the user has to issue this command:
```shell
from cellphonedb.utils import db_utils

db_utils.create_db(cpdb_input_dir) 
```

Result database file is generated in the folder `out` with `cellphonedb_{datetime}.zip`.
Do not change the name of the input files, otherwise CellphoneDB will not recognize them and and error will be thrown.


## Contributing to CellphoneDB

CellphoneDB is an open-source project. If you are interested in contributing to this project, please let us know.
You can check all project documentation in the [docs](https://cellphonedb.readthedocs.io/en/latest/#) section


## Citing CellphoneDB

The first version of CellphoneDB was originally developed at the [Teichmann Lab](http://www.teichlab.org/) in the Wellcome Sanger Institute (Cambridge, UK) by Roser Vento-Tormo and Mirjana Efremova. Currently, it is being further developed and supported by the [Vento-Tormo Lab](https://ventolab.org/) (CellphoneDB ≥v3).

If you use CellphoneDB or CellphoneDB-data, please cite our papers:

- **CellphoneDB v1 (original)**: Single-cell reconstruction of the early maternal-fetal interface in humans. Vento-Tormo R, Efremova M, et al., Nature. 2018 [link](https://www.nature.com/articles/s41586-018-0698-6)

- **CellphoneDB v2**: Inferring cell-cell communication from combined expression of multi-subunit receptor-ligand complexes. Efremova M, Vento-Tormo M, Teichmann S, Vento-Tormo R. Nat Protoc. 2020 [link](https://www.nature.com/articles/s41596-020-0292-x)

- **CellphoneDB v3**: Mapping the temporal and spatial dynamics of the human endometrium in vivo and in vitro. L Garcia-Alonso, L-François Handfield, K Roberts, K Nikolakopoulou et al. Nature Genetics 2021 [link](https://www.nature.com/articles/s41588-021-00972-2)

- **CellphoneDB v4**: Single-cell roadmap of human gonadal development. L Garcia-Alonso, V Lorenzi et al. 2022 Nature [link](https://www.nature.com/articles/s41586-022-04918-4)

- **CellphoneDB v5 (latest)**: CellphoneDB v5: inferring cell-cell communication from single-cell multiomics data. (in preparation)