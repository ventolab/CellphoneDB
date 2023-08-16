## Novel features in v4
1) New python package that can be easily executed in Jupyter Notebook and Collabs. 
2) A new method to ease the query of CellPhoneDB results.
3) Tutorials to run CellPhoneDB (available [here](https://github.com/ventolab/CellphoneDB/tree/master/notebooks))
4) Improved computational efficiency of method 2 `cpdb_statistical_analysis_method`.
5) A new database ([cellphonedb-data v4.1.0](https://github.com/ventolab/cellphonedb-data)) with more manually curated interactions, making up to a total of ~3,000 interactions. This release of CellphoneDB database has three main changes:
    - Integrates new manually reviewed interactions with evidenced roles in cell-cell communication. 
    - Includes non-protein molecules acting as ligands.
    - CellphoneDB does not longer imports interactions from external resources. This is to avoid the inclusion of low-confidence interactions.

See updates from [previous releases here](https://github.com/ventolab/CellphoneDB/blob/master/docs/RESULTS-DOCUMENTATION.md#release-notes).

## New in CellPhoneDB-data v4

This release involves a major [CellphoneDB database](https://github.com/ventolab/CellphoneDB-data) update. We have invested quite some time curating more cell-cell communication interactions validated experimentally. Specifically, we have:

 1. Manually curated more protein-protein interactions involved in cell-cell communication, with special focus on protein acting as heteromeric complexes. The new database includes almost **2,000 high-confidence interactions**, including heteromeric complexes! We believe modelling complexes is key to minimise false positives in the predictions.
 2. Annotated non-peptidic molecules (i.e., not encoded by a gene) acting as ligands. Examples of these include steroid hormones (e.g., estrogen). To do so, we have reconstructed the biosynthetic pathways and used the last representative enzyme as a proxy of ligand abundance. We retrieve this information by manually reviewing and curating relevant literature and peer-reviewed pathway resources such as REACTOME. We include more than **200 interactions involving non-peptidic ligands**!


Check [Garcia-Alonso & Lorenzi et al](https://www.nature.com/articles/s41586-022-04918-4) for an example applying CellphoneDB v4.

## New in CellPhoneDB v3

1. **Incorporate spatial information** CellPhoneDB now allows the incorporation of spatial information of the cells via the `microenvironments` file. This is a two columns file indicating which cell type is in which spatial microenvironment (see [example](https://github.com/ventolab/CellphoneDB/blob/master/in/endometrium_atlas_example/endometrium_example_microenviroments.tsv) ). CellphoneDB will use this information to define possible pairs of interacting cells (i.e. pairs of clusters sharing/coexisting in a microenvironment). You can define microenvironments with prior knowledge, imaging or Visium analysis with [cell2location](https://cell2location.readthedocs.io/en/latest/notebooks/cell2location_short_demo_downstream.html#4.-Identify-groups-of-co-located-cell-types-using-matrix-factorisation).
2. **New analysis method added**,  using Differentially Expressed Genes (`cellphonedb method degs_analysis`) as an alternative to the permutation-based approach (`cellphonedb method statistical_analysis`). This `degs_analysis` approach will select interactions where all the genes are expressed by a fraction of cells above a `--threshold` and at least one gene is a DEG. The user identifies the DEGs using their preferred tool and provides the information to CellphoneDB via text file. The first column should be the cell type/cluster and the second column the associated gene id. The remaining columns are ignored (see [example](https://github.com/ventolab/CellphoneDB/blob/master/in/endometrium_atlas_example/endometrium_example_DEGs.tsv) ). We provide [notebooks](https://github.com/ventolab/CellphoneDB/tree/master/notebooks) for both Seurat and Scanpy users.  
3. **Database update** WNT pathway has been further curated.

Check [Garcia-Alonso et al](https://www.nature.com/articles/s41588-021-00972-2) for an example applying CellphoneDB v3.