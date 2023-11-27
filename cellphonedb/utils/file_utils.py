import os
from datetime import datetime
import pickle
from typing import TextIO, Optional, Tuple
import csv
import scipy.io
import io
import zipfile
from anndata import read_h5ad, AnnData
import pandas as pd

from cellphonedb.src.exceptions.NotADataFrameException import NotADataFrameException
from cellphonedb.src.exceptions.ReadFileException import ReadFileException
from cellphonedb.src.exceptions.ReadFromPickleException import ReadFromPickleException
from cellphonedb.src.core.preprocessors import method_preprocessors, counts_preprocessors

DEBUG = False


def read_data_table_from_file(file: str, index_column_first: bool = False, separator: str = '',
                              dtype=None, na_values=None, compression=None, optional=False) -> pd.DataFrame:
    if os.path.isdir(file):
        return _read_mtx(file)

    filename, file_extension = os.path.splitext(file)

    if file_extension == '.h5ad':
        return _read_h5ad(file)
    if file_extension == '.h5':
        return _read_h5(file)
    if file_extension == '.pickle':
        return _read_pickle(file)

    if not separator:
        separator = _get_separator(file_extension)
    try:
        f = open(file)
    except Exception:
        if not optional:
            raise ReadFileException(file)
        else:
            # Don't raise an exception if the file is optional (i.e. present in only certain DB versions)
            return None
    else:
        with f:
            return _read_data(f, separator, index_column_first, dtype, na_values, compression)


def write_to_file(df: pd.DataFrame, filename: str, output_path: str,
                  output_format: Optional[str] = None, index_label=None, index=False):
    _, file_extension = os.path.splitext(filename)

    if output_format is None:
        if not file_extension:
            default_format = 'txt'
            default_extension = '.{}'.format(default_format)

            separator = _get_separator(default_extension)
            filename = '{}{}'.format(filename, default_extension)
        else:
            separator = _get_separator(file_extension)
    else:
        selected_extension = '.{}'.format(output_format)

        if file_extension != selected_extension:
            separator = _get_separator(selected_extension)
            filename = '{}{}'.format(filename, selected_extension)

            if file_extension:
                print('WARNING: Selected extension missmatches output filename ({}, {}): It will be added => {}'.format(
                        selected_extension, file_extension, filename))
        else:
            separator = _get_separator(selected_extension)

    df.to_csv('{}/{}'.format(output_path, filename), sep=separator, index=index, index_label=index_label)


def _read_mtx(path: str) -> pd.DataFrame:

    mtx_path = os.path.join(path, 'matrix.mtx')
    bc_path = os.path.join(path, 'barcodes.tsv')
    feature_path = os.path.join(path, 'features.tsv')

    df = pd.DataFrame(scipy.io.mmread(mtx_path).toarray())
    with open(bc_path) as bc_file:
        df.columns = [bc[0].strip() for bc in list(csv.reader(bc_file, delimiter="\t"))]
    with open(feature_path) as feature_file:
        df.index = [feat[0].strip() for feat in list(csv.reader(feature_file, delimiter="\t"))]
    df.index.name = 'Gene'

    return df


def _read_h5ad(path: str) -> pd.DataFrame:
    adata = read_h5ad(path)
    df = adata.to_df().T
    return df


def _read_h5(path: str) -> pd.DataFrame:
    df = pd.read_hdf(path)
    return df


def _read_pickle(path: str) -> pd.DataFrame:
    try:
        with open(path, 'rb') as f:
            df = pickle.load(f)
            if isinstance(df, pd.DataFrame):
                return df
            else:
                raise NotADataFrameException(path)
    except Exception:
        raise ReadFromPickleException(path)


def _read_data(file_stream: TextIO, separator: str, index_column_first: bool, dtype=None,
               na_values=None, compression=None) -> pd.DataFrame:
    return pd.read_csv(file_stream, sep=separator, index_col=0 if index_column_first else None, dtype=dtype,
                       na_values=na_values, compression=compression)


def set_paths(output_path, project_name):
    if project_name:
        output_path = os.path.realpath(os.path.expanduser('{}/{}'.format(output_path, project_name)))

    os.makedirs(output_path, exist_ok=True)

    # if _path_is_not_empty(output_path):
    #     print(
    #         'WARNING: Output directory ({}) exist and is not empty. Result can overwrite old results'.format(output_path))

    return output_path


def _path_is_not_empty(path):
    return bool([f for f in os.listdir(path) if not f.startswith('.')])


def _get_separator(mime_type_or_extension: str) -> str:
    extensions = {
        '.csv': ',',
        '.tsv': '\t',
        '.txt': '\t',
        '.tab': '\t',
        'text/csv': ',',
        'text/tab-separated-values': '\t',
    }
    default_separator = ','

    return extensions.get(mime_type_or_extension.lower(), default_separator)


# From interaction_properties.py
def is_cellphonedb_interactor(interaction: pd.Series, suffixes=('_1', '_2')) -> bool:
    if interaction['annotation_strategy'] == 'curated':
        return True

    if interaction['annotation_strategy'] == 'user_curated':
        return True

    if interaction['id_multidata{}'.format(suffixes[0])] == interaction['id_multidata{}'.format(suffixes[1])]:
        return False

    if interaction['annotation_strategy'] == 'guidetopharmacology.org':
        return True

    if can_be_receptor(interaction, suffixes[0]) and \
            can_be_ligand(interaction, suffixes[1]):
        return True

    if can_be_receptor(interaction, suffixes[1]) and \
            can_be_ligand(interaction, suffixes[0]):
        return True

    return False


# From multidata_properties.py
def can_be_receptor(multidata: pd.Series, suffix: str = '') -> bool:
    if multidata['receptor{}'.format(suffix)] and \
            not multidata['other{}'.format(suffix)]:
        return True
    return False


# From multidata_properties.py
def can_be_ligand(multidata: pd.Series, suffix: str = '') -> bool:
    if multidata['secreted_highlight{}'.format(suffix)]:
        return True
    return False


def dbg(*argv):
    if DEBUG:
        for arg in argv:
            print(arg)


def write_to_csv(rows, file_path, delimiter=','):
    with open(file_path, 'w') as f:
        writer = csv.writer(f, delimiter=delimiter, quoting=csv.QUOTE_NONE, escapechar='\\')
        for row in rows:
            writer.writerow(row)


def get_counts_meta_adata(counts_fp, meta_fp) -> AnnData:
    filename, file_extension = os.path.splitext(counts_fp)

    if file_extension == '.h5ad':
        adata = read_h5ad(counts_fp)
    elif file_extension == '.txt':
        df = read_data_table_from_file(counts_fp, index_column_first=True)
        obs = pd.DataFrame()
        obs.index = df.columns
        var = pd.DataFrame(index=df.index)
        adata = AnnData(df.T.values, obs=obs, var=var, dtype='float64')

    raw_meta = read_data_table_from_file(meta_fp, index_column_first=False)
    meta = method_preprocessors.meta_preprocessor(raw_meta)
    adata.obs = meta

    return adata


def estimate_memory_for_analyses(meta_fp, subsampling=False, scoring=False, num_cores=1):
    raw_meta = read_data_table_from_file(meta_fp, index_column_first=False)
    num_cells = raw_meta.shape[0]
    if scoring:
        statistical_analysis_mem = num_cells * 0.5 + num_cores * 1024 * 1.5
    else:
        if subsampling:
            statistical_analysis_mem = num_cells * 0.5
        else:
            statistical_analysis_mem = 0.3 * num_cells + num_cores * 1024 * 0.5
    statistical_analysis_mem_in_gb = round(statistical_analysis_mem/1024, 0)
    basic_deg_analysis_mem_in_gb = round(0.3 * num_cells / 1024, 0)
    print("Basic or DEG analysis: " + str(basic_deg_analysis_mem_in_gb) + " GB")
    print("Statistical analysis: " + str(statistical_analysis_mem_in_gb) + " GB")


def get_timestamp_suffix():
    return datetime.now().strftime("%m_%d_%Y_%H%M%S")


def save_dfs_as_tsv(out, suffix, analysis_name, name2df):
    if suffix is None:
        suffix = get_timestamp_suffix()
    os.makedirs(out, exist_ok=True)
    for name, df in name2df.items():
        if not df.empty:
            file_path = os.path.join(out, "{}_{}_{}.{}".format(analysis_name, name, suffix, "txt"))
            df.to_csv(file_path, sep='\t', index=False)
            print("Saved {} to {}".format(name, file_path))


def save_scored_interactions_as_zip(out, suffix, analysis_name, interaction_scores_dict):
    name = "interaction_scores"
    if suffix is None:
        suffix = get_timestamp_suffix()
    os.makedirs(out, exist_ok=True)
    zip_buffer = io.BytesIO()
    with zipfile.ZipFile(zip_buffer, "a",
                         zipfile.ZIP_DEFLATED, False) as zip_file:
        for ctPair, df in interaction_scores_dict.items():
            zip_file.writestr('{}.csv'.format(ctPair), df.to_csv(index=False, sep=',').encode('utf-8'))
    file_path = os.path.join(out, '{}_{}_{}.{}'.format(analysis_name, name, suffix, "zip"))
    with open(file_path, 'wb') as f:
        f.write(zip_buffer.getvalue())
    print("Saved {} to {}".format(name, file_path))


def _load_microenvs(microenvs_filepath: str, meta: pd.DataFrame) -> pd.DataFrame:
    """Load microenvironment file

    This method reads a microenvironment file into a DataFrame.
    Runs validations to make sure the file has enough columns and
    that all the cell types in the microenvironment file are included in metadata file.

    Parameters
    ----------
    microenvs_filepath
        Path to the microenvironments file.
    meta
        Meta DataFrame.

    Returns
    -------
    pd.DataFrame
        Microenvironments as a DataFrame read from the input file.
    """
    CELL_TYPE = "cell_type"
    MICRO_ENVIRONMENT = "microenvironment"
    microenvs = read_data_table_from_file(os.path.realpath(microenvs_filepath))
    microenvs.drop_duplicates(inplace=True)
    len_columns = len(microenvs.columns)
    if len_columns < 2:
        raise Exception(f"Missing columns in microenvironments: 2 required but {len_columns} provided")
    elif len_columns > 2:
        print(f"WARNING: Microenvironments file expects 2 columns and got {len_columns}. Dropping extra columns.")
    microenvs = microenvs.iloc[:, 0:2]
    if any(~microenvs.iloc[:, 0].isin(meta.iloc[:, 0])):
        raise Exception("Some clusters/cell_types in microenvironments file are not present in metadata")
    microenvs.columns = [CELL_TYPE, MICRO_ENVIRONMENT]
    return microenvs


def _load_active_tfs(active_tfs_filepath: str, meta: pd.DataFrame) -> dict:
    """Load Active TFs file

    This method reads an Active TFs file into a DataFrame.
    Runs validations to make sure the file has enough columns and
    that all the cell types in the active TFs file are included in metadata file.

    Parameters
    ----------
    active_tfs_fi
        Path to the active TFs file.
    meta
        Meta DataFrame.

    Returns
    -------
    pd.DataFrame
        Active TFs per cell type as a DataFrame read from the input file.
    """
    active_tfs = read_data_table_from_file(os.path.realpath(active_tfs_filepath))
    active_tfs.drop_duplicates(inplace=True)
    len_columns = len(active_tfs.columns)
    if len_columns < 2:
        raise Exception(f"Missing columns in Active TFs file: 2 required but {len_columns} provided")
    elif len_columns > 2:
        print(f"WARNING: Active TFs file expects 2 columns and got {len_columns}. Dropping extra columns.")
    active_tfs = active_tfs.iloc[:, 0:2]
    rows_filter = ~active_tfs.iloc[:, 0].isin(meta.iloc[:, 0])
    if any(rows_filter):
        print("WARNING: The following clusters/cell_types in Active TFs file are not present in metadata:")
        print("\n".join(set(active_tfs[rows_filter].iloc[:, 0].tolist())))

    # Convert DataFrame to dict
    active_tf2cell_types = {}
    for cell_type, tf in active_tfs.values:
        active_tf2cell_types.update({tf: active_tf2cell_types.get(tf, []) + [cell_type]})

    return active_tf2cell_types


def _load_degs(degs_filepath: str, meta: pd.DataFrame) -> pd.DataFrame:
    """Load DEGs file

    This method reads a DEGs file into a DataFrame. Runs validations
    to make sure the file has enough columns and that all the clusters
    in DEGs are included in meta.

    Parameters
    ----------
    degs_filepath
        Path to the DEGs file.
    meta
        Meta DataFrame.

    Returns
    -------
    DataFrame
        DEGs as a DataFrame with cluster and genes columns.
    """
    CLUSTER = "cluster"
    GENE = "gene"
    degs_filepath = os.path.realpath(degs_filepath)
    degs = read_data_table_from_file(degs_filepath)
    len_columns = len(degs.columns)
    if len_columns < 2:
        raise Exception(f"Missing columns in DEGs: 2 required but {len_columns} provided")
    elif len_columns > 2:
        print(f"WARNING: DEGs expects 2 columns and got {len_columns}. Dropping extra columns.")
    degs = degs.iloc[:, 0:2]
    if any(~degs.iloc[:, 0].isin(meta.iloc[:, 0])):
        raise Exception("Some clusters/cell_types in DEGs are not present in metadata")
    degs.columns = [CLUSTER, GENE]
    degs.drop_duplicates(inplace=True)
    return degs


def is_anndata_type(obj):
    """ Return true if obj is of type AnnData; false otherwise"""
    return type(obj).__name__ == "AnnData"


def read_counts(obj) -> (pd.DataFrame, str):
    """ Read counts df from obj - obj can be AnnData or a path to a file (various formats allowed) that contains counts """
    if is_anndata_type(obj):
        counts = obj.to_df().T
        counts_label = "counts from AnnData object"
    else:
        counts = read_data_table_from_file(obj, index_column_first=True)
        counts_label = obj
    return counts, counts_label


def get_user_files(counts=None, meta_fp=None, microenvs_fp=None, degs_fp=None, active_tfs_fp=None,
                   gene_synonym2gene_name=None, counts_data=None) \
        -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame, dict]:
    """

    Parameters
    ----------
    counts
        Path to the user's counts file, exemplified by \
        https://github.com/ventolab/CellphoneDB/blob/bare-essentials/example_data/test_counts.txt,
        or an in-memory AnnData object
    meta_fp
        Path to the user's meta file, exemplified by \
        https://github.com/ventolab/CellphoneDB/blob/bare-essentials/example_data/test_meta.txt
    microenvs_fp
        Path to the user's microenvironments file, exemplified by \
        https://github.com/ventolab/CellphoneDB/blob/bare-essentials/example_data/test_microenviroments.txt
    degs_fp
        Path to the user's differentially expresses genes (DEGs) file, exemplified by \
        https://github.com/ventolab/CellphoneDB/blob/bare-essentials/example_data/test_degs.txt

    Returns
    -------
    Tuple
        - counts: pd.DataFrame
        - meta: pd.DataFrame
        - microenvs: pd.DataFrame
        - degs: pd.DataFrame
        - active_tf2cell_types: dict

    """
    loaded_user_files = []
    # Read user files
    print("Reading user files...", flush=True)
    counts, counts_label = read_counts(counts)

    # N.B. The functionality below has been switched off for the time being, on Kevin's request
    # In counts df, replace any gene synonyms not in gene_input.csv to gene names that are in gene_input.
    # if counts_data == "hgnc_symbol" or counts_data == "gene_name":
    #     counts.rename(index=gene_synonym2gene_name, inplace=True)

    loaded_user_files.append(counts_label)
    raw_meta = read_data_table_from_file(meta_fp, index_column_first=False)
    meta = method_preprocessors.meta_preprocessor(raw_meta)
    loaded_user_files.append(meta_fp)
    # Ensure that counts values are of type float32, and that all cells in meta exist in counts
    counts = counts_preprocessors.counts_preprocessor(counts, meta)
    if microenvs_fp:
        microenvs = _load_microenvs(microenvs_fp, meta)
        loaded_user_files.append(microenvs_fp)
    else:
        microenvs = pd.DataFrame()

    if active_tfs_fp:
        active_tf2cell_types = _load_active_tfs(active_tfs_fp, meta)
        loaded_user_files.append(active_tfs_fp)
    else:
        active_tf2cell_types = {}

    if degs_fp:
        degs = _load_degs(degs_fp, meta)
        loaded_user_files.append(degs_fp)
    else:
        degs = pd.DataFrame()

    print("The following user files were loaded successfully:")
    for fn in loaded_user_files:
        print(fn)

    return counts, meta, microenvs, degs, active_tf2cell_types
