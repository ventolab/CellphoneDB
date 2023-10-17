import numpy as np
import pandas as pd
from fbpca import pca
import warnings
from geosketch import gs
from cellphonedb.src.core.core_logger import core_logger

warnings.simplefilter("ignore", FutureWarning)


class Subsampler(object):
    def __init__(self, log: bool, num_pc: int = 100, num_cells: int = None, verbose: bool = None, debug_seed: int = None):
        """
        Parameters
        ----------
        log: bool
            If true, each element of counts array will be converted to its natural logarithm before subsampling
        num_pc: int
            Number of principal components to be used during sub-sampling (for more information see:
            https://github.com/brianhie/geosketch)
        num_cells: int
            Number of samples to obtain from counts
        verbose: bool
        debug_seed: bool
            Set to True to obtain the same sub-sampling across different runs
        """

        self.verbose = verbose
        self.log = log
        self.num_pc = num_pc
        self.num_cells = num_cells
        np.random.seed(debug_seed)

    def subsample(self, counts: pd.DataFrame) -> pd.DataFrame:
        """
        Parameters
        ----------
        counts: pd.DataFrame
            counts DataFrame to be sub-sampled
        Returns
        -------
        pd.DataFrame
            Sub-sampled counts using the parameters passed in __init__ method

        """
        input_cells = counts.shape[1]

        if self.num_cells is None:
            self.num_cells = int(input_cells / 3)

        core_logger.info('Subsampling {} to {}'.format(input_cells, self.num_cells))

        counts_t = counts.T

        if self.log:
            pca_input = np.log1p(counts_t)
        else:
            pca_input = counts_t

        try:
            u, s, vt = pca(pca_input.values, k=self.num_pc)
            x_dimred = u[:, :self.num_pc] * s[:self.num_pc]
            sketch_index = gs(x_dimred, self.num_cells, replace=False)
            x_matrix = counts_t.iloc[sketch_index]
        except Exception as e:
            core_logger.warning('Subsampling failed: ignored.')
            if self.verbose:
                core_logger.warning(str(e))
            return counts

        core_logger.info('Done subsampling {} to {}'.format(input_cells, self.num_cells))

        return x_matrix.T
