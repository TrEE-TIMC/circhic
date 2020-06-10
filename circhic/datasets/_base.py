import os
import warnings
import numpy as np
from scipy import sparse
import pandas as pd


def load_ccrescentus():
    """
    Loads *C. crescentus* contact counts.

    Returns
    -------
    dictionary :
        a dictionary containing:
            - counts: an (n, n) ndarray corresponding to the raw contact
        counts for *C. crescentus*
            - lengths: (l, ) ndarray containing the lengths of all chromosomes.
    """
    module_path = os.path.dirname(__file__)
    counts = _load_counts(
        os.path.join(module_path,
                     "data/ccrescentus/SRX263925_9958.matrix"))
    counts = counts.toarray()
    counts = counts.T + counts
    lengths = np.array([counts.shape[0]])

    results = {"counts": counts,
               "lengths": lengths}
    return results


def load_bsubtilis():
    """
    Loads *B. subtilis* contact counts.

    Returns
    -------
    dictionary :
        a dictionary containing:
            - counts: an (n, n) ndarray corresponding to the raw contact
              counts for *B. subtilis*
            - lengths: (l, ) ndarray containing the lengths of all chromosomes.
    """
    module_path = os.path.dirname(__file__)
    counts = _load_counts(
        os.path.join(module_path,
                     "data/bsubtilis/SRX1014144_9790.matrix"))
    counts = counts.toarray()
    counts = counts.T + counts
    lengths = np.array([counts.shape[0]])
    results = {"counts": counts,
               "lengths": lengths}

    return results


def _load_counts(filename, lengths=None):
    """
    Fast loading of a raw interaction counts file

    Parameters
    ----------
    filename : str
        path to the file to load. The file should be of the following format:
        i, j, counts

    lengths : ndarray
        lengths of each chromosomes

    Returns
    --------
    X : the interaction counts file
    """
    base = 1
    n = None
    if lengths is not None:
        n = lengths.sum()
        shape = (n, n)
    else:
        shape = None
    # This is the interaction count files
    dataframe = pd.read_csv(filename, sep="\t", comment="#", header=None)
    row, col, data = dataframe.values.T

    # If there are NAs remove them
    mask = np.isnan(data)
    if np.any(mask):
        warnings.warn(
            "NAs detected in %s. "
            "Removing NAs and replacing with 0." % filename)
        row = row[np.invert(mask)]
        col = col[np.invert(mask)]
        data = data[np.invert(mask)]

    col -= base
    row -= base

    if shape is None:
        n = max(col.max(), row.max()) + 1
        shape = (int(n), int(n))

    data = data.astype(float)
    counts = sparse.coo_matrix((data, (row, col)), shape=shape)
    return counts
