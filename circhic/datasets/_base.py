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

    Example
    -------

    loading the data

    >> from circhic import datasets
    >> data = datasets.load_bsubtilis()
    >> print(data)
        {'counts': array([[  62.,  469.,  457., ...,  382.,  701., 2311.],
            [ 469., 4908., 1245., ...,  362.,  642., 1227.],
            [ 457., 1245., 5940., ...,  487.,  753., 1180.],
            ...,
            [ 382.,  362.,  487., ..., 2740., 1496., 1210.],
            [ 701.,  642.,  753., ..., 1496., 3778., 2406.],
            [2311., 1227., 1180., ..., 1210., 2406., 5244.]]),
        'lengths': array([412])}
    """
    module_path = os.path.dirname(__file__)
    lengths = _load_lengths(
        os.path.join(module_path,
                     "data/bsubtilis/SRX1014144_9790_abs.bed"))

    counts = _load_counts(
        os.path.join(module_path,
                     "data/bsubtilis/SRX1014144_9790.matrix"),
        lengths=lengths)
    counts = counts.toarray()
    counts = counts.T + counts
    results = {"counts": counts,
               "lengths": lengths}

    return results


def load_ecoli():
    """
    Loads *E. coli* contact counts.

    Returns
    -------
    dictionary :
        a dictionary containing:
            - counts: an (n, n) ndarray corresponding to the raw contact
              counts for *B. subtilis*
            - lengths: (l, ) ndarray containing the lengths of all chromosomes.
    """
    module_path = os.path.dirname(__file__)
    lengths = _load_lengths(
        os.path.join(module_path,
                     "data/ecoli/SRX3451210_9897_abs.bed"))
    counts = _load_counts(
        os.path.join(module_path,
                     "data/ecoli/SRX3451210_9897.matrix"),
        lengths=lengths)
    counts = counts.toarray()
    counts = counts.T + counts
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


def _load_lengths(filename, return_base=False):
    """
    Fast loading of the bed files

    Parameters
    ----------
    filename : str,
        path to the file to load. The file should be a bed file

    return_base : bool, optional, default: False
        whether to return if it is 0 or 1-base

    Returns
    -------
    lengths : the lengths of each chromosomes
    """
    data = pd.read_csv(filename, sep="\t", comment="#", header=None)
    data = data.values
    _, idx, lengths = np.unique(data[:, 0], return_counts=True,
                                return_index=True)
    if return_base:
        return lengths[idx.argsort()], data[0, 3]
    else:
        return lengths[idx.argsort()]
