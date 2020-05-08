import os
import numpy as np
from iced import io


def load_ccrescentus():
    """
    Loads *C. crescentus* contact counts.

    Returns
    -------

    counts, lengths :
        a tuple containing an (n, n) ndarray corresponding to the raw contact
        counts for *C. crescentus*, and an (l, ) ndarray containing the
        lengths of all chromosomes.

    """
    module_path = os.path.dirname(__file__)
    counts = io.load_counts(
        os.path.join(module_path,
                     "data/ccrescentus/SRX263925_10000.matrix"),
        base=1)
    counts = counts.toarray()
    counts = counts.T + counts
    lengths = np.array([counts.shape[0]])
    return counts, lengths
