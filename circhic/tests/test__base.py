import pytest
from circhic import CircHiCFigure
from circhic import datasets
import numpy as np


def test_plot_raxis():
    fig = CircHiCFigure(lengths=np.array([78]))
    with pytest.raises(ValueError):
        fig.plot_raxis()

    data = datasets.load_ccrescentus()
    counts = data["counts"]
    lengths = data["nbins"]

    # Some small flush test
    fig = CircHiCFigure(lengths=lengths)
    fig.plot_hic(counts)
    fig.plot_raxis()
