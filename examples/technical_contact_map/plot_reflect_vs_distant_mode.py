"""
===========================
Reflect versus distant mode
===========================

In this example, we showcase the difference between the reflect versus the
distant mode.
"""
from iced.normalization import ICE_normalization

from circhic import datasets
from circhic._base import CircHiCFigure


# Load the data, compute the cumulative raw counts.
data = datasets.load_ecoli()
counts = data["counts"]
lengths = data["lengths"]

# Normale the data using ICE, and keep the biases
counts, bias = ICE_normalization(counts, output_bias=True)

###############################################################################
# The reflect mode (the default) will plot data corresponding to the upper and
# lower triangular of the original contact count matrix
circhicfig = CircHiCFigure(lengths)
im, ax = circhicfig.plot_hic(counts, mode="reflect", inner_gdis=50,
                             inner_radius=0.3)

###############################################################################
# The distant mode will plot the data between inner_gdis and outer_gdis, thus
# ignoring the section close to the diagonal.
circhicfig = CircHiCFigure(lengths)
im, ax = circhicfig.plot_hic(counts, mode="distant", inner_gdis=50,
                             inner_radius=0.3)
