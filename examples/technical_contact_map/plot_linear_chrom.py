"""
===========================
Plotting linear chromosomes
===========================
"""

from circhic import datasets, CircHiCFigure

data = datasets.load_kbm7()
counts = data["counts"]
nbins = data["nbins"]


circfig = CircHiCFigure(lengths=nbins, chromosome_type="linear")
circfig.plot_hic(counts)
circfig.plot_raxis()
