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
circfig.plot_hic(counts, inner_gdis=80, outer_gdis=80, cmap="bone_r",
                 inner_radius=0.5)
circfig.plot_raxis()
