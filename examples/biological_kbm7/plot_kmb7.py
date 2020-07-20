"""
================================
A human chromosome: KBM7's chr11
================================

Here, we plot the raw contact counts of KBM7's chromosome 11, and the
cumulative raw contact counts.
"""

from circhic import CircHiCFigure, datasets

data = datasets.load_kbm7()
counts = data["counts"]
lengths = data["nbins"]

cum_raw_counts = counts.sum(axis=1)

circhicfig = CircHiCFigure(lengths, chromosome_type="linear")
circhicfig.plot_hic(counts, inner_radius=0.5, outer_radius=0.89, cmap="bone_r",
                    inner_gdis=80, outer_gdis=80)
_, ax = circhicfig.plot_bars(
    cum_raw_counts,
    inner_radius=0.9,
    outer_radius=0.99,
    color="0")
circhicfig.set_genomic_ticklabels()
