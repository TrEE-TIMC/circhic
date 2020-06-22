"""
===============================================
a circHiC figure showing biases and mappability
===============================================

"""
import numpy as np

import matplotlib.pyplot as plt

from circhic import datasets
from circhic import CircHiCFigure


from iced.normalization import ICE_normalization


granularity = 0.5
# Load the data, compute the cumulative raw counts.
data = datasets.load_ccrescentus()
counts = data["counts"]
lengths = data["nbins"]

cumul_raw_counts = counts.sum(axis=0)
# Normale the data using ICE, and keep the biases
counts, bias = ICE_normalization(counts, output_bias=True)

fig = plt.figure(figsize=(6, 6))
circhicfig = CircHiCFigure(lengths, figure=fig)
m, ax = circhicfig.plot_hic(counts, granularity=granularity,
                            outer_radius=0.75, inner_radius=0.1,
                            inner_gdis=120, outer_gdis=60,
                            vmin=77, cmap="bone_r")

circhicfig.plot_raxis()

# Assume you want to plot data from that ranges in a polar plot outside of the
# first one. Then the 0 axis should be at, say, 80% of the axis
bar, _ = circhicfig.plot_bars(
    cumul_raw_counts, inner_radius=0.8, outer_radius=0.9,
    color="0")

# Now, plot a last plot, for the top 10% of the original axes
lines, _ = circhicfig.plot_lines(
    bias, color="0", inner_radius=0.9, outer_radius=1)

cab = circhicfig.set_colorbar(m, orientation="horizontal")
cab.set_label("Normalized contact counts", fontweight="bold",
              fontsize="small")

# Now, try to do a simple "theta-axis" on the outer
ticklabels = [
        "%d kb" % (i * 10) for i in np.arange(0, lengths.sum(), 50)[:-1]]
ticklabels[0] = "ORI"
ax = circhicfig.set_genomic_ticklabels(
    tickpositions=np.arange(0, lengths.sum(), 50)[:-1],
    ticklabels=ticklabels,
    outer_radius=0.95)
ax.tick_params(colors="0.3")
#             verticalalignment="bottom")
fig.legend((bar, lines[0]), ("Mappability", "Bias"), fontsize="x-small",
           bbox_to_anchor=(0.8, 0.1, 0.15, 0.15), frameon=False)
