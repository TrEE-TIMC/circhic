"""
===============================================
a circHiC figure showing biases and mappability
===============================================

"""
import numpy as np
from iced.normalization import ICE_normalization

from circhic import datasets
from circhic._base import CircHiCFigure


counts, lengths = datasets.load_ccrescentus()

cumul_raw_counts = counts.sum(axis=0)

counts, bias = ICE_normalization(counts, output_bias=True)

# Now replace missing data with NA
missing_loci = counts.sum(axis=0) == 0
counts[missing_loci] = np.nan
counts[:, missing_loci] = np.nan


circhicfig = CircHiCFigure(lengths)
circhicfig.plot_hic(counts, outer_radius=0.75)

# Assume you want to plot data from that ranges in a polar plot outside of the
# first one. Then the 0 axis should be at, say, 80% of the axis
circhicfig.plot_bars(cumul_raw_counts, inner_radius=0.8, outer_radius=0.9,
                     color="0")

# Now, plot a last plot, for the top 10% of the original axes
circhicfig.plot_lines(bias, color="0", inner_radius=0.9, outer_radius=1)


# Now, try to do a simple "theta-axis" on the outer
ax = circhicfig.set_genomic_ticklabels()
