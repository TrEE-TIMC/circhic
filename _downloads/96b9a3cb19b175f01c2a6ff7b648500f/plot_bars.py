"""
==============
Plotting bars
==============

Here's a small example showcasing how to plot bars on a circular strip.
"""
import numpy as np
from circhic._base import CircHiCFigure

###############################################################################
# First, simulate some data
lengths = np.array([3500])
random_state = np.random.RandomState(42)
data = random_state.randn(100)

###############################################################################
# Then, create the `circhic` figure and plot the bars
circhicfig = CircHiCFigure(lengths)
_, ax = circhicfig.plot_bars(
    data,
    inner_radius=0.5,
    color="0")
ax.set_title("Plotting bars", fontweight="bold")
