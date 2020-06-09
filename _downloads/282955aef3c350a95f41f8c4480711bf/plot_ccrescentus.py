"""
===============
*C. crescentus*
===============

Loading and plotting *C. crescentus* data.
"""

from circhic import datasets
import matplotlib.pyplot as plt
from matplotlib import colors

counts, lengths = datasets.load_ccrescentus()

fig, ax = plt.subplots()
ax.imshow(counts, norm=colors.SymLogNorm(1), interpolation="none")
