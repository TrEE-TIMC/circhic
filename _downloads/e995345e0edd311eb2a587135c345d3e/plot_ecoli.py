"""
===============
*E. coli*
===============

Loading and plotting *E. coli* data.
"""

from circhic import datasets
import matplotlib.pyplot as plt
from matplotlib import colors

data = datasets.load_ecoli()
counts = data["counts"]

fig, ax = plt.subplots()
ax.imshow(counts, norm=colors.SymLogNorm(1), interpolation="none")
