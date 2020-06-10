"""
============================
Human (KBM7) - chromosome 14
============================

Loading and plotting KBM7's chromosome 14 data.
"""

from circhic import datasets
import matplotlib.pyplot as plt
from matplotlib import colors

data = datasets.load_kbm7()
counts = data["counts"]

fig, ax = plt.subplots()
ax.imshow(counts, norm=colors.SymLogNorm(1), interpolation="none")
