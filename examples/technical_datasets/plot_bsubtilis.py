"""
===============
*B. subtilis*
===============

Loading and plotting *B. subtilis* data.
"""

from circhic import datasets
import matplotlib.pyplot as plt
from matplotlib import colors

data = datasets.load_bsubtilis()
counts = data["counts"]

fig, ax = plt.subplots()
ax.imshow(counts, norm=colors.SymLogNorm(1), interpolation="none")
