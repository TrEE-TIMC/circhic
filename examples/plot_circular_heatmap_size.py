"""
=====================================
Changing the size of the inner circle
=====================================

"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import colors

from circhic import datasets
from circhic.tools import genCircData


counts, lengths = datasets.load_ccrescentus()

ori_ccres = 1  # FIXME check the origin of *C crescentus*
vmax = np.max([counts[i, (i+1) % lengths.sum()] for i in range(lengths.sum())])

vmin = vmax / 50
min_non_zero = np.min(counts[counts != 0])

plt.figure(figsize=(16, 10))

plt.subplot(1, 2, 1)

Circ = genCircData(counts, bin_circ=0.5, r_in=0.2, res=10000,
                   pos0=ori_ccres, s_in=1200000, s_out=600000)
plt.imshow(Circ, norm=colors.SymLogNorm(min_non_zero), vmin=vmin, vmax=vmax,
           cmap='viridis')
plt.xticks([])
plt.yticks([])
plt.axis('off')


plt.subplot(1, 2, 2)
Circ = genCircData(counts, bin_circ=0.5, r_in=0.4, res=10000,
                   pos0=ori_ccres, s_in=600000, s_out=600000)
plt.imshow(Circ, norm=colors.SymLogNorm(min_non_zero), vmin=vmin, vmax=vmax,
           cmap='viridis')
plt.xticks([])
plt.yticks([])
plt.axis('off')

plt.tight_layout()
