"""
=====================================
Plotting genomic data in inner circle
=====================================

"""
import numpy as np
import matplotlib.pyplot as plt

from circhic import datasets
from circhic.tools import plotCirc

counts, lengths = datasets.load_ccrescentus()
cumul_raws = counts.sum(axis=0)

ori_ccres = 1  # FIXME check the origin of *C crescentus*

r_in = 0.2
s_in, s_out = 1200000, 700000
genomic_data = {'data': -cumul_raws, 'r_min': 0.2, 'r_max': 0.5}
marks = {'ori': {'bin': 0, 'marker': 'P', 'ms': 18, 'color': 'r'},
         'ter': {'bin': int(lengths.sum()/2), 'marker': 'X',
                 'ms': 18, 'color': 'm'}}

plotCirc(counts, r_in, s_in, s_out, genomic_data, marks, pos0=ori_ccres)
plt.title(r'\textrm{with genomic data inside (unmapability)}', pad=-10);
