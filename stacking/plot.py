# plot.py
#
# Created by Brett H. Andrews on 15 Jun 2017.


import os
from os.path import join

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

from stacking.spectral_processing import make_grid


def plot_spec(spec, spec_std):
    """spec with std"""
    grid = make_grid()
    
    c0 = sns.color_palette()[0]
    fig, ax = plt.subplots()
    ax.plot(grid, spec, c=c0)
    ax.fill_between(grid, spec - spec_std, spec + spec_std, facecolor=c0, alpha=0.5)


def plot_snr():
    path_snr = join(os.path.expanduser('~'), 'projects', 'mzr', 'stacks', 'dr7_M0.1e', 'results',
                    'snr', 'snr.csv')
    snr = pd.read_csv(path_snr, index_col=0)
    snr.plot()
