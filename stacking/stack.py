# stack.py
#
# Created by Brett H. Andrews on 13 Jun 2017.


import os
from os.path import join
import pickle
import re

import click
import numpy as np
import pandas as pd
import sklearn
import matplotlib.pyplot as plt
import seaborn as sns

from stacking.spectral_processing import make_grid


@click.command()
@click.option('--filelists', default=None, multiple=True)
@click.option('--overwrite', '-o', default=False, is_flag=True)
@click.option('--samples', default=100)
def stack_resample(filelists, overwrite, samples):
    """Stack spectra using bootstrap resampling.

    Parameters:
        filelists (str):
            Filelists of spectra to stack. Default is ``None``, which
            will use all stacks.
        overwrite (bool):
            If ``True``, overwrite existing files. Default is ``False``.
        samples (int):
            Number of samples (with replacement) to draw. Default is
            ``100``.
    """
    path_mzr = join(os.path.expanduser('~'), 'projects', 'mzr')
    path_dr7 = join(path_mzr, 'stacks', 'dr7_M0.1e')

    # path_filelists = join(path_dr7, 'filelists')
    # if not filelists:
    #     filelists = os.listdir(path_filelists)
    #     filelists.sort(key=lambda s: float(s.split('M')[1].split('_')[0]))
    filelists = make_filelists()

    for filelist in filelists:
        click.echo(filelist)

        binpar = filelist.split('.txt')[0]
        path_raw = join(path_dr7, binpar, 'raw_stack')
        path_spec = join(path_raw, binpar + '.pck')

        with open(path_spec, 'rb') as fin:
            spectra = pickle.load(fin)

        stack = np.mean(spectra, axis=0)

        # TODO check for convergence
        # draw 10, 20, 50, 100 times and see how std changes
        # easiest to do over a narrow range (collapse down to one number)
        # stds = np.array(100)

        spec_mean, spec_std = resample_stacks(spectra, samples)

        check_stack(stack, binpar, path_raw, path_dr7)


def resample_stacks(spectra, samples):
    """Resample stacks.

    Parameters:
        spectra (array):
            Spectra to be sampled from.
        samples (int):
            Number of samples to draw.

    Returns:
        array, array: mean and std spectra
    """
    grid = make_grid()
    stacks_resampled = np.zeros((samples, len(grid)))

    for ii in range(samples):
        spectra_resampled = sklearn.utils.resample(spectra)
        stacks_resampled[ii] = np.mean(spectra_resampled, axis=0)

    spec_mean = np.mean(stacks_resampled, axis=0)
    spec_std = np.std(stacks_resampled, ddof=1, axis=0)

    return spec_mean, spec_std


def check_stack(stack, binpar, path_raw, path_dr7):
    """Check that current stack agrees with stack from AM13.

    Parameters:
        stack (array):
            Current stack.
        binpar (str):
            Bin parameters.
        path_raw (str):
            Path to raw_stack directory.
        path_dr7 (str):
            Path to dr7_M0.1e directory.
    """
    for filename in os.listdir(path_raw):
        fm = re.fullmatch(f'{binpar}_n\d+.txt', filename)
        if fm is not None:
            path_stack_comp = join(path_dr7, binpar, 'raw_stack', fm.string)
            break

    comp = pd.read_csv(path_stack_comp, delim_whitespace=True, header=None, names=['wave', 'flux'])
    stack_comp = comp.flux.values

    # print(filelist, str(np.abs(stack - stack_comp).max()))
    assert np.abs(stack - stack_comp).max() < 1e-6, \
        'stack and comparison stack differ by more than 1e-6'



def plot_spec(grid, spec, spec_std):
    """spec with std"""
    c0 = sns.color_palette()[0]
    fig, ax = plt.subplots()
    ax.plot(grid, spec, c=c0)
    ax.fill_between(grid, spec - spec_std, spec + spec_std, facecolor=c0, alpha=0.5)


def make_filelists():
    """stacks"""
    return [
        'M7.0_7.1.txt',
        'M7.1_7.2.txt',
        'M7.2_7.3.txt',
        'M7.3_7.4.txt',
        'M7.4_7.5.txt',
        'M7.5_7.6.txt',
        'M7.6_7.7.txt',
        'M7.7_7.8.txt',
        'M7.8_7.9.txt',
        'M7.9_8.0.txt',
        'M8.0_8.1.txt',
        'M8.1_8.2.txt',
        'M8.2_8.3.txt',
        'M8.3_8.4.txt',
        'M8.4_8.5.txt',
        'M8.6_8.7.txt',
        'M8.7_8.8.txt',
        'M8.8_8.9.txt',
        'M8.9_9.0.txt',
        'M9.0_9.1.txt',
        'M9.1_9.2.txt',
        'M9.2_9.3.txt',
        'M9.3_9.4.txt',
        'M9.4_9.5.txt',
        # 'M9.5_9.6.txt',
        # 'M9.6_9.7.txt',
        # 'M9.7_9.8.txt',
        # 'M9.8_9.9.txt',
        # 'M9.9_10.0.txt',
        # 'M10.0_10.1.txt',
        # 'M10.1_10.2.txt',
        # 'M10.2_10.3.txt',
        # 'M10.3_10.4.txt',
        # 'M10.4_10.5.txt',
        # 'M10.5_10.6.txt',
        # 'M10.6_10.7.txt',
        # 'M10.7_10.8.txt',
        # 'M10.8_10.9.txt',
        # 'M10.9_11.0.txt',
        # 'M11.0_11.1.txt',
        # 'M11.1_11.2.txt',
        # 'M11.2_11.3.txt',
        # 'M11.3_11.4.txt',
        # 'M11.4_11.5.txt'
        ]
