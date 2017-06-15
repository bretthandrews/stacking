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

    path_filelists = join(path_dr7, 'filelists')
    if not filelists:
        filelists = os.listdir(path_filelists)
        filelists.sort(key=lambda s: float(s.split('M')[1].split('_')[0]))

    # filelists = ['M9.4_9.5.txt', 'M9.5_9.6.txt', 'M9.6_9.7.txt']

    path_snr = join(path_dr7, 'results', 'snr', 'snr.csv')

    if os.path.isfile(path_snr) and not overwrite:
        raise ValueError(f'Not written (overwrite with --overwrite): {path_snr}')

    full_snr = []
    window_snr = []
    indices = []

    for filelist in filelists:
        click.echo(filelist)

        binpar = filelist.split('.txt')[0]
        path_raw = join(path_dr7, binpar, 'raw_stack')
        path_spec = join(path_raw, binpar + '.pck')

        with open(path_spec, 'rb') as fin:
            spectra = pickle.load(fin)

        stack = np.mean(spectra, axis=0)

        check_stack(stack, binpar, path_raw, path_dr7)
    
        grid = make_grid()
        
        indices.append(binpar)

        spec_mean, spec_std = resample_stacks(spectra, samples)
        spec_snr = stack / spec_std
        spec_median_snr = np.median(spec_snr)
        window_median_snr = np.median(spec_snr[(grid >=4400) & (grid <= 4450)])
        full_snr.append(spec_median_snr)
        window_snr.append(window_median_snr)

    snr = pd.DataFrame(list(map(list, zip(*(full_snr, window_snr)))), index=indices,
                       columns=['full spec', '4400-4450 A'])
    snr.to_csv(path_snr)


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


snr = pd.read_csv('/Users/andrews/projects/mzr/stacks/dr7_M0.1e/results/snr/snr.csv')
snr.plot()
