# stack.py
#
# Adapted from mzr/scripts/dr7_M0.1e/prestackspec_v2.15.py
#
# Created by Brett H. Andrews on 13 Jun 2017.


import os
from os.path import join

import numpy as np
from scipy import interpolate
import pandas as pd
import sklearn

from stacking.utils.general import deredden, load_spectrum, mean_flux, get_table_index, rest_wave

@click.command()
@click.option('--filelists', default=None, multiple=True)
@click.option('--overwrite', '-o', default=False, is_flag=True)
def bootstrap_stack(filelists, overwrite):
    """
    
    Parameters:
        filelists (str):
            Filelists of spectra to stack. Default is ``None``, which
            will use all stacks.
        overwrite (bool):
            If ``True``, overwrite existing files. Default is ``False``.         
    """

path_mzr = join(os.path.expanduser('~'), 'projects', 'mzr')

# Read in redshift and E(B-V) values
path_meta = join(path_mzr, 'data', 'raw_FITS_extras_dr7')
table = pd.read_csv(join(path_meta, 'master_data_dr7b.csv'))

path_filelists = join(path_mzr, 'stacks', 'dr7_M0.1e', 'filelists')
filelists = None  # Remove............................................................................
filelists = os.listdir(path_filelists) if filelists is None else filelists


grid = np.arange(3700, 7360.1, 1.0, float)

# for filelist in filelists:

filelist = 'M8.2_8.3.txt'  # Remove....................................................................
path_filelist = join(path_filelists, filelist)

with open(path_filelist, 'r') as fin:
    filenames = [line.strip() for line in fin]

# def load_spectra(filenames, redshifts, ebvs):
spectra = np.zeros((len(filenames), len(grid)))

for ii, filename in enumerate(filenames):

    ind = get_table_index(table, filename)
    
    spec_obs, wave_cen_pix1, ang_per_pix = load_spectrum(filename)
    wave_rest = rest_wave(wave_cen_pix1, ang_per_pix, len(spec_obs), table.z[ind])
    spec_dered = deredden(wave_rest, spec_obs, table.ebv[ind])
    spec_raw = spec_dered * (1 + table.z[ind])

    linear_interp = interpolate.interp1d(wave_rest, spec_raw, bounds_error=False, fill_value=0.)
    spec_regrid = linear_interp(grid)
    
    mean_cont_flux = mean_flux(spec_regrid, grid, wave_low=4400, wave_upp=4500)
    spec_normalized = spec_regrid / mean_cont_flux
    
    spectra[ii] = spec_normalized


n_samples = 50  # Remove................................................................................
stacks = np.zeros((n_samples, len(grid)))

for ii in range(n_samples):
    spectra_resampled = sklearn.utils.resample(spectra)
    stacks[ii] = np.mean(spectra_resampled, axis=0)
    
spec_mean = np.mean(stacks, axis=0)
spec_std = np.std(stacks, ddof=1, axis=0)

# Plot
import matplotlib.pyplot as plt
import seaborn as sns

fig, ax = plt.subplots()
ax.plot(grid, spec_mean)
ax.plot(grid, spec_mean + spec_std)
ax.plot(grid, spec_mean - spec_std)

