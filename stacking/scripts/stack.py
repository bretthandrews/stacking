# stack.py
#
# Adapted from mzr/scripts/dr7_M0.1e/prestackspec_v2.15.py
#
# Created by Brett H. Andrews on 13 Jun 2017.


import os
from os.path import join
from functools import reduce

import numpy as np
import pandas as pd
import sklearn
from astropy.io import fits

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
filelists = os.listdir(path_filelists) if filelists is None else filelists


grid = np.arange(3700, 7360.1, 1.0, float)

# for filelist in filelists:

filelist = 'M8.2_8.3.txt'  # Remove....................................................................
path_filelist = join(path_filelists, filelist)

with open(path_filelist, 'r') as fin:
    filenames = [line.strip() for line in fin]

# def load_spectra(filenames, redshifts, ebvs):
spectra = np.zeros((len(filenames), len(grid)))
    
# for ii, filename in enumerate(filenames):
filename = filenames[0]  # Remove......................................................................


# get index of table corresponding to spectrum
mjd, pid, fid = [int(it) for it in filename.split('spSpec-')[-1].strip('.fit').split('-')]
ind_mjd = np.where(table.mjd == mjd)
ind_pid = np.where(table.pid == pid)
ind_fid = np.where(table.fid == fid)
ind = reduce(np.intersect1d, (ind_mjd, ind_pid, ind_fid))[0]

# load spectrum
fitsobj = fits.open(filename)
spec_obs = fitsobj[0].data[0]
wave_cen_pix1 = fitsobj[0].header['COEFF0']  # central wavelength of first pixel
ang_per_pix = fitsobj[0].header['COEFF1']    # Angstroms per pixel
fitsobj.close()

# transform from vacuum to air to rest wavelengths
vac = 10.**(wave_cen_pix1 + ang_per_pix * np.arange(len(spec_obs))
air = vac / (1.0 + 2.735182E-4 + 131.4182 / vac**2 + 2.76249E8 / vac**4)
wave_rest = air / (1 + z[ind])

spec_dered = lf.ccm_deredden(wave_rest, data[0], ebv[ind])  # correct for MW reddening
spec_raw = spec_dered * (1 + z[ind])                        # shift to rest frame

# linearly interpolate spectrum
interp_func = interpolate.interp1d(wave_rest, spec_raw, bounds_error=False, fill_value=0.)
spec_regrid = interp_func(grid)

# compute mean continuum flux from w1 to w2 (4400--4500 Angstroms)
w1 = np.where(grid == 4400)
w2 = np.where(grid == 4450)
mean_cont_flux = np.mean(spec_regrid[w1[0][0]:w2[0][0]+1])

spec_normalized = spec_regrid / mean_cont_flux
spectra[ii] = spec_normalized

sklearn.utils.resample(spectra, n_samples=n_samples, random_state=1234)


