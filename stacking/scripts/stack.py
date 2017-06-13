# stack.py
#
# Adapted from mzr/scripts/dr7_M0.1e/prestackspec_v2.15.py
#
# Created by Brett H. Andrews on 13 Jun 2017.


import os
from os.path import join
from functools import reduce

import numpy as np
from scipy import interpolate
import pandas as pd
import sklearn
from astropy.io import fits

from stacking.utils.general import deredden

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

for ii, filename in enumerate(filenames):
    
    ind = get_table_index(table, filename)
    spec_obs, wave_cen_pix1, ang_per_pix = load_spectrum(filename)
    wave_rest = rest_wavelenghts(wave_cen_pix1, ang_per_pix, len(spec_obs))


# correct for MW reddening and shift to rest frame
spec_dered = deredden(wave_rest, spec_obs, table.ebv[ind])
spec_raw = spec_dered * (1 + table.z[ind])

linearly_interpolate = interpolate.interp1d(wave_rest, spec_raw, bounds_error=False, fill_value=0.)
spec_regrid = linearly_interpolate(grid)

# compute mean continuum flux from w1 to w2 (4400--4500 Angstroms)
w1 = np.where(grid == 4400)
w2 = np.where(grid == 4450)
mean_cont_flux = np.mean(spec_regrid[w1[0][0]:w2[0][0]+1])

spec_normalized = spec_regrid / mean_cont_flux
spectra[ii] = spec_normalized


stacks = np.zeros((n_samples, len(grid)))
for ii in range(n_samples):
    spectra_resampled = sklearn.utils.resample(spectra)
    stack_spec = np.mean(spectra_resampled, axis=1)



def get_table_index(table, filename):
    """Get index of table corresponding to spectrum.
    
    Parameters:
        table (DataFrame):
            Table of galaxy properties.
        filename (str):
            Name of spectrum FITS file.

    Returns:
        int
    """
    mjd, pid, fid = [int(it) for it in filename.split('spSpec-')[-1].strip('.fit').split('-')]
    ind_mjd = np.where(table.mjd == mjd)
    ind_pid = np.where(table.pid == pid)
    ind_fid = np.where(table.fid == fid)
    return reduce(np.intersect1d, (ind_mjd, ind_pid, ind_fid))[0]


def load_spectrum(filename):
    """Load spectrum from FITS file.
    
    Parameters:
        filename (str):
            Name of spectrum FITS file.

    Returns:
        array, float, float
    """
    fitsobj = fits.open(filename)
    spec_obs = fitsobj[0].data[0]
    wave_cen_pix1 = fitsobj[0].header['COEFF0']  # central wavelength of first pixel
    ang_per_pix = fitsobj[0].header['COEFF1']    # Angstroms per pixel
    fitsobj.close()
    return spec_obs, wave_cen_pix1, ang_per_pix


def rest_wavelengths(wave_cen_pix1, ang_per_pix, spec_len, redshift):
    """Create rest wavelength array.

    Parameters:
        wave_cen_pix1 (float):
            Central wavelength of first pixel.
        ang_per_pix (float):
            Angstroms per pixel.
        spec_len (int):
            Length of spectrum.
        redshift (float):
            Redshift of galaxy.
    Returns:
        array
    """
    vac = 10.**(wave_cen_pix1 + ang_per_pix * np.arange(spec_len))
    air = vac / (1.0 + 2.735182E-4 + 131.4182 / vac**2 + 2.76249E8 / vac**4)
    return air / (1 + redshift)
