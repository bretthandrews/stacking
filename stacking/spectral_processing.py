# spectral_processing.py
#
# Adapted from mzr/scripts/dr7_M0.1e/prestackspec_v2.15.py
#
# deredden() adapted from ccm_deredden() in lineflux_v1a.py
#
# Created by Brett H. Andrews on 13 Jun 2017.


import os
from os.path import join
import pickle

import numpy as np
from scipy.interpolate import interp1d
import pandas as pd
from astropy.io import fits
import click

from stacking.paths import get_table_index


@click.command()
@click.option('--filelists', default=None, multiple=True)
@click.option('--overwrite', '-o', default=False, is_flag=True)
def process_spectra(filelists, overwrite):
    """Deredden, deredshift, regrid, and normalize spectra.
    
    Parameters:
        filelists (str):
            Filelists of spectra to stack. Default is ``None``, which
            will use all stacks.
        overwrite (bool):
            If ``True``, overwrite existing files. Default is ``False``.         
    """

    path_mzr = join(os.path.expanduser('~'), 'projects', 'mzr')
    path_dr7 = join(path_mzr, 'stacks', 'dr7_M0.1e')

    # Read in redshift and E(B-V) values
    path_meta = join(path_mzr, 'data', 'raw_FITS_extras_dr7')
    table = pd.read_csv(join(path_meta, 'master_data_dr7b.csv'))

    path_filelists = join(path_dr7, 'filelists')
    if not filelists:
        filelists = os.listdir(path_filelists)
        filelists.sort(key=lambda s: float(s.split('M')[1].split('_')[0]))

    for filelist in filelists:
        click.echo(filelist)
        path_filelist = join(path_filelists, filelist)

        with open(path_filelist, 'r') as fin:
            filenames = [line.strip() for line in fin]
        
        binpar = filelist.split('.txt')[0]
        path_spec_out = join(path_dr7, binpar, 'raw_stack', binpar + '.pck')

        if not os.path.isfile(path_spec_out) or overwrite:
            spectra = homogenize_spectra(filenames=filenames, table=table)

            with open(path_spec_out, 'wb') as fout:
                pickle.dump(spectra, fout)

            click.echo(f'{path_spec_out}')
        else:
            click.echo(f'Not written (overwrite with --overwrite): {path_spec_out}')


def make_grid(wave_low=3700., wave_upp=7360.1, dlam=1.):
    """Create wavelength grid.
    
    Parameters:
        wave_low (float):
            Lower limit of wavelength grid in Angstroms. Default is
            ``3700.``.
        wave_upp (float):
            Upper limit of wavelength grid in Angstroms. Default is
            ``7360.1``.
        dlam (float):
            Pixel size of wavelength grid in Angstroms. Default is
            ``1.``.
    """
    return np.arange(wave_low, wave_upp, dlam, float)


def homogenize_spectra(filenames, table):
    """Deredden, deredshift, regrid, and normalize spectra.
    
    Parameters:
        filenames (list):
            Names of spectra FITS files.
        table (DataFrame):
            Data table including redshift and ebv.

    Returns:
        array
    """
    grid = make_grid()
    spectra = np.zeros((len(filenames), len(grid)))

    for ii, filename in enumerate(filenames):

        ind = get_table_index(table, filename)

        spec_obs, wave_cen_pix1, ang_per_pix = load_spectrum(filename)
        wave_rest = rest_wave(wave_cen_pix1, ang_per_pix, len(spec_obs), table.z[ind])
        spec_dered = deredden(wave_rest, spec_obs, table.ebv[ind])
        spec_raw = spec_dered * (1 + table.z[ind])

        linear_interp = interp1d(wave_rest, spec_raw, bounds_error=False, fill_value=0.)
        spec_regrid = linear_interp(grid)

        mean_cont_flux = mean_flux(spec_regrid, grid, wave_low=4400, wave_upp=4450)
        spec_normalized = spec_regrid / mean_cont_flux

        spectra[ii] = spec_normalized
    
    return spectra


def deredden(wave, flux_obs, ebv, R_V=3.1):
    """Deredden using Cardelli, Clayton, & Mathis (1989) extinction law.
    
    Parameters:
        wave (array):
            Wavelength.
        flux_obs (array):
            Oberved flux.
        ebv (float):
            E(B-V).
        R_V (float):
            Default is ``3.1``.
    
    Returns:
        array
    """
    xx = 1. / (wave / 10000.)
    yy = xx - 1.82
    
    a = (1. + (0.17699 * yy) - (0.50447 * yy**2.) - (0.02427 * yy**3.) + (0.72085 * yy**4.) +
         (0.01979 * yy**5.) - (0.77530 * yy**6.) + (0.32999 * yy**7.))
    
    b = ((1.41338 * yy) + (2.28305 * yy**2.) + (1.07233*yy**3.) - (5.38434 * yy**4.) -
         (0.62251 * yy**5.) + (5.30260 * yy**6.) - (2.09002 * yy**7.))
    
    A_V = R_V * ebv
    A_lam = A_V * (a + (b / R_V))
    
    flux_dered = flux_obs * 10.**(0.4 * A_lam)
    
    return flux_dered


def mean_flux(spectrum, grid, wave_low, wave_upp):
    """Compute mean flux in a wavelength window.
    
    Parameters:
        spectrum (array):
            Spectrum.
        grid (array):
            Wavelength grid.
        wave_low (float):
            Blue limit of wavelength window in Angstroms.
        wave_upp (float):
            Red limit of wavelength window in Angstroms.
    
    Returns:
        array
    """
    ind1 = np.where(grid == wave_low)[0][0]
    ind2 = np.where(grid == wave_upp)[0][0]
    return np.mean(spectrum[ind1:ind2+1])



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


def rest_wave(wave_cen_pix1, ang_per_pix, spec_len, redshift):
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
