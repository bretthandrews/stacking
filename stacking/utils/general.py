# general.py
#
# deredden() adapted from ccm_deredden() in lineflux_v1a.py
#
# Created by Brett H. Andrews on 13 Jun 2017.


import numpy as np
from astropy.io import fits


def deredden(wave, flux_obs, ebv, R_V=3.1):
    """Deredden using Cardelli, Clayton, & Mathis (1989) extinction law."""
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


def mean_flux(spectrum, wave_low, wave_upp):
    """Compute mean flux in a wavelength window.
    
    Parameters:
        spectrum (array):
            Spectrum.
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
