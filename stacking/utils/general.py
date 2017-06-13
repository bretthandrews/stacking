# general.py
#
# deredden() adapted from ccm_deredden() in lineflux_v1a.py
#
# Created by Brett H. Andrews on 13 Jun 2017.

import numpy as np

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

