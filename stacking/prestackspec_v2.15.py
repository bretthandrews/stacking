#!/usr/local/python-2.7.1/bin/python2.7 -W ignore::DeprecationWarning

"""
FILE
    prestackspec_v2.15.py

DESCRIPTION

    Remove duplicates.

    Remove AGN using BPT diagram and emission line cuts as described in Tremonti
    et al. (2004) and using the separation criteria from Kauffmann et
    al. (2003).  These galaxy ID numbers are in starforming.pck.

    Select galaxies using the database file with Plate Id, MJD, Fiber ID,
    redshift, Mstar, SFR, 12+log(O/H), and filepath.

    Select galaxies that are flagged as galaxies and were NOT flagged as
    DEBLEND_NOPEAK or DEBLENDED_AT_EDGE.

    Then stack the galaxies.  It also outputs the results of the selection
    procedure (i.e., it writes the binlist.txt and binpar.txt files).

    Galaxies are individually corrected for foreground Milky Way reddening
    according to the E(B-V) value from Schlegel, Finkbeiner, & Davis (1998).

    Grid:
    starting wavelength w1 = 3700
    ending wavelength w2 = 7360
    bin size dw = 1.0
    length of spectrum nw = 3661

    Similar to prestackspec_v2.0.py but bins in log(Mstar).  Also limits the
    redshift range to 0.027 <= z <= 0.25 to ensure that [OII]7330 is in the
    spectrum.

BHA 21/Jan/2013
"""

import os, pyfits, pickle
import numpy as np
from scipy import interpolate
import stack_algorithm as sa
import matplotlib.pyplot as plt
import lineflux_v1a as lf

##### EDIT #####
stackname = 'dr7_M0.1e'
###############

stem_data = '/home/pulsar/andrews/projects/mass-metallicity/data/' + \
            'raw_FITS_extras_dr7/'
stem_out = '/home/pulsar/andrews/projects/mass-metallicity/stacks/'
stdir = stem_out + stackname + '/'

##### Make Directories for Output Part I #####
os.mkdir(stdir)
os.mkdir(stdir + 'auxiliary')
os.mkdir(stdir + 'plots')
os.mkdir(stdir + 'results')
os.mkdir(stdir + 'results/metallicity')
os.mkdir(stdir + 'results/temperature')
os.mkdir(stdir + 'results/ionic_abundance')
os.mkdir(stdir + 'results/density')
os.mkdir(stdir + 'results/line_fluxes')
##############################################

##### Read in Data from MPA/JHU Catalog #####
fin = open(stem_data + 'master_data_dr7b.pck', 'r')
pid, mjd, fid, z, z_err, sn_median, ebv, mstar, sfr, ssfr, oh, clss = \
     pickle.load(fin)
fin.close()
fnames = np.loadtxt(stem_data + 'filenames_dr7.txt', dtype='S105')

# Starforming galaxies
fSF = open(stem_data + 'starforming.pck', 'r')
indSF = pickle.load(fSF)
fSF.close()
#############################################


##### Remove Duplicates #####
fdup = open(stem_data + 'duplicates_dr7.txt', 'r')
orig = [] # original index
dup = [] # index of duplicates
for line in fdup:
    if line.strip().split()[0] == '#':
        print ''
    else:
        index, duplicates = line.strip().split()[0], line.strip().split()[1:]
        orig.append(index)
        dup.append(duplicates)    

fdup.close()
# if duplicate is greater than index (effectively keeps the first spectrum that
# is encountered), replace duplicate index with -1
orig2 = np.copy(orig).astype(int)
n_ind = len(orig)
for i in xrange(n_ind):
    for j in xrange(len(dup[i])):
        dint = int(dup[i][j])
        if (dint >= 0) & (dint > i):
            orig2[dint] = -1

# orig3 contains the indices of the non-duplicate galaxies (868492 non-duplicate
# galaxies)
orig3 = orig2[np.where(orig2 >= 0)]
#############################



##### Read in Non-Flagged Galaxies #####
'''Keep galaxies NOT flagged as DEBLEND_NOPEAK and DEBLENDED_AT_EDGE.'''
ffl = open(stem_data + '/ind_fl_arr.pck', 'r')
ind_fl_arr = pickle.load(ffl)
ffl.close()
###################################



##### Removal #####
'''Cross-match the indices of the star-forming galaxies, the non-duplicate
galaxies, and the non-flagged galaxies.'''
ind_sforig = np.intersect1d(indSF, orig3)
ind_keep = np.intersect1d(ind_sforig, ind_fl_arr)

###################





##### Make Bins #####
bins = np.arange(7.0, 11.51, 0.1)
binnames = ['M' + str(bins[i]) + '_' + str(bins[i+1]) for i in xrange(len(bins) - 1)]
#####################

##### Make Directories for Output Part II #####
for i in xrange(len(bins)-1):
    os.mkdir(stdir + binnames[i])
    os.mkdir(stdir + binnames[i] + '/filelist')
    os.mkdir(stdir + binnames[i] + '/final_spec')
    os.mkdir(stdir + binnames[i] + '/grid_files')
    os.mkdir(stdir + binnames[i] + '/raw_stack')
    os.mkdir(stdir + binnames[i] + '/STARLIGHT')
    os.mkdir(stdir + binnames[i] + '/final_spec/sfdb')

##############################################



##### Look at Individual Galaxies #####

'''For galaxies between 10^7 and 10^8.6 Msun, inspect galaxies by eye and throw
out those that are spurious. ind_bad are the spurious or unusable galaxies, and
ind_good are the usable galaxies in this mass range.'''

mrange = (mstar > 7.) & (mstar < 8.6)
zrange = (z >= 0.027) & (z <= 0.25)
accurate_z = (np.abs(z_err) < 0.001)
criterion = mrange & zrange & accurate_z
ind11 = np.where(criterion)[0]
ind_view = np.intersect1d(ind11, ind_keep)

# ind_bad1 are the indices of the bad galaxies within the ind_view array
ind_bad1 = np.array([1, 9, 10, 16, 19, 28, 41, 54, 59, 63, 71, 74, 101, 115,
                     132, 136, 151, 163, 166, 167, 172, 189, 190, 201, 230, 242,
                     265, 275, 294, 295, 302, 308, 311, 313, 314, 321, 335, 336,
                     339, 341, 356, 361, 362, 377, 380, 382, 396, 414, 424, 430,
                     441, 446, 451, 453, 454, 466, 476, 478, 480, 493, 496, 501,
                     518, 551, 554, 559, 565, 575, 578, 580, 591, 597, 600, 608,
                     612, 621, 637, 639, 654, 663, 664, 669, 681, 690, 693, 705,
                     709, 714, 715, 720, 722, 746, 751, 758, 778, 795, 808])


# ind_bad and ind_good are the indices of the bad/good galaxies within the
# overall galaxy numbering system
ind_bad = ind_view[ind_bad1]
ind_good = np.delete(ind_view, np.s_[ind_bad1])

# ind_kill are the indices of the bad galaxies in the ind_keep array
ind_kill = []
for i in xrange(len(ind_bad)):
    ind_kill.append(np.where(ind_keep == ind_bad[i])[0][0])

ind_kill = np.array(ind_kill)

# ind_keep2 are the indices of ALL of the good galaxies (i.e., the good galaxies
# selected by hand for 7.0 < logM < 8.6 and all of the galaxies that met the
# other requirements at logM > 8.6
ind_keep2 = np.delete(ind_keep, np.s_[ind_kill])



############################################




##### Stack Galaxies #####
'''Stack the galaxies.  Python returns a MemoryError for n_spectra > 18400.'''

names = []
median_mstar = []
median_sfr = []
for m in xrange(len(bins) - 1):
    # Select galaxies
    mrange = (mstar >= bins[m]) & (mstar < bins[m+1])
    zrange = (z >= 0.027) & (z <= 0.25)
    accurate_z = (np.abs(z_err) < 0.001)
    criterion = mrange & zrange & accurate_z
    ind1 = np.where(criterion)[0]
    ind = np.intersect1d(ind1, ind_keep2)
    filenames = fnames[ind]
    ind_missing = []
    for i in xrange(len(filenames)):
        if not os.path.isfile(filenames[i]):
            ind_missing.append(i)
    filenames = np.delete(filenames, np.s_[ind_missing])
    ind = np.delete(ind, np.s_[ind_missing])
    n_spectra = len(ind)
    # stack
    if n_spectra > 0:
        median_mstar.append(np.median(mstar[ind]))
        median_sfr.append(np.median(sfr[ind]))
        grid = np.arange(3700, 7360.1, 1.0, float)
        spectra = np.zeros((n_spectra, len(grid))) # 2D array with each row a spectrum
        for i in xrange(n_spectra):
            hdu = pyfits.open(filenames[i]) # open fits file
            primary_data = hdu[0].data # load data
            COEFF0 = hdu[0].header['COEFF0'] # central wavelength of first pixel
            COEFF1 = hdu[0].header['COEFF1'] # Angstroms per pixel
            j = np.arange(0, len(primary_data[0]), 1, int) # array of bin numbers
            vac = 10.**(COEFF0 + COEFF1 * j) # vacuum wavelengths
            air = vac / (1.0 + 2.735182E-4 + 131.4182 / vac**2 + \
                         2.76249E8 / vac**4) # air wavelengths
            rest_wave = air / (1 + z[ind[i]]) # rest wavelengths
            dered_spec = lf.ccm_deredden(rest_wave, primary_data[0], ebv[ind[i]]) # corrected for MW reddening
            raw_spec = dered_spec * (1 + z[ind[i]]) # shifted to rest frame
            interp_func = interpolate.interp1d(rest_wave, raw_spec, \
                                               bounds_error=False, fill_value=0.)
            regrid_spec = interp_func(grid) # linear interpolation
            w1 = np.where(grid == 4400) # grab array element number for 4400 Angstroms
            w2 = np.where(grid == 4450) # grab array element number for 4450 Angstroms
            mean_cont_flux = np.mean(regrid_spec[w1[0][0]:w2[0][0]+1]) # mean continuum flux from w1-w2
            regrid_spec = regrid_spec / mean_cont_flux # normalize spectrum
            spectra[i] = regrid_spec
            hdu.close()
        # write final filelist
        stname = binnames[m] + '_n' + str(n_spectra)
        names.append(stname)
        f1 = open(stdir + binnames[m] + '/filelist/' + stname + '_filenames.txt', 'w')
        f1_write = f1.write
        for j in xrange(n_spectra):
            f1_write('%-6s %-8s %-5s %-8s %-10s %-12s %-8s %-6s %s \n' % \
                     (pid[ind[j]], mjd[ind[j]], fid[ind[j]], z[ind[j]], \
                      mstar[ind[j]], sfr[ind[j]], oh[ind[j]], clss[ind[j]], \
                      fnames[ind[j]]))
        f1.close()
        stack_spec = sa.stack(grid, spectra, method='mean')
        ASCIIoutput = stdir + binnames[m] + '/raw_stack/' + stname + '.txt'
        np.savetxt(ASCIIoutput, np.transpose((grid, stack_spec)), fmt='%5.1f     %12.8f')





##### Print Bin Names and Bin Parameters to a File #####
bnms = open(stdir + 'auxiliary/' + 'binnames.txt', 'w')
bpar = open(stdir + 'auxiliary/' + 'binpar.txt', 'w')
bpar.write('#%-8s%-8s%-8s%-16s%-16s\n' %
               ('m_low', 'm_up', 'n_gal', 'median_mstar', 'median_sfr'))
for p in xrange(len(names)):
    m_low = names[p].split('_')[0].split('M')[1]
    m_up = names[p].split('_')[1]
    num = names[p].split('_')[2].split('n')[1]
#    bnms.write(names[p] + '\n')
    bpar.write('%-8s%-8s%-8s%-16.4f%-16.4f\n' % 
               (m_low, m_up, num, median_mstar[p], median_sfr[p]))

bnms.close()
bpar.close()
########################################################
