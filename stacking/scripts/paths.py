# copy_spectra.py
#
# Created by Brett H. Andrews on 12 Jun 2017.

"""
1. Build file paths to spectra on passport
2. Access spectra and load into array
3. Bootstrap resample from spectra array
4. Stack spectra
5. Put each stack into array of stacks
6. Calculate mean and stddev of array of stacks
"""

import os
from os.path import join

import click
import pandas as pd


@click.command()
@click.option('--stackname', default='dr7_M0.1e')
@click.option('--binnames', default=None, multiple=True)
@click.option('--path_mzr', default=None)
def generate_filepaths(stackname, binnames, path_mzr):
    """Build file paths to spectra on passport external hard drive.
    
    Parameters:
        stackname (str):
            Name of the set of stacks. Default is ``dr7_M0.1e``.
        binnames (str):
            Name of bins to copy. Default is ``None``, which copies
            all bins.
        path_mzr (str):
            Path to the parent directory containing the stack data,
            scripts, etc. Default is ``None``.            
    """
    # data directory on passport external hard drive
    path_data = join('//', 'Volumes', 'My Passport', 'osu', 'andrews', 'projects',
                     'mass-metallicity', 'data', 'raw_FITS_dr7')
    assert os.path.isdir(path_data), f'``path_data`` does not exist: {path_data}'

    # mass-metallicity relation project directory on MacBook Air
    path_mzr_default = join(os.path.expanduser('~'), 'projects', 'mzr')
    path_mzr = path_mzr if path_mzr is not None else path_mzr_default
    assert os.path.isdir(path_mzr), f'``path_mzr`` does not exist: {path_mzr}'

    # stacks directory on MacBok Air
    path_stack = join(path_mzr, 'stacks', stackname)
    assert os.path.isdir(path_stack), f'``path_stack`` does not exist: {path_stack}'

    # read in bin masses (and SFRs) and number of galaxies per bin
    if not binnames:
        with open(join(path_stack, 'auxiliary', 'binnames.txt'), 'r') as fin:
            binnames = [line.strip() for line in fin]

    
    for binname in binnames:
        binpar = binname.split('_n')[0]
        path_filelist = join(path_stack, binpar, 'filelist', binname + '_filenames.txt')

        # paths to files on pulsar (OSU desktop)
        files_pulsar = pd.read_csv(path_filelist, delim_whitespace=True, header=None, usecols=(8,))
        files_pulsar.columns = ['paths']
        
        # paths to files on passport external hard drive
        files_passport = [join(path_data, *it.split('/')[-2:]) for it in files_pulsar.paths]
        
        # write paths to files on passport to a local file
        path_filelists = join(path_stack, filelists)
        if not os.path.isdir(path_filelists):
            os.mkdir(path_filelists)

        filelist_out = join(path_filelist, binpar, '.txt')
        with open(filelist_out, 'w') as fout:
            for it in files_passport:
                fout.write(it + '\n')


    

