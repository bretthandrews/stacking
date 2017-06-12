# copy_spectra.py
#
# Created by Brett H. Andrews on 12 Jun 2017.

# Read in filelists
# Create output dirs
# copy files

import os
from os.path import join

import click


@click.command()
@click.option('--path_mzr', default=None)
@click.option('--stackname', default='dr7_M0.1e')
@click.option('--binnames', default=None, multiple=True)
def copy_spectra(stackname, binnames, path_mzr):
    path_data = join('//', 'Volumes', 'My Passport', 'osu', 'andrews', 'projects',
                     'mass-metallicity', 'data', 'raw_FITS_dr7')
    assert os.path.isdir(path_data), f'``path_data`` does not exist: {path_data}'

    path_mzr_default = join(os.path.expanduser('~'), 'projects', 'mzr', 'stacks')
    path_mzr = path_mzr if path_mzr is not None else path_mzr_default
    assert os.path.isdir(path_mzr), f'``path_mzr`` does not exist: {path_mzr}'

    path_stack = join(path_mzr, stackname)
    assert os.path.isdir(path_stack), f'``path_stack`` does not exist: {path_stack}'
    
    if not binnames:
        binnames = [it for it in os.listdir(path_stack) if it[0] == 'M']
    
    for binname in binnames:
        # read filelist
        # reconstruct filepaths
        pass
    
    
    

    

