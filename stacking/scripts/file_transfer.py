# copy_spectra.py
#
# Created by Brett H. Andrews on 12 Jun 2017.

# Read in filelists
# Create output dirs
# copy files


import os
from os.path import join

import click
import pandas as pd


@click.command()
@click.option('--stackname', default='dr7_M0.1e')
@click.option('--binnames', default=None, multiple=True)
@click.option('--path_mzr', default=None)
def copy_spectra(stackname, binnames, path_mzr):
    """Copy spectra from long term storage.
    
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
    path_data = join('//', 'Volumes', 'My Passport', 'osu', 'andrews', 'projects',
                     'mass-metallicity', 'data', 'raw_FITS_dr7')
    # assert os.path.isdir(path_data), f'``path_data`` does not exist: {path_data}'

    path_mzr_default = join(os.path.expanduser('~'), 'projects', 'mzr')
    path_mzr = path_mzr if path_mzr is not None else path_mzr_default
    assert os.path.isdir(path_mzr), f'``path_mzr`` does not exist: {path_mzr}'

    path_stack = join(path_mzr, 'stacks', stackname)
    assert os.path.isdir(path_stack), f'``path_stack`` does not exist: {path_stack}'
    
    if not binnames:
        with open(join(path_stack, 'auxiliary', 'binnames.txt'), 'r') as fin:
            binnames = [line.strip() for line in fin]

    
    for binname in binnames:
        if 'M7.2' in binname:
            path_filelist = join(path_stack, binname.split('_n')[0], 'filelist',
                                 binname + '_filenames.txt')

            import ipdb; ipdb.set_trace()
            files = pd.read_csv(path_filelist, delim_whitespace=True, header=None, usecols=(8,))
            files.columns = ['paths']
            dir_ext_hd, filename = files['paths'][0].split('/')[-2:]
            pid, mjd = dir_ext_hd.split('-')[-2:]
            
            path_ext_hd = join(path_data, dir_ext_hd, filename)
            
            path_dr7 = join(path_mzr, 'data', 'raw_fits_dr7')
            path_mjd_pid = join(path_dr7, mjd, pid)
            if not os.path.isdir(path_mjd_pid):
                os.makedirs(path_mjd_pid)
            
            path_local = join(path_mjd_pid, filename)
            
            # shutil.copy(path_ext_hd, path_local)

    
    
    

    

