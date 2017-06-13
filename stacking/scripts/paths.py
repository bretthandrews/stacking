# paths.py
#
# Created by Brett H. Andrews on 12 Jun 2017.


import os
from os.path import join

import click
import pandas as pd


@click.command()
@click.option('--stackname', default='dr7_M0.1e')
@click.option('--binnames', default=None, multiple=True)
@click.option('--path_mzr', default=None)
@click.option('--overwrite', '-o', default=False, is_flag=True)
def generate_filepaths(stackname, binnames, path_mzr, overwrite):
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
        overwrite (bool):
            If ``True``, overwrite existing files. Default is ``False``.         
    """
    # data directory on passport external hard drive
    path_data = join('//', 'Volumes', 'My Passport', 'osu', 'andrews', 'projects',
                     'mass-metallicity', 'data', 'raw_FITS_dr7')
    assert os.path.isdir(path_data), f'``path_data`` does not exist: {path_data}'

    # mass-metallicity relation project directory on MacBook Air
    path_mzr_default = join(os.path.expanduser('~'), 'projects', 'mzr')
    path_mzr = path_mzr if path_mzr is not None else path_mzr_default
    assert os.path.isdir(path_mzr), f'``path_mzr`` does not exist: {path_mzr}'

    # stacks directory on MacBook Air
    path_stack = join(path_mzr, 'stacks', stackname)
    assert os.path.isdir(path_stack), f'``path_stack`` does not exist: {path_stack}'

    # read in bin masses (and SFRs) and number of galaxies per bin
    if not binnames:
        with open(join(path_stack, 'auxiliary', 'binnames.txt'), 'r') as fin:
            binnames = [line.strip() for line in fin]

    # create output directory
    path_filelists = join(path_stack, 'filelists')
    if not os.path.isdir(path_filelists):
        os.mkdir(path_filelists)
        click.echo(f'Created directory: {path_filelists}')

    click.echo('Files written:')
    for binname in binnames:
        binpar = binname.split('_n')[0]
        filelist_in = join(path_stack, binpar, 'filelist', binname + '_filenames.txt')

        with open(filelist_in, 'r') as fin:
            paths = [line.split('raw_FITS_dr7/')[1].strip() for line in fin]

        # paths to files on passport external hard drive
        files_passport = [join(path_data, it) for it in paths]

        # write paths to files on passport to a local file
        filelist_out = join(path_filelists, binpar + '.txt')

        if not os.path.isfile(filelist_out) or overwrite:
            
            with open(filelist_out, 'w') as fout:
                for it in files_passport:
                    fout.write(it + '\n')
            
            click.echo(f'{filelist_out}')
        
        else:
            click.echo(f'Not written (overwrite with --overwrite): {filelist_out}')

    

