# stack.py
#
# Created by Brett H. Andrews on 13 Jun 2017.


import os
from os.path import join

import numpy as np
import sklearn


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
path_dr7 = join(path_mzr, 'stacks', 'dr7_M0.1e')

"""Read in regridded spectra from file"""


n_samples = 50  # Remove................................................................................
stacks = np.zeros((n_samples, len(grid)))

for ii in range(n_samples):
    spectra_resampled = sklearn.utils.resample(spectra)
    stacks[ii] = np.mean(spectra_resampled, axis=0)
    
spec_mean = np.mean(stacks, axis=0)
spec_std = np.std(stacks, ddof=1, axis=0) 

# Plot
import matplotlib.pyplot as plt
import seaborn as sns

fig, ax = plt.subplots()
ax.plot(grid, spec_mean)
ax.plot(grid, spec_mean + spec_std)
ax.plot(grid, spec_mean - spec_std)

