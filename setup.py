from setuptools import setup, find_packages

setup(
    name='stacking',
    version='0.1',
    packages=find_packages(),
    include_package_data=True,
    install_requires=[
        'click',
        'numpy',
        'scipy',
        'pandas',
        'sklearn',
        'astropy'
    ],
    entry_points='''
        [console_scripts]
        generate_filepaths=stacking.paths:generate_filepaths
        process_spectra=stacking.spectral_processing:process_spectra
    ''',
)