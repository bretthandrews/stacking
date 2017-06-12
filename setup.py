from setuptools import setup, find_packages

setup(
    name='stacking',
    version='0.1',
    packages=find_packages(),
    include_package_data=True,
    install_requires=[
        'numpy',
        'Click',
        'astropy',
    ],
    entry_points='''
        [console_scripts]
        copy_spectra=stacking.scripts.file_transfer:copy_spectra
    ''',
)