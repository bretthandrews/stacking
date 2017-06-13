from setuptools import setup, find_packages

setup(
    name='stacking',
    version='0.1',
    packages=find_packages(),
    include_package_data=True,
    install_requires=[
        'numpy',
        'Click',
        'pandas',
    ],
    entry_points='''
        [console_scripts]
        generate_filepaths=stacking.scripts.paths:generate_filepaths
    ''',
)