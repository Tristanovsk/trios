from setuptools import setup, find_packages
from grs import __version__

setup(
    name='trios',
    version=__version__,
    packages=find_packages(),
    package_data={'': ['*.so']},
    #     # If any package contains *.txt files, include them:
    #     '': ['*.txt'],
    #     'lut': ['data/lut/*.nc'],
    #     'aux': ['data/aux/*']
    # },
    include_package_data=True,

    url='https://gitlab.irstea.fr/ETL-TELQUEL/etl/tree/dev/preprocessing/trios',
    license='MIT',
    author='T. Harmel',
    author_email='tristan.harmel@ntymail.com',
    description='Package to help process TriOS radiometer data for various above-water or in-water setups',
    # TODO update Dependent packages (distributions)
    install_requires=['dash','dash_core_components','dash_html_components','pandas', 'scipy', 'numpy', 'netCDF4', 'matplotlib', 'docopt', 'GDAL', 'python-dateutil'],

    entry_points={
        'console_scripts': [
            'trios = main:main',
            'visu_trios = visu.data_visu:main'
        ]}
)
