from setuptools import setup, find_packages
from trios import __version__

setup(
    name='trios',
    version=__version__,
    packages=find_packages(exclude=['build']),
    package_data={'': ['*.so'],
    #     # If any package contains *.txt files, include them:
    #     '': ['*.txt'],
    #     'lut': ['data/lut/*.nc'],
    #     'aux': ['aux/*']
     },
    include_package_data=True,

    url='https://gitlab.irstea.fr/telquel-obs2co/insitu/trios',
    license='MIT',
    author='T. Harmel',
    author_email='tristan.harmel@gmail.com',
    description='Package to help trios TriOS radiometer data for various above-water or in-water setups',
    # TODO update Dependent packages (distributions)
    install_requires=['cmocean','dash','dash_core_components','dash_html_components','pandas', 'scipy', 'numpy',
                      'pyodbc', 'netCDF4', 'matplotlib', 'docopt', 'GDAL', 'python-dateutil','plotly'],

    entry_points={
        'console_scripts': [
            'trios_processing = trios.main:main',
            'trios_visual = visu.data_visu:main'
        ]}
)
