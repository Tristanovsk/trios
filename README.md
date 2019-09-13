# Water color TriOS package

Package to process TriOS-like radiometer data for various acquisition protocols:

- **Above-Water Radiometry** (_awr_): R<sub>rs</sub> (sr<sup>-1</sup>)

<p align="center">
    <img src="images/above_water_system.png" width="400">
</p>

- **In-Water Radiometry** (_iwr_): R<sub>rs</sub> (sr<sup>-1</sup>), K<sub>d</sub> (m<sup>-1</sup>), K<sub>Lu</sub> (m<sup>-1</sup>)

<p align="center">
    <img src="images/in_water_system.png" width="400">
</p>

- **Surface-Water Radiometry** (_swr_): R<sub>rs</sub> (sr<sup>-1</sup>)

<p align="center">
    <img src="images/surface_water_radiometry.png" width="400">
</p>



This package also contains tools for interactive visualization of the radiometric data: 

![animated1](images/visu_trios_data.gif)


## Getting Started

These instructions will get you a copy of the project up and running on your local machine for development and testing purposes.

### Prerequisites

What things you need to install the software and how to install them

```
python3 -m pip install --user --upgrade setuptools
```

### Installing

First, clone `the repository <https://gitlab.irstea.fr/telquel-obs2co/insitu/trios>`__ and execute the following command in the
local copy:

::

    $ python setup.py install

This will install the package into the system's Python path. If another
destination directory is preferred, it can be set by

::

    $ python setup.py install --prefix=<where-to-install>

This installation is supposed to download
and compile all the associated packages as well as prepare the executables `trios_processing` and `trios_visual`.

## Running the tests

```
trios_processing ./test/data/ 150 awr --lat 42.30351823 --lon 9.462897398 --odir ./test/results --no_clobber
trios_processing ./test/data/ 150 awr --lat 42.30351823 --lon 9.462897398 --odir ./test/results --method M99 --name _M99 --plot --figdir ./test/fig
trios_processing ./test/data/ 150 awr --lat 42.30351823 --lon 9.462897398 --odir ./test/results --method osoaa --name _osoaa --plot
```


## Authors

* **Tristan Harmel** - *Initial work* -

See also the list of [contributors](https://github.com/your/project/contributors) who participated in this project.

## License

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details

## Acknowledgments

This work has been partially supported by the _Programme National de Télédétection Spatiale_ (PNTS,
http://www.insu.cnrs.fr/pnts ), grant n°PNTS-2019-13