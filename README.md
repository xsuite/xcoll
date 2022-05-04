# xcoll

<!---![PyPI - Python Version](https://img.shields.io/pypi/pyversions/xcoll?logo=PyPI?style=plastic) ![PyPI - Wheel](https://img.shields.io/pypi/wheel/xcoll?logo=PyPI?style=plastic)-->

![GitHub release (latest by date)](https://img.shields.io/github/v/release/xsuite/xcoll?style=plastic) ![GitHub](https://img.shields.io/github/license/xsuite/xcoll?style=plastic) ![GitHub all releases](https://img.shields.io/github/downloads/xsuite/xcoll/total?logo=GitHub&style=plastic) ![GitHub issues](https://img.shields.io/github/issues/xsuite/xcoll?logo=GitHub&style=plastic) ![GitHub pull requests](https://img.shields.io/github/issues-pr/xsuite/xcoll?logo=GitHub&style=plastic) ![GitHub repo size](https://img.shields.io/github/repo-size/xsuite/xcoll?logo=GitHub&style=plastic)

Collimation in xtrack simulations

## Description

## Getting Started

### Dependencies

* python >= 3.8
    * numpy
    * pandas
    * xsuite (in particular xobjects, xdeps, xtrack, xpart)
* to use K2:
    * gfortran 

### Installing
`xcoll` is packaged using `poetry`, and can be easily installed with `pip`:
```bash
pip install xcoll
```
For a local installation, clone and install in editable mode (need to have `pip` >22):
```bash
git clone git@github.com:xsuite/xcoll.git
pip install -e xcoll
```

### Using K2
To be able to use the K2 scattering algorithms, these need to be compiled from source.
There is a small script that does this (you can ignore the warnings):
```bash
cd xcoll
./compile_K2.sh
```
This installs a shared library in the package that is tailored to your current python installation.

Without compilation, K2 Collimators can be installed in a `Line`, but not tracked.

### Example

## Features

## Authors

* [Frederik Van der Veken](https://github.com/freddieknets) (frederik@cern.ch)
* [Despina Demetriadou](https://github.com/ddemetriadou)
* [Andrey Abramov](https://github.com/anabramo)
* [Giovanni Iadarola](https://github.com/giadarol)


## Version History

* 0.1
    * Initial Release

## License

This project is licensed under the  Apache License 2.0 - see the LICENSE file for details
