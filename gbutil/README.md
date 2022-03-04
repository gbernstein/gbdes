# gbutil
Gary's C++ and Python utility routines

## Description
Miscellaneous classes and functions of wide utility.  Python routines include some specifications of the DECam focal plane properties. 

## Prerequisites
* A C++ compiler that is compliant with C++-11 standards.
* [FFTW](www.fftw.org) Fourier transform library, if you want to use the `fft` class in C++.
* [yaml-cpp](https://github.com/jbeder/yaml-cpp) YAML-format reader by Jay Bedel, if you want to build the `Lookup1d` and `Poly2d` classes.  As of this writing (Mar 2017), the releases all depend on Boost but the master branch does not, so I use the master branch.
* A linear algebra library, either
 * [TMV](https://github.com/rmjarvis/tmv) by Mike Jarvis, -_OR_-
 * [Eigen](https://eigen.tuxfamily.org).  
 * Either of these can make good use of the [Intel Math Kernel Library](https://software.intel.com/en-us/intel-mkl)
 * The Python libraries assume that [astropy](http://www.astropy.org) is available, and the routines to acquire Gaia data make use of [astroquery](https://astroquery.readthedocs.io/en/latest).

## Installation:
### C++
You of course need to install the external libraries TMV or Eigen, and (optionally) MKL. And the FFTW and yaml-cpp libraries, should you want to use the utilities that depend on them.
The Makefile requires these environment variables to be set:
* `CXX`: path to the C++-11-compliant compiler/linker.
* `CXXFLAGS`: options given to the compiler/linker.  Be sure to specify C++-11 (or higher) compliance if it's not the default, and an optimizer (the linear algebra libraries rely heavily on compile-time optimization).
* `FFTW_DIR`: (optional) path to root of the FFTW installation, should have `/include` and `/lib` subdirectories.
* `YAML_DIR`: (optional) path to root of the yaml-cpp installation
  * `TMV_DIR`: path to root of the TMV installation
  *  -*OR*-
  * `EIGEN_DIR`: path to root of Eigen installation. This is header-only library; there should be a `${EIGEN_DIR}/Eigen` subdirectory where the headers are found.  If both `TMV_DIR` and `EIGEN_DIR` are given, TMV will be used.
* `MKL_DIR`: Path to the MKL.  This is optional, and should be given if you want to use Eigen with MKL.  TMV usage of MKL will have been determined when you built it.

Once these are all set you should be able to just run `make` to build the C++ programs. The few undocumented test routines are built by `make tests`.

### Python
Either `python setup.py install` or just `make python` will install the Python package `gbutil` (which include that subpackage `decam`).

## Use
See the individual utilities' codes for documentation on what they do.
