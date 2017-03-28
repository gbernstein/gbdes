# gbdes
Gary's DECam instrumental signature fitting and processing programs

## Description
This repo contains code for deriving photometric and astrometric calibration solutions for complex multi-detector astronomical imagers.  The methods and uses of these codes are quite complex and currently poorly documented.  It's coming along.  My apologies in advance for the unprofessional nature of this repository!

The `src/subs` directory contains C++ code for shared classes and functions.  The `src` contains C++ code for executable programs and Python code meant to be called at the command line.  Right now there are no importable Python modules.  

## Prerequisites
* A C++ compiler that is compliant with C++-11 standards.  The time-intensive codes are multithreaded with OpenMP, so your compiler will need this capability for solutions of typical size.  Note the MacOS clang compiler does not do the latter.
* [FFTW](www.fftw.org) Fourier transform library
* [yaml-cpp](https://github.com/jbeder/yaml-cpp) YAML-format reader by Jay Bedel.  As of this writing (Mar 2017), the releases all depend on Boost but the master branch does not, so I use the master branch.
* [CFITSIO](https://heasarc.gsfc.nasa.gov/fitsio/fitsio.html) FITS library.
* A linear algebra library, either
 * [TMV](https://github.com/rmjarvis/tmv) by Mike Jarvis, -_OR_-
 * [Eigen](https://eigen.tuxfamily.org).  
 * Either of these can make good use of the [Intel Math Kernel Library](https://software.intel.com/en-us/intel-mkl)
* These other repos of mine are compiled into `gbdes` so you need to have them cloned on your machine for builds:
 * [`gbutil`](https://github.com/gbernstein/gbutil) general utilities
 * [`gbfits`](https://github.com/gbernstein/gbfits) C++ interface to CFITSIO
 * [`astrometry`[(https://github.com/gbernstein/astrometry) routines for astrometric transformations and model construction
 * [`photometry`](https://github.com/gbernstein/photometry) routines for photometric model construction
 * [`gmbpy`](https://github.com/gbernstein/gmbpy) Python module has a few utilities called by some of the Python routines in `gbdes`.  You may not need this, which is good since that repo does not yet have proper setup.py and such.

## Installation
You of course need to install the external libraries FFTW, yaml-cpp, CFITSIO, TMV/Eigen, and (optionally) MKL.  The Makefile for `gbdes` will call the Makefiles for the other gbernstein repos as needed, so you don't need to build them explicitly, and they do not (yet) have library files to keep track of.
The Makefile requires these environment variables to be set:
* `CXX`: path to the C++-11-compliant compiler/linker.
* `CXXFLAGS`: options given to the compiler/linker.  Be sure to specify C++-11 (or higher) compliance if it's not the default, OpenMP usage if desired, and an optimizer (the linear algebra libraries rely heavily on compile-time optimization).
* `FFTW_DIR`: path to root of the FFTW installation, should have `/include` and `/lib` subdirectories.
* `YAML_DIR`: path to root of the yaml-cpp installation
* `CFITSIO_DIR`: path to root of the CFITSIO installation
* `TMV_DIR`: path to root of the TMV installation
*  -*OR*-
* `EIGEN_DIR`: path to root of Eigen installation. This is header-only library; there should be a `${EIGEN_DIR}/Eigen` subdirectory where the headers are found.  If both `TMV_DIR` and `EIGEN_DIR` are given, TMV will be used.
* `MKL_DIR`: Path to the MKL.  This is optional, and should be given if you want to use Eigen with MKL.  TMV usage of MKL will have been determined when you built it.
* `GBUTIL_DIR`, `GBFITS_DIR`, `ASTROMETRY_DIR`, and `PHOTOMETRY_DIR` are the directories where you've cloned these repos.

Once these are all set you should be able to just run `make cpp` to build the C++ programs.  `make python` will copy over the Python executables into the bin/ directory.  `make` should do both.

## Use
When the codes are built, the executables of the C++, as well as copies of Python executables, are in the `bin/` directory.  Put this in your path or move them where you please - there is no `make install` yet.
The `LD_LIBRARY_PATH` environment variable will need to be set to reach your FFTW, yaml-cpp, TMV, Eigen, and MKL libraries.


## Notes
* The `tests` programs are not functional yet, don't try to build them.
* So far I have tested on Gnu V6 compiler, with both TMV and Eigen.
* Let me know about the inevitable failures.
