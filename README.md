# gbdes
Gary's DECam instrumental signature fitting and processing programs

## Description
This repo contains code for deriving photometric and astrometric calibration solutions for complex multi-detector astronomical imagers.  The methods and uses of these codes are quite complex and currently poorly documented.  It's coming along.  My apologies in advance for the unprofessional nature of this repository!

The astrometric and photometric solutions derived and used by this code are stored in YAML format. Python code that can read the astrometric solution files and execute the transformations they specify is available in the [pixmappy repo](https://github.com/gbernstein/pixmappy).

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

See the [Starflat cookbook](https://github.com/gbernstein/gbdes/wiki/Starflat-cookbook) in the wiki for details about using the code to derive astrometric and photometric solutions.  In brief, a
typical workflow for deriving photometric and astrometric homogenization models for a load of exposures would use these programs:
1. Create object catalogs from your exposures with your favorite software and store them into FITS files with one binary table extension for the objects in each detector.  These codes will also look for critical information in the headers of these catalogs, such as a starting WCS (e.g. from _SCAMP_) good enough to permit matching.  Catalogs of reference objects (e.g. Gaia positions) can be included.
2. For DECam images, I run the `compresscat.py` script to filter the catalogs down to the useful stellar objects, get rid of the "LDAC" formatting output by _SExtractor,_ and to install in the headers various useful quantities like aperture corrections and the MJD at the midpoint of the exposure.
3. `configure.py` collects all the meta-data from the catalogs and creates a `.config` file holding FITS binary tables describing all of the exposures, instruments, fields, etc. in the data.
4. `WCSFoF` reads the `.config` file then uses a friends-of-friends matching algorithm to link together all detections of common objects found in distinct exposures.  This information is combined with the `.config` input into an `.fof` file of FITS tables.
5. For DECam data I run `decamDCR.py` on the `.fof` file to calculate airmasses and parallactic angles for each exposure, and to calculate and save the expected differential chromatic refraction (DCR) needed when doing precision astrometry later.
6. I also run `recenter.py` which recalculates the optic-axis RA and Dec for each exposure and saves this into the `.fof` file.
7. `PhotoFit` does the big job of optimizing the parameters of a photometric model to maximize agreement between magnitudes measured in different exposures of the same source.  There are many configuration inputs.  The main output is a `.photo` file specifying the resultant model.
8. `MagColor` can be run once there are `.photo` models for multiple filters. It calculates mean magnitudes and colors for all of the objects by combining the data from all exposures.  This step can be iterated with `PhotoFit` to allow chromatic photometric models.
9. `WCSFit` does the next big job of optimizing the parameters of an astrometric model to maximize agreement among the exposures and any reference catalogs.  Again, very complicated configuration, and the output is an `.astro` file.  

The `allfit.py` script is meant to orchestrate this process but is not yet generally useful.

Once the models have been derived, we can use them in various ways:
* `ApplyPhoto` and `ApplyWCS` take as input an ASCII list of object pixel positions / magnitudes and will output calibrated mags / positions according to the derived `.photo` and `.astro` models, respectively.
* `Photo2DESDM` turns a `.photo` model into a full-resolution image in a format to be used as a flat-field correction ("star flat") by DES Data Managment.
* `DrawAstro` likewise produced pixelized versions of the astrometric solutions.  It can also produce FITS-format WCS headers of the TPV type that approximate the solutions.
* `DrawPhoto` will make an image of what the photometric solution looks like when all the DECam detectors are combined.

Other programs here of some utility for DECam processing:
* `DECamMosaic` will paste together individual-CCD images from DECam into a lower-resolution sky-coordinate FITS image for the full focal plane.
* `ListClipped` can be used to save a list of detections discarded as outliers by `PhotoFit` or `WCSFit`, which can then be used as given to other fitting steps as objects to be ignored.
* `Tpv2Pixmap` converts ASCII versions of the FITS TPV WCS solutions into the more general YAML format used by the `.astro` files.
* `UpdateHeaders` is something I use to add/replace header values in DECam images with new values before processing them.

## Notes
* The `tests` programs are not functional yet, don't try to build them.
* So far I have tested on Gnu V6 compiler, with both TMV and Eigen.
* Let me know about the inevitable failures.
