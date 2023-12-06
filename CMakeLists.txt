cmake_minimum_required(VERSION 3.14)

# set the project name
project(gbdes)

include(CheckLibraryExists)

# specify the C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

find_library(CFITSIO cfitsio REQUIRED)

find_package(YAML-CPP REQUIRED)
add_compile_definitions(USE_YAML=TRUE)

find_library(FFTW3 fftw3 REQUIRED)

# find astshim library using EUPS package environent variable as hint
find_library(ASTSHIM_LIBRARIES
astshim REQUIRED
HINTS $ENV{ASTSHIM_DIR}/lib
)
#find astshim include directory
find_path(
ASTSHIM_INCLUDES
astshim.h REQUIRED
HINTS $ENV{ASTSHIM_DIR}/include
)

# TMV or Eigen must be available, but TMV is preferred
find_package(TMV QUIET)
if(TMV_FOUND)
    add_compile_definitions(USE_TMV=TRUE)
else()
    find_package(Eigen3 QUIET)
    if(NOT Eigen3_FOUND)
        message(FATAL_ERROR "Either TMV or Eigen3 must be available")
    endif()
    add_compile_definitions(USE_EIGEN=TRUE)
endif()   

add_library(gbdes SHARED
        gbutil/src/BinomFact.cpp gbutil/src/Expressions.cpp gbutil/src/Lookup1d.cpp gbutil/src/Poly2d.cpp
        gbutil/src/Shear.cpp gbutil/src/Table.cpp gbutil/src/odeint.cpp gbutil/src/Interpolant.cpp
        gbutil/src/Poisson.cpp gbutil/src/Pset.cpp gbutil/src/StringStuff.cpp gbutil/src/fft.cpp
        gbfits/FITS.cpp gbfits/FTableExpression.cpp gbfits/FitsTable.cpp gbfits/Header.cpp gbfits/Image.cpp
        gbfits/FTable.cpp gbfits/FitsImage.cpp gbfits/Hdu.cpp gbfits/HeaderFromStream.cpp
        astrometry/src/Astrometry.cpp astrometry/src/PiecewiseMap.cpp astrometry/src/PixelMapCollection.cpp
        astrometry/src/SubMap.cpp astrometry/src/Wcs.cpp astrometry/src/PixelMap.cpp
        astrometry/src/PolyMap.cpp astrometry/src/TemplateMap.cpp astrometry/src/YAMLCollector.cpp
        photometry/src/PhotoMap.cpp photometry/src/PhotoMapCollection.cpp photometry/src/PhotoPiecewise.cpp
        photometry/src/PhotoTemplate.cpp photometry/src/SubMap.cpp
        src/subs/Accum.cpp src/subs/ExtensionAttribute.cpp src/subs/MapFit.cpp src/subs/PhotoPrior.cpp
        src/subs/WCSFitRoutine.cpp src/subs/FitSubroutines.cpp src/subs/Match.cpp src/subs/ReadPhotoPriors.cpp
        src/subs/WCSFoFRoutine.cpp src/subs/DECamInfo.cpp src/subs/MapDegeneracies.cpp src/subs/PhotoMatch.cpp
        src/subs/TPVMap.cpp src/subs/WcsSubs.cpp
        )

target_include_directories(gbdes PUBLIC
        ${PROJECT_SOURCE_DIR}/gbutil/include
        ${PROJECT_SOURCE_DIR}/gbfits
        ${PROJECT_SOURCE_DIR}/astrometry/include
        ${ASTSHIM_INCLUDES}
        ${PROJECT_SOURCE_DIR}/photometry/include
        ${PROJECT_SOURCE_DIR}/include
        )

target_link_libraries(gbdes PUBLIC Eigen3::Eigen ${ASTSHIM_LIBRARIES} ${FFTW3} yaml-cpp cfitsio )

add_executable(WCSFit src/WCSFit.cpp)
target_link_libraries(WCSFit PUBLIC gbdes)
install(TARGETS WCSFit DESTINATION bin)

add_subdirectory(pydir)

