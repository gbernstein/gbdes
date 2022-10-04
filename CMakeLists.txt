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

add_subdirectory(gbutil)
add_subdirectory(gbfits)
add_subdirectory(astrometry)
add_subdirectory(photometry)

add_subdirectory(src/subs)
add_subdirectory(src)
add_subdirectory(pydir)

