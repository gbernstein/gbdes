cmake_minimum_required(VERSION 3.14)

# specify the C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# set the project name
project(photometry)

find_package(yaml-cpp REQUIRED)
add_compile_definitions(USE_YAML=TRUE)

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

find_package(gbutil)

add_subdirectory(src)
