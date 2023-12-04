cmake_minimum_required(VERSION 3.14)
project(gbdes_py)

find_package(pybind11 REQUIRED)

set(PYTHON_MODULE_EXTENSION ".so")

pybind11_add_module(wcsfit wcsfit.cpp)

target_link_libraries(wcsfit PUBLIC subs cfitsio ${ASTSHIM_LIBRARIES})

target_include_directories(wcsfit PUBLIC ${PROJECT_SOURCE_DIR}/include ${ASTSHIM_INCLUDES})

install(TARGETS wcsfit DESTINATION "${PROJECT_SOURCE_DIR}")
