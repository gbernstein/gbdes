cmake_minimum_required(VERSION 3.14)
project(gbdes_py)

find_package(pybind11 REQUIRED)

set(PYTHON_MODULE_EXTENSION ".so")

#pybind11_add_module(example testpybind.cpp)
pybind11_add_module(wcsfit wcsfit.cpp)

#target_link_libraries(example PUBLIC subs gbutil::gbutil astrometry::astrometry ${FFTW3_LIBRARIES}
target_link_libraries(wcsfit PUBLIC subs gbutil::gbutil astrometry::astrometry ${FFTW3_LIBRARIES}
                      gbfits::gbfits cfitsio yaml-cpp)
target_include_directories(wcsfit PUBLIC ${PROJECT_SOURCE_DIR}/include)

#install(TARGETS example DESTINATION "${PROJECT_SOURCE_DIR}")
install(TARGETS wcsfit DESTINATION "${PROJECT_SOURCE_DIR}")
