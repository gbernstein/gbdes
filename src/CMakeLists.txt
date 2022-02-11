# add the executable
add_executable(WCSFit WCSFit.cpp)

target_link_libraries(WCSFit PUBLIC subs gbutil::gbutil astrometry::astrometry gbfits::gbfits 
                      cfitsio ${FFTW3_LIBRARIES} yaml-cpp)

install(TARGETS WCSFit DESTINATION bin)
