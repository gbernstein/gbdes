# add the executable
add_executable(WCSFit WCSFit.cpp)

target_link_libraries(WCSFit PUBLIC subs gbfits cfitsio)

install(TARGETS WCSFit DESTINATION bin)
