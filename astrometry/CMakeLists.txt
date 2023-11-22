add_library(astrometry src/Astrometry.cpp src/PiecewiseMap.cpp src/PixelMapCollection.cpp src/SubMap.cpp 
            src/Wcs.cpp src/PixelMap.cpp src/PolyMap.cpp src/TemplateMap.cpp src/YAMLCollector.cpp)

target_include_directories(astrometry PUBLIC ${PROJECT_SOURCE_DIR}/astrometry/include ${ASTSHIM_INCLUDES})
target_link_libraries(astrometry PUBLIC Eigen3::Eigen ${FFTW3_LIBRARIES} yaml-cpp gbutil)
