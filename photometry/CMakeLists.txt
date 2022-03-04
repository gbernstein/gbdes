add_library(photometry src/PhotoMap.cpp src/PhotoMapCollection.cpp src/PhotoPiecewise.cpp src/PhotoTemplate.cpp
            src/SubMap.cpp)

target_include_directories(photometry PUBLIC ${PROJECT_SOURCE_DIR}/photometry/include)
target_link_libraries(photometry PUBLIC Eigen3::Eigen ${FFTW3_LIBRARIES} yaml-cpp gbutil)
