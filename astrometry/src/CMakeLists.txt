add_library(${PROJECT_NAME} SHARED Astrometry.cpp PiecewiseMap.cpp PixelMapCollection.cpp SubMap.cpp Wcs.cpp
            PixelMap.cpp PolyMap.cpp TemplateMap.cpp YAMLCollector.cpp)

target_include_directories(${PROJECT_NAME} INTERFACE 
    $<INSTALL_INTERFACE:include>
    $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>)
target_include_directories(${PROJECT_NAME} PRIVATE ${PROJECT_SOURCE_DIR}/include)
target_link_libraries(${PROJECT_NAME} Eigen3::Eigen ${FFTW3_LIBRARIES} yaml-cpp gbutil::gbutil)

# installation
include(GNUInstallDirs)

install(TARGETS ${PROJECT_NAME} EXPORT ${PROJECT_NAME}-config)

install(
    EXPORT ${PROJECT_NAME}-config 
    NAMESPACE ${PROJECT_NAME}::                                                         
    DESTINATION ${CMAKE_INSTALL_DATAROOTDIR}/cmake/${PROJECT_NAME})

install(FILES
        ${PROJECT_SOURCE_DIR}/include/Astrometry.h
        ${PROJECT_SOURCE_DIR}/include/LeapSeconds.h
        ${PROJECT_SOURCE_DIR}/include/PiecewiseMap.h
        ${PROJECT_SOURCE_DIR}/include/PixelMap.h
        ${PROJECT_SOURCE_DIR}/include/PixelMapCollection.h
        ${PROJECT_SOURCE_DIR}/include/PolyMap.h
        ${PROJECT_SOURCE_DIR}/include/TemplateMap.h
        ${PROJECT_SOURCE_DIR}/include/Wcs.h
        ${PROJECT_SOURCE_DIR}/include/YAMLCollector.h
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
