add_library(${PROJECT_NAME} SHARED PhotoMap.cpp PhotoMapCollection.cpp PhotoPiecewise.cpp PhotoTemplate.cpp
            SubMap.cpp)

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
        ${PROJECT_SOURCE_DIR}/include/PhotoMap.h
        ${PROJECT_SOURCE_DIR}/include/PhotoMapCollection.h
        ${PROJECT_SOURCE_DIR}/include/PhotoPiecewise.h
        ${PROJECT_SOURCE_DIR}/include/PhotoTemplate.h
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
