cmake_minimum_required(VERSION 3.14)

# specify the C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# set the project name
project(gbfits)

find_library(CFITSIO cfitsio REQUIRED)

find_package(Eigen3)
find_package(yaml-cpp)
find_package(gbutil REQUIRED)

add_library(${PROJECT_NAME} SHARED FITS.cpp FTableExpression.cpp FitsTable.cpp Header.cpp Image.cpp 
            FTable.cpp FitsImage.cpp Hdu.cpp HeaderFromStream.cpp)

target_include_directories(${PROJECT_NAME} INTERFACE 
    $<INSTALL_INTERFACE:include>
    $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}>)
target_include_directories(${PROJECT_NAME} PRIVATE ${PROJECT_SOURCE_DIR})
target_link_libraries(${PROJECT_NAME} gbutil::gbutil)

# installation
include(GNUInstallDirs)

install(TARGETS ${PROJECT_NAME} EXPORT ${PROJECT_NAME}-config)

install(
    EXPORT ${PROJECT_NAME}-config 
    NAMESPACE ${PROJECT_NAME}::                                                         
    DESTINATION ${CMAKE_INSTALL_DATAROOTDIR}/cmake/${PROJECT_NAME})

install(FILES
        ${PROJECT_SOURCE_DIR}/FITS.h
        ${PROJECT_SOURCE_DIR}/FTable.h
        ${PROJECT_SOURCE_DIR}/FTable2.h
        ${PROJECT_SOURCE_DIR}/FitsImage.h
        ${PROJECT_SOURCE_DIR}/FitsTable.h
        ${PROJECT_SOURCE_DIR}/FitsTypes.h
        ${PROJECT_SOURCE_DIR}/Hdu.h
        ${PROJECT_SOURCE_DIR}/Header.h
        ${PROJECT_SOURCE_DIR}/HeaderFromStream.h
        ${PROJECT_SOURCE_DIR}/Image.h
        ${PROJECT_SOURCE_DIR}/Image2.h
        ${PROJECT_SOURCE_DIR}/ImageAlgorithms.h
        ${PROJECT_SOURCE_DIR}/ReadLdacHeader.h
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})