add_library(${PROJECT_NAME} SHARED BinomFact.cpp Expressions.cpp Lookup1d.cpp Poly2d.cpp Shear.cpp 
    Table.cpp odeint.cpp Interpolant.cpp Poisson.cpp Pset.cpp StringStuff.cpp fft.cpp)

target_include_directories(${PROJECT_NAME} INTERFACE 
    $<INSTALL_INTERFACE:include>
    $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>)
target_include_directories(${PROJECT_NAME} PRIVATE ${PROJECT_SOURCE_DIR}/include)
target_link_libraries(${PROJECT_NAME}  Eigen3::Eigen ${FFTW3_LIBRARIES} yaml-cpp)

# installation
include(GNUInstallDirs)

install(TARGETS ${PROJECT_NAME} EXPORT ${PROJECT_NAME}-config)

install(
    EXPORT ${PROJECT_NAME}-config 
    NAMESPACE ${PROJECT_NAME}::                                                         
    DESTINATION ${CMAKE_INSTALL_DATAROOTDIR}/cmake/${PROJECT_NAME})

#install(DIRECTORY ${PROJECT_SOURCE_DIR}/include
#        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}../ )
install(FILES
        ${PROJECT_SOURCE_DIR}/include/AndersonDarling.h
        ${PROJECT_SOURCE_DIR}/include/Bounds.h
        ${PROJECT_SOURCE_DIR}/include/Interpolant.h 
        ${PROJECT_SOURCE_DIR}/include/NameIndex.h  
        ${PROJECT_SOURCE_DIR}/include/Shear.h       
        ${PROJECT_SOURCE_DIR}/include/Stopwatch.h         
        ${PROJECT_SOURCE_DIR}/include/odeint.h
        ${PROJECT_SOURCE_DIR}/include/Array.h                  
        ${PROJECT_SOURCE_DIR}/include/Brent.h         
        ${PROJECT_SOURCE_DIR}/include/Legendre.h       
        ${PROJECT_SOURCE_DIR}/include/Poisson.h    
        ${PROJECT_SOURCE_DIR}/include/Simpson.h     
        ${PROJECT_SOURCE_DIR}/include/StringStuff.h
        ${PROJECT_SOURCE_DIR}/include/Array3d.h                
        ${PROJECT_SOURCE_DIR}/include/Expressions.h   
        ${PROJECT_SOURCE_DIR}/include/LinearAlgebra.h  
        ${PROJECT_SOURCE_DIR}/include/Poly2d.h     
        ${PROJECT_SOURCE_DIR}/include/Solve.h       
        ${PROJECT_SOURCE_DIR}/include/SymmetricUpdater.h
        ${PROJECT_SOURCE_DIR}/include/AstronomicalConstants.h  
        ${PROJECT_SOURCE_DIR}/include/Expressions2.h  
        ${PROJECT_SOURCE_DIR}/include/Lookup1d.h       
        ${PROJECT_SOURCE_DIR}/include/Pset.h       
        ${PROJECT_SOURCE_DIR}/include/Statistics.h  
        ${PROJECT_SOURCE_DIR}/include/Table.h
        ${PROJECT_SOURCE_DIR}/include/BinomFact.h              
        ${PROJECT_SOURCE_DIR}/include/Function1d.h    
        ${PROJECT_SOURCE_DIR}/include/Marquardt.h      
        ${PROJECT_SOURCE_DIR}/include/Random.h     
        ${PROJECT_SOURCE_DIR}/include/Std.h         
        ${PROJECT_SOURCE_DIR}/include/fft.h
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
