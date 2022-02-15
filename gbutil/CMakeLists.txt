add_library(gbutil src/BinomFact.cpp src/Expressions.cpp src/Lookup1d.cpp src/Poly2d.cpp src/Shear.cpp 
            src/Table.cpp src/odeint.cpp src/Interpolant.cpp src/Poisson.cpp src/Pset.cpp src/StringStuff.cpp
            src/fft.cpp)

target_include_directories(gbutil PUBLIC ${PROJECT_SOURCE_DIR}/gbutil/include)
target_link_libraries(gbutil PUBLIC Eigen3::Eigen ${FFTW3_LIBRARIES} yaml-cpp)
