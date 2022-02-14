add_library(subs Accum.cpp ExtensionAttribute.cpp MapFit.cpp PhotoPrior.cpp WCSFitRoutine.cpp
            CMakeLists.txt FitSubroutines.cpp Match.cpp ReadPhotoPriors.cpp WCSFoFRoutine.cpp
            DECamInfo.cpp MapDegeneracies.cpp PhotoMatch.cpp TPVMap.cpp WcsSubs.cpp)

target_link_libraries(subs PUBLIC gbutil::gbutil astrometry::astrometry photometry::photometry 
                      gbfits::gbfits)
target_include_directories(subs PUBLIC ${PROJECT_SOURCE_DIR}/include)
