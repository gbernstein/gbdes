add_library(gbfits SHARED FITS.cpp FTableExpression.cpp FitsTable.cpp Header.cpp Image.cpp 
            FTable.cpp FitsImage.cpp Hdu.cpp HeaderFromStream.cpp)

target_include_directories(gbfits PUBLIC ${PROJECT_SOURCE_DIR}/gbfits)
target_link_libraries(gbfits gbutil cfitsio)
