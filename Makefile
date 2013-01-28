# $Id: Makefile,v 1.20 2012/01/02 19:25:58 garyb Exp $

# First are site-dependent items - expect these to be defined in environment!
#CXX = g++-4 -fopenmp
#TMV_DIR
#CFITSIO_DIR
#FFTW_DIR
#BOOST_DIR

export CXX
export TMV_DIR
export CFITSIO_DIR
export BOOST_DIR
export FFTW_DIR

# OPTFLAGS will be exported for subdir makes
OPTFLAGS = -O3 -DASSERT
export OPTFLAGS


# ABS_INCLUDES are absolute paths 
ABS_INCLUDES = -I $(TMV_DIR)/include -I $(CFITSIO_DIR)/include \
	-I $(FFTW_DIR)/include -I $(BOOST_DIR)/include

# LIB_DIRS 
LIB_DIRS = -L $(CFITSIO_DIR)/lib -L $(TMV_DIR)/lib -L $(FFTW_DIR)/lib \
	-L $(BOOST_DIR)/lib

TMV_LINK := $(shell cat $(TMV_DIR)/share/tmv/tmv-link)
#TMV_LINK := $(shell cat $(TMV_DIR)/share/tmv-link)

##### Below here should be site-independent
SUBDIRS = utilities images astrometry

# INCLUDES can be relative paths, and will not be exported to subdirectory makes.
INCLUDES = -I utilities -I images -I astrometry

CXXFLAGS = $(OPTFLAGS) $(ABS_INCLUDES) $(INCLUDES)
SRC = $(shell ls *.cpp)

LIBS = -lm $(LIB_DIRS) -lboost_regex -lfftw3 -lcfitsio -ltmv_symband $(TMV_LINK)

SUBOBJ =utilities/BinomFact.o utilities/StringStuff.o utilities/Interpolant.o \
	utilities/fft.o utilities/Table.o utilities/Pset.o utilities/Poly2d.o \
	utilities/Expressions.o \
	images/FITS.o images/Header.o images/Hdu.o images/FitsTable.o \
	images/FTable.o images/FTableExpression.o \
	astrometry/PixelMap.o astrometry/Astrometry.o astrometry/PolyMap.o \
	astrometry/SubMap.o astrometry/Wcs.o astrometry/PixelMapCollection.o \
	astrometry/SerializeProjection.o

#images/FitsTable.o images/Image.o images/FITSImage.o \

#OBJ =  SCAMPMap.o $(SUBOBJ)
OBJ =  $(SUBOBJ)

all: depend subs

testSerialize: testSerialize.o TPVMap.o $(OBJ)
	$(CXX) $(CXXFLAGS) $^  $(LIBS) -o $@
ApplyWCS: ApplyWCS.o TPVMap.o $(OBJ)
	$(CXX) $(CXXFLAGS) $^  $(LIBS) -o $@
FitsGlue: FitsGlue.o $(OBJ)
	$(CXX) $(CXXFLAGS) $^  $(LIBS) -o $@
WCSFoF: WCSFoF.o TPVMap.o $(OBJ)
	$(CXX) $(CXXFLAGS) $^  $(LIBS) -o $@
WCSFit: WCSFit.o TPVMap.o Match.o $(OBJ)
	$(CXX) $(CXXFLAGS) $^  $(LIBS) -o $@
DumpLDACHeader: DumpLDACHeader.o $(OBJ)
	$(CXX) $(CXXFLAGS) $^  $(LIBS) -o $@
testPrep: testPrep.o $(OBJ)
	$(CXX) $(CXXFLAGS) $^  $(LIBS) -o $@
test: test.o $(OBJ)
	$(CXX) $(CXXFLAGS) $^  $(LIBS) -o $@
testFT: testFT.o $(OBJ)
	$(CXX) $(CXXFLAGS) $^  $(LIBS) -o $@
testFT2: testFT2.o $(OBJ)
	$(CXX) $(CXXFLAGS) $^  $(LIBS) -o $@
testFT3: testFT3.o $(OBJ)
	$(CXX) $(CXXFLAGS) $^  $(LIBS) -o $@
testFT4: testFT4.o $(OBJ)
	$(CXX) $(CXXFLAGS) $^  $(LIBS) -o $@
testFT5: testFT5.o $(OBJ)
	$(CXX) $(CXXFLAGS) $^  $(LIBS) -o $@
testFT6: testFT6.o $(OBJ)
	$(CXX) $(CXXFLAGS) $^  $(LIBS) -o $@
testFT7: testFT7.o $(OBJ)
	$(CXX) $(CXXFLAGS) $^  $(LIBS) -o $@
testFT8: testFT8.o $(OBJ)
	$(CXX) $(CXXFLAGS) $^  $(LIBS) -o $@
###############################################################
## Standard stuff:
###############################################################


subs:
	for dir in $(SUBDIRS); do (cd $$dir && $(MAKE)); done

depend:
	for dir in $(SUBDIRS); do (cd $$dir && $(MAKE) depend); done
	$(CXX) $(CXXFLAGS) -MM $(SRC) > .$@

clean:
	for dir in $(SUBDIRS); do (cd $$dir && $(MAKE) clean); done
	rm -f *.o *~ *.dvi *.aux core .depend

ifeq (.depend, $(wildcard .depend))
include .depend
endif

.PHONY: all install dist depend clean 
