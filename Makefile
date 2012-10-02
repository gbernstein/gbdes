# $Id: Makefile,v 1.20 2012/01/02 19:25:58 garyb Exp $

# First are site-dependent items
CXX = g++-4 -fopenmp

# OPTFLAGS will be exported for subdir makes
OPTFLAGS = -O3 -DASSERT

# ABS_INCLUDES are absolute paths that will be exported to subdirectories
ABS_INCLUDES = -I /sw/include -I /usr/local/cfitsio/include -I /usr/local/tmv/include \
	-I /usr/local/fftw/include

# LIB_DIRS 
LIB_DIRS = -L /sw/lib -L /usr/local/cfitsio/lib -L /usr/local/tmv/lib -L /usr/local/fftw/lib

TMV_LINK := $(shell cat /usr/local/tmv/share/tmv/tmv-link)

##### Below here should be site-independent
SUBDIRS = utilities images astrometry

# INCLUDES can be relative paths, and will not be exported to subdirectory makes.
INCLUDES = -I utilities -I images -I astrometry

CXXFLAGS = $(OPTFLAGS) $(ABS_INCLUDES) $(INCLUDES)
SRC = $(shell ls *.cpp)

LIBS = -lm $(LIB_DIRS) -lfftw3 -lcfitsio -ltmv_symband $(TMV_LINK)

SUBOBJ =utilities/BinomFact.o utilities/StringStuff.o utilities/Interpolant.o \
	utilities/fft.o utilities/Table.o utilities/Pset.o utilities/Poly2d.o \
	images/FITS.o images/Header.o images/Hdu.o \
	astrometry/PixelMap.o astrometry/Astrometry.o astrometry/PolyMap.o \
	astrometry/PixelMapCollection.o 

#images/FITSTable.o images/Image.o images/FITSImage.o images/HeaderFromStream.o \

#OBJ =  SCAMPMap.o $(SUBOBJ)
OBJ =  $(SUBOBJ)

all: depend subs

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
###############################################################
## Standard stuff:
###############################################################

export CXX
export OPTFLAGS
export ABS_INCLUDES

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
