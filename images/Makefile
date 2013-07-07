# These should come from environment:
# CXX
# TMV_DIR
# FFTW_DIR
# BOOST_DIR
# CFITSIO_DIR
# OPTFLAGS

export CXX
export OPTFLAGS

ABS_INCLUDES = -I $(FFTW_DIR)/include -I $(TMV_DIR)/include \
	-I $(CFITSIO_DIR)/include -I $(BOOST_DIR)/include

SUBDIRS = 

INCLUDES = -I ../utilities

CXXFLAGS = $(OPTFLAGS) $(ABS_INCLUDES) $(INCLUDES)

SRC = $(shell ls *.cpp)

OBJ = FITS.o Header.o Hdu.o FitsTable.o FTable.o FTableExpression.o Image.o FitsImage.o

all: $(OBJ)

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
	rm -f *.o *~ core .depend

ifeq (.depend, $(wildcard .depend))
include .depend
endif

.PHONY: all install dist depend clean 
