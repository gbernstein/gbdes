#   $Id: Makefile,v 1.6 2012/02/02 23:34:25 garyb Exp $ 
SUBDIRS = 

INCLUDES = -I ../utilities

CXXFLAGS = $(OPTFLAGS) $(ABS_INCLUDES) $(INCLUDES)

SRC = $(shell ls *.cpp)

OBJ = Astrometry.o Ephemeris.o PixelMap.o PolyMap.o PixelMapCollection.o

all: depend $(OBJ)

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
	rm -f *.o *~ *.aux *.log *.dvi core .depend

ifeq (.depend, $(wildcard .depend))
include .depend
endif

export

.PHONY: all install dist depend clean 
