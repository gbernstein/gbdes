#   $Id: Makefile,v 1.7 2011/02/14 15:22:13 garyb Exp $ 
SUBDIRS = 

INCLUDES = -I ../utilities

CXXFLAGS = $(OPTFLAGS) $(ABS_INCLUDES) $(INCLUDES)

SRC = $(shell ls *.cpp)

OBJ = FITS.o Image.o FITSImage.o HeaderFromStream.o

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
	rm -f *.o *~ core .depend

ifeq (.depend, $(wildcard .depend))
include .depend
endif

export
# don't export INCLUDES??

.PHONY: all install dist depend clean 
