#   $Id: Makefile,v 1.3 2012/01/01 18:56:10 garyb Exp $ 
SUBDIRS = 

INCLUDES = 

CXXFLAGS = $(OPTFLAGS) $(ABS_INCLUDES) $(INCLUDES)

SRC = $(shell ls *.cpp)

OBJ = BinomFact.o fft.o StringStuff.o Poisson.o Table.o Pset.o odeint.o Poly2d.o Interpolant.o

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

.PHONY: all install dist depend clean 
