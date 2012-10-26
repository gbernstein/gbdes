# These should come from environment:
# CXX
# TMV_DIR
# FFTW_DIR
# OPTFLAGS

export CXX
export OPTFLAGS

ABS_INCLUDES = -I $(FFTW_DIR)/include -I $(TMV_DIR)/include -I $(BOOST_DIR)/include

SUBDIRS = 

INCLUDES = 

CXXFLAGS = $(OPTFLAGS) $(ABS_INCLUDES) $(INCLUDES)

SRC = $(shell ls *.cpp)

OBJ = BinomFact.o fft.o StringStuff.o Poisson.o Table.o Pset.o odeint.o Poly2d.o \
	Interpolant.o Expressions.o

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

export

.PHONY: all install dist depend clean 
