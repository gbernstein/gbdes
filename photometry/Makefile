#Following should be in environment:
#CXX = g++-4 -fopenmp
#TMV_DIR
#OPTFLAGS
export CXX
export OPTFLAGS

# ABS_INCLUDES are absolute paths 
ABS_INCLUDES = -I $(TMV_DIR)/include -I $(BOOST_DIR)/include

SUBDIRS = 

INCLUDES = -I ../utilities

CXXFLAGS = $(OPTFLAGS) $(ABS_INCLUDES) $(INCLUDES)

SRC = $(shell ls *.cpp)

OBJ = PhotoMap.o PhotoMapCollection.o SubMap.o PhotoPrior.o PhotoMatch.o

# LIB_DIRS 
LIB_DIRS = -L $(TMV_DIR)/lib -L $(BOOST_DIR)/lib
TMV_LINK := $(shell cat $(TMV_DIR)/share/tmv/tmv-link)
LIBS = -lm $(LIB_DIRS) -lboost_regex -ltmv_symband $(TMV_LINK)

all:  $(OBJ)

testTemplate: testTemplate.o PhotoTemplate.o ../utilities/StringStuff.o ../utilities/Poly2d.o $(OBJ)
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
	rm -f *.o *~ *.aux *.log *.dvi core .depend

ifeq (.depend, $(wildcard .depend))
include .depend
endif

export

.PHONY: all install dist depend clean 
