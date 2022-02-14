# These should come from environment:
# CXX
# CXXFLAGS
CXXFLAGS = -fPIC
# CFITSIO_DIR

INCLUDES := 

LIBS := -lm

# Collect the includes and libraries we need
ifdef CFITSIO_DIR
INCLUDES += -I $(CFITSIO_DIR)/include
LIBS += -L $(CFITSIO_DIR)/lib -lcfitsio
else
$(error Require CFITSIO_DIR in environment)
endif

# Collect the includes and libraries we need
ifdef GBUTIL_DIR
INCLUDES += -I $(GBUTIL_DIR)/include
GBUTIL_OBJ := $(GBUTIL_DIR)/obj
else
$(error Require GBUTIL_DIR in environment)
endif

# Rule for compilation:
%.o: %.cpp
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c $< -o $@

# External directories where we'll need to clean/build dependents
EXTDIRS = $(GBUTIL_DIR)

OBJ = FITS.o Header.o Hdu.o FitsTable.o FTable.o FTableExpression.o \
	Image.o FitsImage.o HeaderFromStream.o
SRC = $(OBJ:%.o=%.cpp)

EXTOBJ := $(GBUTIL_OBJ)/StringStuff.o $(GBUTIL_OBJ)/Expressions.o

all: $(OBJ)

###############################################################
## Standard stuff:
###############################################################

exts:
	for dir in $(EXTDIRS); do (cd $$dir && $(MAKE)); done

local-depend:
	$(CXX) $(CXXFLAGS) $(INCLUDES) -MM $(SRC) > .depend

depend: local-depend
	for dir in $(EXTDIRS); do (cd $$dir && $(MAKE) depend); done

local-clean:
	rm -f *.o *~ core .depend tests/*.o

clean: local-clean
	for dir in $(EXTDIRS); do (cd $$dir && $(MAKE) clean); done

ifeq (.depend, $(wildcard .depend))
include .depend
endif

export

.PHONY: all install dist depend clean 
