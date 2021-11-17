# These site-dependent items should be defined in environment:

# CXX = g++-6 -fopenmp
CXXFLAGS = -std=c++11 -fPIC

# TMV_DIR
# CFITSIO_DIR
# FFTW_DIR -or- EIGEN_DIR (former takes precedence)
# YAML_DIR
# GBUTILS_DIR
# GBFITS_DIR
# ASTROMETRY_DIR
# PHOTOMETRY_DIR

# MKL_DIR (optional, used with Eigen)

INCLUDES := 

LIBS := -lm

EXTDIRS := 

# Collect the includes and libraries we need
ifdef FFTW_DIR
INCLUDES += -I $(FFTW_DIR)/include
LIBS += -L $(FFTW_DIR)/lib -lfftw3
else
$(error Require FFTW_DIR in environment)
endif

ifdef YAML_DIR
INCLUDES += -I $(YAML_DIR)/include -D USE_YAML
LIBS += -L $(YAML_DIR)/lib -lyaml-cpp
else
$(error Require YAML_DIR in environment)
endif

ifdef CFITSIO_DIR
INCLUDES += -I $(CFITSIO_DIR)/include
LIBS += -L $(CFITSIO_DIR)/lib -lcfitsio
else
$(error Require CFITSIO_DIR in environment)
endif

ifneq ($(origin CFITSIO_NEEDS_CURL),undefined)
LIBS += -lcurl
endif

ifdef GBUTIL_DIR
INCLUDES += -I $(GBUTIL_DIR)/include
EXTDIRS += $(GBUTIL_DIR)
GBUTIL_OBJ = $(GBUTIL_DIR)/obj
else
$(error Require GBUTIL_DIR in environment)
endif

ifdef GBFITS_DIR
INCLUDES += -I $(GBFITS_DIR)
EXTDIRS += $(GBFITS_DIR)
else
$(error Require GBFITS_DIR in environment)
endif

ifdef ASTROMETRY_DIR
INCLUDES += -I $(ASTROMETRY_DIR)/include
EXTDIRS += $(ASTROMETRY_DIR)
ASTROMETRY_OBJ := $(ASTROMETRY_DIR)/obj
else
$(error Require ASTROMETRY_DIR in environment)
endif

ifdef PHOTOMETRY_DIR
INCLUDES += -I $(PHOTOMETRY_DIR)/include
EXTDIRS += $(PHOTOMETRY_DIR)
PHOTOMETRY_OBJ := $(PHOTOMETRY_DIR)/obj
else
$(error Require PHOTOMETRY_DIR in environment)
endif

# Use either TMV or EIGEN, not both (prefer TMV)
ifdef TMV_DIR
INCLUDES += -I $(TMV_DIR)/include -D USE_TMV
LIBS += $(shell cat $(TMV_DIR)/share/tmv/tmv-link) -ltmv_symband 
else 
ifdef EIGEN_DIR
INCLUDES += -I $(EIGEN_DIR) -D USE_EIGEN
endif
endif

# Check that either TMV or EIGEN are available (ok to have both)
$(if $(or $(TMV_DIR),$(EIGEN_DIR)), , $(error Need either TMV_DIR or EIGEN_DIR))

ifdef MKL_DIR
INCLUDES += -I $(MKL_DIR)/include -D USE_MKL
LIBS += -L$(MKL_DIR)/lib/intel64 -Wl,--no-as-needed -lmkl_intel_lp64 -lmkl_gnu_thread -lmkl_core -lgomp -lpthread -lm -ldl
CXXFLAGS += -m64
endif

# Object files found in external packages:
EXTOBJS =$(GBUTIL_OBJ)/BinomFact.o $(GBUTIL_OBJ)/StringStuff.o $(GBUTIL_OBJ)/Interpolant.o \
	$(GBUTIL_OBJ)/fft.o $(GBUTIL_OBJ)/Table.o $(GBUTIL_OBJ)/Pset.o \
	$(GBUTIL_OBJ)/Poly2d.o $(GBUTIL_OBJ)/Expressions.o $(GBUTIL_OBJ)/Shear.o \
	$(GBUTIL_OBJ)/Lookup1d.o \
	$(GBFITS_DIR)/FITS.o $(GBFITS_DIR)/Header.o $(GBFITS_DIR)/Hdu.o \
	$(GBFITS_DIR)/FitsTable.o $(GBFITS_DIR)/FTable.o $(GBFITS_DIR)/FTableExpression.o \
	$(GBFITS_DIR)/Image.o $(GBFITS_DIR)/FitsImage.o $(GBFITS_DIR)/HeaderFromStream.o \
	$(ASTROMETRY_OBJ)/PixelMap.o $(ASTROMETRY_OBJ)/Astrometry.o \
	$(ASTROMETRY_OBJ)/PolyMap.o $(ASTROMETRY_OBJ)/SubMap.o \
	$(ASTROMETRY_OBJ)/Wcs.o $(ASTROMETRY_OBJ)/PixelMapCollection.o \
	$(ASTROMETRY_OBJ)/TemplateMap.o $(ASTROMETRY_OBJ)/PiecewiseMap.o \
	$(ASTROMETRY_OBJ)/YAMLCollector.o \
	$(PHOTOMETRY_OBJ)/PhotoMap.o $(PHOTOMETRY_OBJ)/PhotoMapCollection.o \
	$(PHOTOMETRY_OBJ)/SubMap.o $(PHOTOMETRY_OBJ)/PhotoTemplate.o \
	$(PHOTOMETRY_OBJ)/PhotoPiecewise.o

##### 
BINDIR = bin
OBJDIR = obj
SRCDIR = src
SUBDIR = src/subs
INCLUDEDIR = include
TESTDIR = tests
TESTBINDIR = testbin

# INCLUDES can be relative paths, and will not be exported to subdirectory makes.
INCLUDES += -I $(INCLUDEDIR)

# Executable C++ programs
EXECS :=  $(wildcard $(SRCDIR)/*.cpp)
TARGETS := $(EXECS:$(SRCDIR)/%.cpp=$(BINDIR)/%)
OBJS := $(EXECS:$(SRCDIR)/%.cpp=$(OBJDIR)/%.o)
SHARED := $(EXECS:$(SRCDIR)/%.cpp=$(OBJDIR)/%.so)

# Python executables
PYEXECS :=  $(wildcard $(SRCDIR)/*.py)
PYTARGETS :=  $(PYEXECS:$(SRCDIR)/%.py=$(BINDIR)/%.py)
# C++ subroutines
SUBS :=  $(wildcard $(SUBDIR)/*.cpp)
SUBOBJS := $(SUBS:$(SUBDIR)/%.cpp=$(OBJDIR)/%.o)

CP = /bin/cp -p
RM = /bin/rm -f

#######################
# Rules - ?? dependencies on INCLUDES ??
#######################

all: cpp python

cpp: exts $(TARGETS) $(SHARED)

python: $(PYTARGETS)

# No setup.py to do here ... python ./setup.py install

# Compilation
$(OBJS):  $(OBJDIR)/%.o : $(SRCDIR)/%.cpp
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c $< -o $@

$(SUBOBJS): $(OBJDIR)/%.o : $(SUBDIR)/%.cpp 
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c $< -o $@

# Linking
$(TARGETS): $(BINDIR)/% : $(OBJDIR)/%.o $(SUBOBJS) $(EXTOBJS)
	$(CXX) $(CXXFLAGS) $^  $(LIBS) -o $@

$(SHARED): $(OBJDIR)/%.so : $(SUBOBJS) $(EXTOBJS)
	$(CXX) -shared -fPIC -o $(OBJDIR)/libgbdes.so $(SUBOBJS) $(EXTOBJS)


# Python executables - copy into bin directory
$(PYTARGETS): $(BINDIR)/% : $(SRCDIR)/%
	$(CP) $^ $@

# External objects - call external makefiles.  ?? Conditional??
# Semicolon prevents calling default rule for .o/.cpp
$(EXTOBJS): exts ;

######### Test programs

TESTSRC := $(wildcard $(TESTDIR)/*.cpp)
TESTINCLUDE := -I $(TESTDIR)
TESTOBJS := $(TESTSRC:$(TESTDIR)/%.cpp=$(OBJDIR)/%.o)
TESTTARGETS := $(TESTSRC:$(TESTDIR)/%.cpp=$(TESTBINDIR)/%)
TESTSPY := $(wildcard $(TESTDIR)/*.py)

tests: $(TESTTARGETS)

$(TESTOBJS):  $(OBJDIR)/%.o : $(TESTDIR)/%.cpp
	$(CXX) $(CXXFLAGS) $(INCLUDES) $(TESTINCLUDE) -c $^ -o $@

$(TESTTARGETS): $(TESTBINDIR)/% : $(OBJDIR)/%.o $(SUBOBJS) $(EXTOBJS)
	$(CXX) $(CXXFLAGS) $^  $(LIBS) -o $@

###############################################################
## Standard stuff:
###############################################################

exts:
	for dir in $(EXTDIRS); do (cd $$dir && $(MAKE)); done

depend: local-depend
	for dir in $(EXTDIRS); do (cd $$dir && $(MAKE) depend); done

local-depend:
	$(RM) .depend
	for src in $(SUBS:$(SUBDIR)/%.cpp=%); \
	 do $(CXX) $(CXXFLAGS) $(INCLUDES) -MM $(SUBDIR)/$$src.cpp -MT obj/$$src.o >> .depend; \
	done
	for src in $(EXECS:%.cpp=%); \
	 do $(CXX) $(CXXFLAGS) $(INCLUDES) -MM $$src.cpp -MT obj/$$src.o >> .depend; \
        done

clean: local-clean
	for dir in $(EXTDIRS); do (cd $$dir && $(MAKE) clean); done

local-clean:
	rm -rf $(OBJDIR)/*.o $(OBJDIR)/*.so $(BINDIR)/* $(TESTBINDIR)/* *~ *.dvi *.aux core .depend

ifeq (.depend, $(wildcard .depend))
include .depend
endif

.PHONY: all install dist depend clean local-clean local-depend exts tests cpp python

