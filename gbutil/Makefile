# These should come from environment:
# CXX
CXXFLAGS = -fPIC

# TMV_DIR
# FFTW_DIR
# YAML_DIR
# MKL_DIR (optionally)

OBJDIR = obj
SRCDIR = src
INCLUDEDIR = include
TESTDIR = tests
TESTBINDIR = testbin

INCLUDES := -I $(INCLUDEDIR)

LIBS := -lm

# Here are (most of) the object libraries we can compile for external use:
OBJS := $(OBJDIR)/BinomFact.o $(OBJDIR)/StringStuff.o $(OBJDIR)/Poisson.o \
        $(OBJDIR)/Table.o $(OBJDIR)/Pset.o $(OBJDIR)/odeint.o \
        $(OBJDIR)/Interpolant.o $(OBJDIR)/Expressions.o $(OBJDIR)/Shear.o \
	$(OBJDIR)/Poly2d.o $(OBJDIR)/Lookup1d.o

# Collect the includes and libraries we need
ifdef FFTW_DIR
INCLUDES += -I $(FFTW_DIR)/include
LIBS += -L $(FFTW_DIR)/lib -lfftw3
OBJS += $(OBJDIR)/fft.o 
else
$(warning: WARNING: No FFTW_DIR in environment, skipping fft.cpp compilation)
endif

ifdef YAML_DIR
INCLUDES += -I $(YAML_DIR)/include -D USE_YAML
LIBS += -L $(YAML_DIR)/lib -lyaml-cpp
else
$(warning WARNING: No YAML_DIR in environment, skipping Lookup1d, some Poly2d functionality)
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

# External directories where we'll need to clean/build dependents
EXTDIRS = 

all: $(OBJS)

python:
	python setup.py install

# Rule for compilation:
$(OBJS):  $(OBJDIR)/%.o : $(SRCDIR)/%.cpp
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c $< -o $@

######### Test programs

TESTSRC := $(wildcard $(TESTDIR)/*.cpp)
TESTINCLUDE := -I $(TESTDIR)
TESTOBJS := $(TESTSRC:$(TESTDIR)/%.cpp=$(OBJDIR)/%.o)
TESTTARGETS := $(TESTSRC:$(TESTDIR)/%.cpp=$(TESTBINDIR)/%)
TESTSPY := $(wildcard $(TESTDIR)/*.py)

tests: $(TESTTARGETS)

$(TESTOBJS):  $(OBJDIR)/%.o : $(TESTDIR)/%.cpp
	$(CXX) $(CXXFLAGS) $(INCLUDES) $(TESTINCLUDE) -c $^ -o $@

$(TESTTARGETS): $(TESTBINDIR)/% : $(OBJDIR)/%.o $(OBJS)
	$(CXX) $(CXXFLAGS) $^  $(LIBS) -o $@

###############################################################
## Standard stuff:
###############################################################

exts:
	for dir in $(EXTDIRS); do (cd $$dir && $(MAKE)); done

local-depend:
	rm -f .depend
	for src in $(OBJS:$(OBJDIR)/%.o=%) ; \
	   do $(CXX) $(CXXFLAGS) $(INCLUDES) -MM $(SRCDIR)/$$src.cpp -MT $(OBJDIR)/$$src.o >> .depend; \
	done

depend: local-depend
	for dir in $(EXTDIRS); do (cd $$dir && $(MAKE) depend); done

local-clean:
	rm -rf $(OBJDIR)/* $(TESTBINDIR)/* *~ core .depend

clean: local-clean
	for dir in $(EXTDIRS); do (cd $$dir && $(MAKE) clean); done

ifeq (.depend, $(wildcard .depend))
include .depend
endif

export

.PHONY: all install dist depend clean local-clean local-depend exts python tests
