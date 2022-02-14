# Makefile for Astrometry package.
# The only executables to be made are the tests, there are no "subs"
# 
# These site-dependent items should be defined in environment:

# CXX
# CXXFLAGS
CXXFLAGS = -fPIC
# TMV_DIR -or- EIGEN_DIR
# GBUTILS_DIR
# YAML_DIR
# MKL_DIR (optional, used w/EIGEN)

INCLUDES := 

LIBS := -lm

EXTDIRS := 

# Collect the includes and libraries we need
ifdef YAML_DIR
INCLUDES += -I $(YAML_DIR)/include -D USE_YAML
LIBS += -L $(YAML_DIR)/lib -lyaml-cpp
else
$(error Require YAML_DIR in environment)
endif

ifdef GBUTIL_DIR
INCLUDES += -I $(GBUTIL_DIR)/include
EXTDIRS += $(GBUTIL_DIR)
GBUTIL_OBJ = $(GBUTIL_DIR)/obj
else
$(error Require GBUTIL_DIR in environment)
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

# Object files found in external packages, :
EXTOBJS =$(GBUTIL_OBJ)/StringStuff.o $(GBUTIL_OBJ)/Poly2d.o $(GBUTIL_OBJ)/Lookup1d.o

##### 
BINDIR = bin
OBJDIR = obj
SRCDIR = src
SUBDIR = src
INCLUDEDIR = include
TESTDIR = tests
TESTBINDIR = testbin


# INCLUDES can be relative paths, and will not be exported to subdirectory makes.
INCLUDES += -I $(INCLUDEDIR)

# Executable C++ programs - none
# C++ subroutines
SUBS :=  $(wildcard $(SUBDIR)/*.cpp)
SUBOBJS := $(SUBS:$(SUBDIR)/%.cpp=$(OBJDIR)/%.o)

CP = /bin/cp -p
RM = /bin/rm -f

#######################
# Rules - ?? dependencies on INCLUDES ??
#######################

all: $(SUBOBJS)

# Compilation
$(SUBOBJS): $(OBJDIR)/%.o : $(SUBDIR)/%.cpp 
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c $< -o $@

# Remake external packages to update external objects
# Shut off default rule for this; can I do it only when out of date??
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
	rm -rf $(OBJDIR)/*.o $(BINDIR)/* $(TESTBINDIR)/* *~ *.dvi *.aux core .depend

ifeq (.depend, $(wildcard .depend))
include .depend
endif

.PHONY: all install dist depend clean local-clean local-depend exts tests
