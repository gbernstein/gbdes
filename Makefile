# First are site-dependent items - expect these to be defined in environment!
#CXX = g++-4 -fopenmp
#TMV_DIR
#CFITSIO_DIR
#FFTW_DIR
#BOOST_DIR
# GBTOOLS_DIR (defaults to ../gbtools if not in environment)

# GBTOOLS_DIR = ../gbtools

export CXX
export TMV_DIR
export CFITSIO_DIR
export BOOST_DIR
export FFTW_DIR
export YAML_DIR

# OPTFLAGS will be exported for subdir makes
OPTFLAGS = -O3 -DASSERT -DUSE_YAML
export OPTFLAGS


# Paths to the parts of the GBTOOLS we will use here
UTILITIES := $(GBTOOLS_DIR)/utilities
IMAGES := $(GBTOOLS_DIR)/images
ASTROMETRY := $(GBTOOLS_DIR)/astrometry
PHOTOMETRY := $(GBTOOLS_DIR)/photometry

# ABS_INCLUDES are absolute paths 
ABS_INCLUDES = -I $(TMV_DIR)/include -I $(CFITSIO_DIR)/include \
	-I $(FFTW_DIR)/include -I $(BOOST_DIR)/include -I $(YAML_DIR)/include \
	-I $(UTILITIES) -I $(IMAGES) -I $(ASTROMETRY) -I $(PHOTOMETRY)

LIB_DIRS = -L $(CFITSIO_DIR)/lib -L $(TMV_DIR)/lib -L $(FFTW_DIR)/lib \
	-L $(BOOST_DIR)/lib -L $(YAML_DIR)/lib

TMV_LINK := $(shell cat $(TMV_DIR)/share/tmv/tmv-link)

##### Below here should be site-independent
SUBDIRS = $(GBTOOLS_DIR)

# INCLUDES can be relative paths, and will not be exported to subdirectory makes.
INCLUDES = 

CXXFLAGS = $(OPTFLAGS) $(ABS_INCLUDES) $(INCLUDES)
LIBS = -lm $(LIB_DIRS) -lyaml-cpp -lboost_regex -lfftw3 -lcfitsio -ltmv_symband $(TMV_LINK)

SUBOBJ =$(UTILITIES)/BinomFact.o $(UTILITIES)/StringStuff.o $(UTILITIES)/Interpolant.o \
	$(UTILITIES)/fft.o $(UTILITIES)/Table.o $(UTILITIES)/Pset.o $(UTILITIES)/Poly2d.o \
	$(UTILITIES)/Expressions.o $(UTILITIES)/Lookup1d.o \
	$(IMAGES)/FITS.o $(IMAGES)/Header.o $(IMAGES)/Hdu.o $(IMAGES)/FitsTable.o \
	$(IMAGES)/FTable.o $(IMAGES)/FTableExpression.o \
	$(IMAGES)/Image.o $(IMAGES)/FitsImage.o \
	$(IMAGES)/HeaderFromStream.o \
	$(ASTROMETRY)/PixelMap.o $(ASTROMETRY)/Astrometry.o $(ASTROMETRY)/PolyMap.o \
	$(ASTROMETRY)/SubMap.o $(ASTROMETRY)/Wcs.o $(ASTROMETRY)/PixelMapCollection.o \
	$(ASTROMETRY)/TemplateMap.o $(ASTROMETRY)/PiecewiseMap.o \
	$(PHOTOMETRY)/PhotoMap.o $(PHOTOMETRY)/PhotoMapCollection.o $(PHOTOMETRY)/SubMap.o \
	$(PHOTOMETRY)/PhotoMatch.o $(PHOTOMETRY)/PhotoPrior.o $(PHOTOMETRY)/PhotoTemplate.o \
	$(PHOTOMETRY)/PhotoPiecewise.o

SRC = $(shell ls *.cpp)

OBJ =  TPVMap.o FitSubroutines.o $(SUBOBJ)

all: depend subs

bfcorrect: bfcorrect.o $(OBJ)
	$(CXX) $(CXXFLAGS) $^  $(LIBS) -o $@
Photo2DESDM: Photo2DESDM.o DECamInfo.o $(OBJ)
	$(CXX) $(CXXFLAGS) $^  $(LIBS) -o $@
CorrectFlat: CorrectFlat.o DECamInfo.o $(OBJ)
	$(CXX) $(CXXFLAGS) $^  $(LIBS) -o $@
DrawFlat: DrawFlat.o DECamInfo.o $(OBJ)
	$(CXX) $(CXXFLAGS) $^  $(LIBS) -o $@
DrawFlat2: DrawFlat2.o DECamInfo.o $(OBJ)
	$(CXX) $(CXXFLAGS) $^  $(LIBS) -o $@
DrawPhoto: DrawPhoto.o DECamInfo.o $(OBJ)
	$(CXX) $(CXXFLAGS) $^  $(LIBS) -o $@
DrawRegnault: DrawRegnault.o $(OBJ)
	$(CXX) $(CXXFLAGS) $^  $(LIBS) -o $@
ApplyPhoto: ApplyPhoto.o $(OBJ)
	$(CXX) $(CXXFLAGS) $^  $(LIBS) -o $@
ApplyWCS: ApplyWCS.o $(OBJ)
	$(CXX) $(CXXFLAGS) $^  $(LIBS) -o $@
TabulatePixmap: TabulatePixmap.o DECamInfo.o $(OBJ)
	$(CXX) $(CXXFLAGS) $^  $(LIBS) -o $@
UnFlattenCatalog: UnFlattenCatalog.o DECamInfo.o $(OBJ)
	$(CXX) $(CXXFLAGS) $^  $(LIBS) -o $@
UpdateHeaders: UpdateHeaders.o $(OBJ)
	$(CXX) $(CXXFLAGS) $^  $(LIBS) -o $@
FitsGlue: FitsGlue.o $(OBJ)
	$(CXX) $(CXXFLAGS) $^  $(LIBS) -o $@
WCSFoF: WCSFoF.o ExtensionAttribute.o $(OBJ)
	$(CXX) $(CXXFLAGS) $^  $(LIBS) -o $@
WCSFit: WCSFit.o Match.o MapFit.o $(OBJ)
	$(CXX) $(CXXFLAGS) $^  $(LIBS) -o $@
PhotoFit: PhotoFit.o ReadPhotoPriors.o $(OBJ)
	$(CXX) $(CXXFLAGS) $^  $(LIBS) -o $@
MagColor: MagColor.o $(OBJ)
	$(CXX) $(CXXFLAGS) $^  $(LIBS) -o $@
DumpLDACHeader: DumpLDACHeader.o $(OBJ)
	$(CXX) $(CXXFLAGS) $^  $(LIBS) -o $@
ListClipped: ListClipped.o $(OBJ)
	$(CXX) $(CXXFLAGS) $^  $(LIBS) -o $@
#
testYAML:  testYAML.o $(OBJ)
	$(CXX) $(CXXFLAGS) $^  $(LIBS) -o $@
AndresYAML1:  AndresYAML1.o $(OBJ)
	$(CXX) $(CXXFLAGS) $^  $(LIBS) -o $@

testDecam: testDecam.o DECamInfo.o $(OBJ)
	$(CXX) $(CXXFLAGS) $^  $(LIBS) -o $@
testSerialize: testSerialize.o $(OBJ)
	$(CXX) $(CXXFLAGS) $^  $(LIBS) -o $@
testTemplateMap: testTemplateMap.o $(OBJ)
	$(CXX) $(CXXFLAGS) $^  $(LIBS) -o $@
testPrep: testPrep.o $(OBJ)
	$(CXX) $(CXXFLAGS) $^  $(LIBS) -o $@
test: test.o $(OBJ)
	$(CXX) $(CXXFLAGS) $^  $(LIBS) -o $@
testFT: testFT.o $(OBJ)
	$(CXX) $(CXXFLAGS) $^  $(LIBS) -o $@
testFT2: testFT2.o $(OBJ)
	$(CXX) $(CXXFLAGS) $^  $(LIBS) -o $@
testFT3: testFT3.o $(OBJ)
	$(CXX) $(CXXFLAGS) $^  $(LIBS) -o $@
testFT4: testFT4.o $(OBJ)
	$(CXX) $(CXXFLAGS) $^  $(LIBS) -o $@
testFT5: testFT5.o $(OBJ)
	$(CXX) $(CXXFLAGS) $^  $(LIBS) -o $@
testFT6: testFT6.o $(OBJ)
	$(CXX) $(CXXFLAGS) $^  $(LIBS) -o $@
testFT7: testFT7.o $(OBJ)
	$(CXX) $(CXXFLAGS) $^  $(LIBS) -o $@
testFT8: testFT8.o $(OBJ)
	$(CXX) $(CXXFLAGS) $^  $(LIBS) -o $@
testImage: testImage.o $(OBJ) 
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
	rm -f *.o *~ *.dvi *.aux core .depend

ifeq (.depend, $(wildcard .depend))
include .depend
endif

.PHONY: all install dist depend clean 
