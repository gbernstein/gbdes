// Subroutines used by PhotoFit.cpp
#ifndef PHOTOSUBS_H
#define PHOTOSUBS_H

#include <list>
#include <set>
#include "StringStuff.h"

// Load the classes needed
#include "PhotoMapCollection.h"
#include "PhotoMatch.h"
#include "Instrument.h"
#include "FitSubroutines.h"

typedef Photo::Extension Extension;
typedef Photo::ColorExtension ColorExtension;

using namespace photometry;

// Function to produce a list of PhotoPriors from an input file
list<PhotoPrior*>
readPriors(string filename, 
	   const vector<Instrument*>& instruments, 
	   const vector<Exposure*>& exposures, 
	   const vector<Extension*>& extensions);
//???	   const vector<long>& detectionsPerExposure);

#endif
