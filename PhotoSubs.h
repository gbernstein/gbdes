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
typedef ExtensionBase<photometry::SubMap, photometry::Detection> Extension;
typedef ColorExtensionBase<photometry::Match> ColorExtension;

using namespace photometry;

// Function to produce a list of PhotoPriors from an input file
list<PhotoPrior*>
readPriors(string filename, 
	   const vector<Instrument*>& instruments, 
	   const vector<Exposure*>& exposures, 
	   const vector<Extension*>& extensions, 
	   const vector<long>& detectionsPerExposure);


// Class to accumulate residual statistics
class Accum {
public:
  double sum_m;
  double sum_mw;
  double sum_mm;
  double sum_mmw;
  double sum_w;
  int n;
  double sum_x;
  double sum_y;
  double chisq;
  Accum();
  void add(const Detection* d, double magoff=0., double dof=1.);
  double rms() const;
  double reducedChisq() const;
  string summary() const;
};


#endif
