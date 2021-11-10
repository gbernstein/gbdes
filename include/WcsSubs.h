// Subroutines for use by WCSFit.cpp
#ifndef WCSSUBS_H
#define WCSSUBS_H

#include <list>
#include <set>
#include "StringStuff.h"
#include "NameIndex.h"
#include "Astrometry.h"
#include "FitsTable.h"

// Load the classes needed
#include "FitSubroutines.h"
#include "PixelMapCollection.h"
#include "Match.h"
#include "Instrument.h"
typedef Astro::Extension Extension;
typedef Astro::ColorExtension ColorExtension;

// Function that will using starting WCS to fit all of the defaulted
// maps used by the selected extensions.  Then will put the
// initialized parameters back into the PMC and clear the defaulted flag.
void fitDefaulted(astrometry::PixelMapCollection& pmc,
		  set<Extension*> useThese,
		  const vector<unique_ptr<Instrument>>& instruments,
		  const vector<Exposure*>& exposures,
		  bool logging=true);

// Define and issue WCS for each extension in use, and set projection to
// field coordinates.
void setupWCS(const vector<unique_ptr<astrometry::SphericalCoords>>& fieldProjections,
	      const vector<unique_ptr<Instrument>>& instruments,
	      const vector<Exposure*>& exposures,
	      vector<Extension*>& extensions,
	      astrometry::PixelMapCollection& pmc);

// Analyze the PixelMap to find list of exposures that we can
// initialize first to set up all defaulted device maps.
list<int>
pickExposuresToInitialize(const vector<unique_ptr<Instrument>>& instruments,
			  const vector<Exposure*>& exposures,
			  const vector<Extension*>& extensions,
			  astrometry::PixelMapCollection& pmc);

#endif
