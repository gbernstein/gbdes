// Code for fitting parameters of defaulted pixel maps to the staring WCS solution.
#include "Std.h"
#include "PixelMapCollection.h"
#include "Instrument.h"
#include "Match.h"
#include "WcsSubs.h"
#include "FitSubroutines.h"

using namespace astrometry;

void
fitDefaulted(PixelMapCollection& pmc,
	     set<Extension*> extensions,
	     const vector<Instrument*>& instruments,
	     const vector<Exposure*>& exposures) {

  // Make a new pixel map collection that will hold only the maps
  // involved in this fit.
  PixelMapCollection pmcFit;

  // Take all wcs's used by these extensions and copy them into
  // pmcFit
  for (auto extnptr : extensions) {
    auto pm = pmc.cloneMap(extnptr->mapName);
    pmcFit.learnMap(*pm);
    delete pm;
  }

  // Find all the atomic map components that are defaulted.
  // Fix the parameters of all the others
  set<string> defaultedAtoms;
  set<string> fixAtoms;
  for (auto mapname : pmcFit.allMapNames()) {
    if (!pmc.isAtomic(mapname))
      continue;
    if (pmc.getDefaulted(mapname))
      defaultedAtoms.insert(mapname);
    else
      fixAtoms.insert(mapname);
  }
  // We are done if nothing is defaults
  if (defaultedAtoms.empty())
    return;
  
  // If the exposure has any defaulted maps, initialize them
  /**/if (true) {
    cerr << "Initializing maps " ;
    for (auto name : defaultedAtoms) cerr << name << " ";
    cerr << endl;
  }
  pmcFit.learnMap(IdentityMap());  // Make sure we know this one
  pmcFit.setFixed(fixAtoms);
  pmcFit.rebuildParameterVector();

  auto identityMap = pmcFit.issueMap(IdentityMap().getName());
      
  // Make Matches at a grid of points on each extension's device,
  // matching pix coords to the coordinates derived from startWCS.
  list<Match*> matches;
  for (auto extnptr : extensions) {
    const Exposure& expo = *exposures[extnptr->exposure];
    // Get projection used for this extension, and set up
    // startWCS to reproject into this system.
    extnptr->startWcs->reprojectTo(*expo.projection);

    // Get a realization of the extension's map
    auto map = pmcFit.issueMap(extnptr->mapName);
    
    // Get the boundaries of the device it uses
    Bounds<double> b=instruments[expo.instrument]->domains[extnptr->device];

    // Generate a grid of matched Detections
    const int nGridPoints=512;	// Number of test points for map initialization

    // "errors" on world coords of test points 
    const double testPointSigma = 0.01*ARCSEC/DEGREE;
    const double fitWeight = pow(testPointSigma,-2.);
    // Put smaller errors on the "reference" points.  Doesn't really matter.
    const double refWeight = 10. * fitWeight;
    
    // Distribute points equally in x and y, but shuffle the y coords
    // so that the points fill the rectangle
    vector<int> vx(nGridPoints);
    for (int i=0; i<vx.size(); i++) vx[i]=i;
    vector<int> vy = vx;
    std::random_shuffle(vy.begin(), vy.end());
    double xstep = (b.getXMax()-b.getXMin())/nGridPoints;
    double ystep = (b.getYMax()-b.getYMin())/nGridPoints;
    for (int i=0; i<vx.size(); i++) {
      double xpix = b.getXMin() + (vx[i]+0.5)*xstep;
      double ypix = b.getXMin() + (vy[i]+0.5)*ystep;
      double xw, yw;
      extnptr->startWcs->toWorld(xpix, ypix, xw, yw);
      Detection* dfit = new Detection;
      Detection* dref = new Detection;
      dfit->xpix = xpix;
      dfit->ypix = ypix;
      dref->xpix = xw;
      dref->ypix = yw;
      dref->xw = xw;
      dref->yw = yw;
      
      map->toWorld(xpix,ypix,xw,yw);
      dfit->xw = xw;
      dfit->yw = yw;

      // Set up errors and maps for these matched "detections"
      dfit->wtx = fitWeight;
      dfit->wty = fitWeight;
      dref->wtx = refWeight;
      dref->wty = refWeight;
      dref->map = identityMap;
      dfit->map = map;
      
      // ??? dfit->color = 0.;
      // ??? dref->color = 0.; ??Problem if reference color is not zero?
      matches.push_back(new Match(dfit));
      matches.back()->add(dref);
    }

  }

  // Build CoordAlign object and solve for defaulted parameters
  CoordAlign ca(pmcFit, matches);
  ca.setRelTolerance(0.01);  // weaker tolerance for fit convergence
  ca.fitOnce(false);

  // Copy defaulted parameters back into the parent pmc.
  for (auto mapname : defaultedAtoms) {
    auto pm = pmcFit.cloneMap(mapname);
    pmc.copyParamsFrom(*pm);
    delete pm;
  }

  // Delete the Matches and Detections
  for (auto m : matches) {
    m->clear(true);  // Flush the detections
    // And get rid of match itself.
    delete m;
  } 

  return;
}


// Define and issue WCS for each extension in use, and set projection to
// field coordinates.
void setupWCS(const vector<SphericalCoords*>& fieldProjections,
	      const vector<Instrument*>& instruments,
	      const vector<Exposure*>& exposures,
	      vector<Extension*>& extensions,
	      PixelMapCollection& pmc) {
  for (auto extnptr : extensions) {
    if (!extnptr) continue; // Not in use
    Exposure& expo = *exposures[extnptr->exposure];
    int ifield = expo.field;
    if ( expo.instrument < 0) {
      // Tag & reference exposures have no instruments and no fitting
      // being done.  Coordinates are fixed to xpix = xw.
      // Build a simple Wcs that takes its name from the field
      SphericalICRS icrs;
      if (!pmc.wcsExists(extnptr->wcsName)) {
	// If we are the first reference/tag exposure in this field:
	pmc.defineWcs(extnptr->wcsName, icrs, 
				extnptr->mapName,
				DEGREE);
	auto wcs = pmc.issueWcs(extnptr->wcsName);
	// And have this Wcs reproject into field coordinates, learn as map
	wcs->reprojectTo(*fieldProjections[ifield]);
	pmc.learnMap(*wcs, false, false);
      }
      extnptr->wcs = pmc.issueWcs(extnptr->wcsName);
      extnptr->map = pmc.issueMap(extnptr->wcsName);
    } else {
      // Real instrument, WCS goes into its exposure coordinates
      pmc.defineWcs(extnptr->wcsName, *expo.projection,
			      extnptr->mapName,
			      DEGREE);
      extnptr->wcs = pmc.issueWcs(extnptr->wcsName);

      // Reproject this Wcs into the field system and get a SubMap including reprojection:
      extnptr->wcs->reprojectTo(*fieldProjections[ifield]);
      pmc.learnMap(*extnptr->wcs, false, false);
      extnptr->map = pmc.issueMap(extnptr->wcsName);
    }
  } // end extension loop
}

list<int>
pickExposuresToInitialize(const vector<Instrument*>& instruments,
			  const vector<Exposure*>& exposures,
			  const vector<Extension*>& extensions,
			  PixelMapCollection& pmc) {
  list<int> exposuresToInitialize;
  for (int iInst=0; iInst < instruments.size(); iInst++) {
    if (!instruments[iInst]) continue; // Not in use
    auto& instr = *instruments[iInst];
    // Classify the device maps for this instrument
    set<int> defaultedDevices;
    set<int> initializedDevices;
    set<int> unusedDevices;

    // And the exposure maps as well:
    set<int> defaultedExposures;
    set<int> initializedExposures;
    set<int> unusedExposures;
    set<int> itsExposures;  // All exposure numbers using this instrument
      
    for (int iDev=0; iDev<instr.nDevices; iDev++) {
      string mapName = instr.mapNames[iDev];
      if (!pmc.mapExists(mapName)) {
	unusedDevices.insert(iDev);
	continue;
      }
      if (pmc.getDefaulted(mapName)) {
	defaultedDevices.insert(iDev);
      } else {
	initializedDevices.insert(iDev);
      }
    }
	

    for (int iExpo=0; iExpo < exposures.size(); iExpo++) {
      if (!exposures[iExpo]) continue; // Not in use
      if (exposures[iExpo]->instrument==iInst) {
	itsExposures.insert(iExpo);
	string mapName = exposures[iExpo]->name;
	if (!pmc.mapExists(mapName)) {
	  unusedExposures.insert(iExpo);
	} else if (pmc.getDefaulted(mapName)) {
	  defaultedExposures.insert(iExpo);
	} else {
	  initializedExposures.insert(iExpo);
	}
      }
    }

    /**/cerr << "Done collecting initialized/defaulted" << endl;

    // No need for any of this if there are no defaulted devices
    if (defaultedDevices.empty())
      continue;
      
    // Now take an inventory of all extensions to see which device
    // solutions are used in coordination with which exposure solutions
    vector<set<int>> exposuresUsingDevice(instr.nDevices);
    for (auto extnptr : extensions) {
      if (!extnptr) continue; // Not in use
      int iExpo = extnptr->exposure;
      if (itsExposures.count(iExpo)>0) {
	// Extension is from one of the instrument's exposures
	int iDev = extnptr->device;
	if (pmc.dependsOn(extnptr->mapName, exposures[iExpo]->name)
	    && pmc.dependsOn(extnptr->mapName, instr.mapNames[iDev])) {
	  // If this extension's map uses both the exposure and device
	  // maps, then enter this dependence into our sets
	  exposuresUsingDevice[iDev].insert(iExpo);
	  if (unusedExposures.count(iExpo)>0 ||
	      unusedDevices.count(iDev)>0) {
	    cerr << "Logic problem: extension map "
		 << extnptr->mapName
		 << " is using allegedly unused exposure or device map"
		 << endl;
	    exit(1);
	  }
	}
      }
    }

    /**/cerr << "Done building exposure/device graph" << endl;

    // Find a non-defaulted exposure using all of the defaulted devices
    int exposureForInitializing = -1;
    for (auto iExpo : initializedExposures) {
      bool hasAllDevices = true;
      for (auto iDev : defaultedDevices) {
	if (exposuresUsingDevice[iDev].count(iExpo)==0) {
	  hasAllDevices = false;
	  break;
	}
      }
      if (hasAllDevices) {
	exposureForInitializing = iExpo;
	break;
      }
    }

    if (exposureForInitializing < 0) {
      cerr << "Could not find an exposure that can initialize defaulted devices \n"
	   << "for instrument " << instr.name
	   << endl;
      cerr << "Write more code if you want to exploit more complex situations \n"
	   << "where some non-defaulted devices can initialize a defaulted exposure."
	   << endl;
      exit(1);
    } else {
      exposuresToInitialize.push_back(exposureForInitializing);
      cerr << "Using exposure " << exposures[exposureForInitializing]->name
	   << " to initialize defaulted devices on instrument " << instr.name
	   << endl;
    }
  } // end instrument loop

  return exposuresToInitialize;
}
