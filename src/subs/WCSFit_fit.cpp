// Program to fit astrometric solutions to detection catalogs already matched by WCSFoF.

#include "WCSFit_fit.h"

#define PROGRESS(val, msg) if (verbose>=val) cerr << "-->" <<  #msg << endl

//FitClass::FitClass(string inputMaps) {
FitClass::FitClass() {
    
    int test = 1;
    minMatches = 2;
    minFitExposures = 200;
    maxError = 100.0;
    reserveFraction = 0.2;
    randomNumberSeed = 0;
    minimumImprovement = 0.02;
    clipThresh = 5.;
    fixMapList = {};
    fieldNames = {};
    fieldProjections = vector<SphericalCoords*>(0);
    fieldEpochs = vector<double>(0);
    instruments = vector<Instrument*>(0);
    
    /*cerr << "setup inputYAML" << endl;
    YAMLCollector inputYAML(inputMaps, PixelMapCollection::magicKey);
    //inputYAML = tmpYAML;

    // Make sure inputYAML knows about the Identity transformation:
    {
      istringstream iss("Identity:\n  Type:  Identity\n");
      inputYAML.addInput(iss);
    }*/
    
    // Tried setting this but not accessible. WHY??
    //set<string> degenerateTypes = {"Poly","Linear","Constant"};
    
    //UNCOMMENTING THIS LINE CAUSES SOME ISSUE WITH CONFUSION BETWEEN VARIABLES?
    //PixelMapCollection mapCollection;

    //PixelMapCollection* pmcInit = new PixelMapCollection;
    //PixelMapCollection mapCollection = PixelMapCollection(;
    // Teach PixelMapCollection about new kinds of PixelMaps:
    loadPixelMapParser();
    
}

/*void FitClass::addInputYAML(string inputMaps) {
  inputYAML = astrometry::YAMLCollector(inputMaps, PixelMapCollection::magicKey);

  // Make sure inputYAML knows about the Identity transformation:
  {
    istringstream iss("Identity:\n  Type:  Identity\n");
    inputYAML.addInput(iss);
  }
}*/

void FitClass::setExposures(vector<shared_ptr<Exposure>> expos, double sysErr, double refSysErr) {

  SPexposures = expos;
  //exposures.reserve(SPexposures.size());
  vector<Exposure*> tmp_expos;
  for (auto e : SPexposures) {
    //exposures.push_back(e.get());
    tmp_expos.push_back(e.get());
  }
  setExposures(tmp_expos, sysErr, refSysErr);
}

void FitClass::setExposures(vector<Exposure*> expos, double sysErr, double refSysErr) {

  /*SPexposures = expos;
  exposures.reserve(SPexposures.size());
  for (auto e : SPexposures) {
    exposures.push_back(e.get());
  }*/
  exposures = expos;

  // Convert error parameters from I/O units to internal
  float referenceSysError = refSysErr * RESIDUAL_UNIT/WCS_UNIT;
  float sysError = sysErr * RESIDUAL_UNIT/WCS_UNIT;


  if (sysError > 0.) {
    cerr << "sysError: " << to_string(sysError) << endl;
    Matrix22 astrometricCovariance(0.);
    astrometricCovariance(0,0) = sysError*sysError;
    astrometricCovariance(1,1) = sysError*sysError;
    for (auto e : exposures) {
      if (e && e->instrument >= 0) {
      e->astrometricCovariance += astrometricCovariance;
      cerr << "add to e" << endl;
      }
    }
  }
  if (referenceSysError > 0.) {
    cerr << "rsysError: " << to_string(referenceSysError) << endl;
    Matrix22 astrometricCovariance(0.);
    astrometricCovariance(0,0) = referenceSysError*referenceSysError;
    astrometricCovariance(1,1) = referenceSysError*referenceSysError;
    for (auto e : exposures) {
      if (e && (e->instrument == REF_INSTRUMENT || e->instrument== PM_INSTRUMENT))
        e->astrometricCovariance += astrometricCovariance;
      // Note that code in FitSubroutines::makePMDetection() will
      // add this systematic error only to the Detection::invCov 2d
      // covariance, not the full 5d PMDetection::pmInvCov calculation,
      // since we'll assume these 5d projects (Gaia!) have treated errors well.
      // For single-epoch reference catalogs, PM is probably the largest sys error.
    }
  }
}

void FitClass::setExtensions(vector<shared_ptr<Extension>> extens) {

  SPextensions = extens;
  extensions.reserve(SPextensions.size());
  for (auto e : SPextensions) {
    //extensions.push_back(e.get());
    typename Astro::Extension* extn = new typename Astro::Extension;

    extn->exposure = e->exposure;
    extn->device = e->device;
    extn->airmass = e->airmass;
    extn->magshift = e->magshift;
    extn->wcsName = e->wcsName;
    extn->mapName = e->mapName;
    extn->startWcs = e->startWcs;
    //test reproject:
    //if (extn->device > 0) {
    const Exposure& expo = *exposures[extn->exposure];
    //astrometry::SphericalICRS testCoords(2.62, 0.039);
    extn->startWcs->reprojectTo(*expo.projection);
    //}
    //extn->startWcs->reprojectTo(testCoords);
    extensions.push_back(extn);
  }

}

void FitClass::setRefWCSNames() {
  PROGRESS(2,Setting reference wcsNames);

  // A special loop here to set the wcsname of reference extensions to the
  // name of their field.
  for (auto extnptr : extensions) {
    if (!extnptr) continue;
    const Exposure& expo = *exposures[extnptr->exposure];
    if ( expo.instrument >= 0) continue;
    int ifield = expo.field;
    extnptr->wcsName = fieldNames.nameOf(ifield); // ??? mapName instead???
  }
}

void FitClass::addMap(YAMLCollector& inputYAML, string mapName, vector<string> mapParams) {
//void FitClass::addMap(string mapName, vector<string> mapParams) {
  astrometry::YAMLCollector::Dictionary d;
  d["INSTRUMENT"] = mapParams[0];
  d["DEVICE"] = mapParams[1];
  d["EXPOSURE"] = mapParams[2];
  d["BAND"] = mapParams[3];
  inputYAML.addMap(mapName,d);
}

void FitClass::setupMaps(YAMLCollector& inputYAML) {//, PixelMapCollection& mapCollection) {
//void FitClass::setupMaps() {

    // Make sure inputYAML knows about the Identity transformation:
    //{
    //  istringstream iss("Identity:\n  Type:  Identity\n");
    //  inputYAML.addInput(iss);
    //}
    

      /////////////////////////////////////////////////////
    //  Create and initialize all maps
    /////////////////////////////////////////////////////
    
    PROGRESS(1,Building initial PixelMapCollection);
    cerr << "check 4" << endl;
    // Now build a preliminary set of pixel maps from the configured YAML files
    PixelMapCollection* pmcInit = new PixelMapCollection;
    assert (pmcInit);
    cerr << "pmcInit maps: " << to_string(pmcInit->nMaps()) << endl;
    //cerr << "pmcInit len: " << to_string(pmcInit->allMapNames().size()) << endl;
    pmcInit->learnMap(IdentityMap(), false, false);
    cerr << "check 4.01" << endl;
    cerr << "pmcInit len: " << to_string(pmcInit->allMapNames().size()) << endl;
    {
      istringstream iss(inputYAML.dump());
      //string tmp;
      //iss >> tmp;
      //cerr << iss.str()  << endl;
      cerr << " check 4.02" << endl;
      if (!pmcInit->read(iss)) {
        cerr << "Failure parsing the initial YAML map specs" << endl;
        exit(1);
      }
    }
    cerr << "read done" << endl;
    cerr << "pmcInit len: " << to_string(pmcInit->allMapNames().size()) << endl;
    cerr << "check 4.1" << endl;
    // Check every map name against the list of those to fix.
    // Includes all devices of any instruments on the fixMapList.
    fixMapComponents<Astro>(*pmcInit,
                            fixMapList,
                            instruments);
    cerr << "pmcInit len: " << to_string(pmcInit->allMapNames().size()) << endl;
    cerr << "check 4.2" << endl;
    /////////////////////////////////////////////////////
    // First sanity check: Check for maps that are defaulted but fixed
    /////////////////////////////////////////////////////
    PROGRESS(2,Checking for fixed and defaulted maps);
    {
      bool done = false;
      for (auto iName : pmcInit->allMapNames()) {
        cerr << "iName: " << iName << " " << to_string(pmcInit->getFixed(iName)) << endl;
        if (pmcInit->getFixed(iName) && pmcInit->getDefaulted(iName)) {
          cerr << "ERROR: Map element " << iName
            << " is frozen at defaulted parameters"
            << endl;
          done = true;
	    }
      }
      if (done) exit(1);
    }
    cerr << "check 4.3" << endl;
    /////////////////////////////////////////////////////
    // Second sanity check: in each field that has *any* free maps,
    // do we have at least one map that is *fixed*?
    /////////////////////////////////////////////////////

    PROGRESS(2,Checking for field degeneracy);

    {
      vector<bool> fieldHasFree(fieldNames.size(), false);
      vector<bool> fieldHasFixed(fieldNames.size(), false);
      for (auto extnptr : extensions) {
        if (!extnptr) continue; // Not in use
        int field = exposures[extnptr->exposure]->field;
        if (pmcInit->getFixed(extnptr->mapName)) {
          fieldHasFixed[field] = true;
        }
        else {
          fieldHasFree[field] = true;
        }
      }
      bool done = false;
      for (int i=0; i<fieldHasFree.size(); ++i) {
        if (fieldHasFree[i] && !fieldHasFixed[i]) {
        cerr << "ERROR: No data in field "
          << fieldNames.nameOf(i)
          << " have fixed maps to break shift degeneracy"
          << endl;
        done = true;
        }
      }
      if (done) exit(1);
    }
    cerr << "check 5" << endl;
    /////////////////////////////////////////////////////
    // Next degeneracy: if linear/poly maps are compounded there
    // must be data for which one of them is fixed.  See if it
    // is needed to set some exposure maps to Identity to break
    // such degeneracies.
    /////////////////////////////////////////////////////

    PROGRESS(2,Checking for polynomial degeneracies);
    set<string> degenerateTypes={"Poly","Linear","Constant"};
    {
      MapDegeneracies<Astro> degen(extensions,
                                   *pmcInit,
                                   degenerateTypes,
                                   false);  // Any non-fixed maps are examined
      // All exposure maps are candidates for setting to Identity
      // (the code will ignore those which already are Identity)
      set<string> exposureMapNames;
      for (auto expoPtr : exposures) {
        if (expoPtr && !expoPtr->name.empty())
            exposureMapNames.insert(expoPtr->name);
      }

      auto replaceThese = degen.replaceWithIdentity(exposureMapNames);
      
      // Supersede their maps if there are any
      if (!replaceThese.empty()) {
        PixelMapCollection pmcAltered;
        for (auto mapname : replaceThese) {
          cerr << "..Setting map <" << mapname << "> to Identity" << endl;
          pmcAltered.learnMap(IdentityMap(mapname));
        }
        // Add the altered maps specs to the input YAML specifications
        ostringstream oss;
        pmcAltered.write(oss);
        istringstream iss(oss.str());
        inputYAML.addInput(iss, "", true); // Prepend these specs to others
      }
    } // End of poly-degeneracy check/correction
    
    PROGRESS(1,Making final mapCollection);
    cerr << "check 6" << endl;
    // Do not need the preliminary PMC any more.
    delete pmcInit;
    pmcInit = 0;
    // And clean out any old maps stored in the YAMLCollector
    inputYAML.clearMaps();
    cerr << "check 7 " << endl;

    // Make final map collection from the YAML
    //PixelMapCollection mapCollection;
    createMapCollection<Astro>(instruments,
                               exposures,
                               extensions,
                               inputYAML,
                               mapCollection);


    // Add WCS for every extension, and reproject into field coordinates
    PROGRESS(2,Defining all WCSs);
    cerr << "check 8 wcs" << endl;
    setupWCS(fieldProjections, instruments, exposures, extensions, mapCollection);
    cerr << "check 9 wcs" << endl;

    /////////////////////////////////////////////////////
  //  Initialize map components that were created with default
  //  parameters by fitting to fake data created from the
  //  starting WCS.  
  /////////////////////////////////////////////////////

  PROGRESS(2,Initializing defaulted maps);
  //set<string> degenerateTypes={"Poly","Linear","Constant"};
  cerr << "check 7" << endl;
  // This routine figures out an order in which defaulted maps can
  // be initialized without degeneracies
  list<set<int>> initializeOrder;
  set<int> initializedExtensions;
  {
    MapDegeneracies<Astro> degen(extensions,
                                  mapCollection,
                                  degenerateTypes,
                                  true);  // Only defaulted maps are used
    // Get the recommended initialization order:
    initializeOrder = degen.initializationOrder();
  }
  
  for (auto extnSet : initializeOrder) {
    for (auto s : extnSet) {
        cerr << to_string(s) << " " ;
      }
      cerr << endl;
    // Fit set of extensions to initialize defaulted map(s)
    set<Extension*> defaultedExtensions;
    for (auto iextn : extnSet) {
      defaultedExtensions.insert(extensions[iextn]);
      initializedExtensions.insert(iextn);
    }
    for (auto s : defaultedExtensions) {
        cerr << s->mapName << " " ;
      }
    fitDefaulted(mapCollection,
                  defaultedExtensions,
                  instruments,
                  exposures,
                  verbose>=2);
  }

  cerr << "check 8" << endl;
  // Try to fit on every extension not already initialized just to
  // make sure that we didn't miss any non-Poly map elements.
  // The fitDefaulted routine will just return if there are no
  // defaulted parameters for the extension.
  for (int iextn=0; iextn<extensions.size(); iextn++) {
    // Skip extensions that don't exist or are already initialized
    if (!extensions[iextn] || initializedExtensions.count(iextn))
      continue;
    set<Extension*> defaultedExtensions = {extensions[iextn]};
    fitDefaulted(mapCollection,
                  defaultedExtensions,
                  instruments,
                  exposures,
                  verbose>=2);
    initializedExtensions.insert(iextn);
  }
  cerr << "Initialization order made" << endl;
    
  // As a check, there should be no more defaulted maps
  bool defaultProblem=false;
  for (auto iName : mapCollection.allMapNames()) {
    if (mapCollection.getDefaulted(iName)) {
      cerr << "Logic error: after all intializations, still have map "
        << iName
        << " as defaulted."
        << endl;
      defaultProblem = true;
    }
  }
  if (defaultProblem) exit(1);

  // Now fix all map elements requested to be fixed
  fixMapComponents<Astro>(mapCollection,
                          fixMapList,
                          instruments);
  
  // Recalculate all parameter indices - maps are ready to roll!
  mapCollection.rebuildParameterVector();

  cout << "# Total number of free map elements " << to_string(mapCollection.nFreeMaps())
    << " with " << to_string(mapCollection.nParams()) << " free parameters."
    << endl;
  
    cerr << "who needs color" << endl;
  // Figure out which extensions' maps require a color entry to function
  whoNeedsColor<Astro>(extensions);
  
  PROGRESS(2,Reprojecting startWcs);

  // Before reading objects, we want to set all starting WCS's to go into
  // field coordinates.
  for (auto extnptr : extensions) {
    if (!extnptr) continue;
    if (extnptr->exposure < 0) continue;
    int ifield = exposures[extnptr->exposure]->field;
    extnptr->startWcs->reprojectTo(*fieldProjections[ifield]);
  }
}

void FitClass::defaultMaps() {
  /////////////////////////////////////////////////////
  //  Initialize map components that were created with default
  //  parameters by fitting to fake data created from the
  //  starting WCS.  
  /////////////////////////////////////////////////////

  PROGRESS(2,Initializing defaulted maps);
  set<string> degenerateTypes={"Poly","Linear","Constant"};
  cerr << "check 7" << endl;
  // This routine figures out an order in which defaulted maps can
  // be initialized without degeneracies
  list<set<int>> initializeOrder;
  set<int> initializedExtensions;
  {
    MapDegeneracies<Astro> degen(extensions,
                                  mapCollection,
                                  degenerateTypes,
                                  true);  // Only defaulted maps are used
    // Get the recommended initialization order:
    initializeOrder = degen.initializationOrder();
  }
  
  for (auto extnSet : initializeOrder) {
    for (auto s : extnSet) {
        cerr << s << " " ;
      }
      cerr << endl;
    // Fit set of extensions to initialize defaulted map(s)
    set<Extension*> defaultedExtensions;
    for (auto iextn : extnSet) {
      defaultedExtensions.insert(extensions[iextn]);
      initializedExtensions.insert(iextn);
    }
    for (auto s : defaultedExtensions) {
        cerr << s->mapName << " " ;
      }
    fitDefaulted(mapCollection,
                  defaultedExtensions,
                  instruments,
                  exposures,
                  verbose>=2);
  }

  cerr << "check 8" << endl;
  // Try to fit on every extension not already initialized just to
  // make sure that we didn't miss any non-Poly map elements.
  // The fitDefaulted routine will just return if there are no
  // defaulted parameters for the extension.
  for (int iextn=0; iextn<extensions.size(); iextn++) {
    // Skip extensions that don't exist or are already initialized
    if (!extensions[iextn] || initializedExtensions.count(iextn))
      continue;
    set<Extension*> defaultedExtensions = {extensions[iextn]};
    fitDefaulted(mapCollection,
                  defaultedExtensions,
                  instruments,
                  exposures,
                  verbose>=2);
    initializedExtensions.insert(iextn);
  }
  
    
  // As a check, there should be no more defaulted maps
  bool defaultProblem=false;
  for (auto iName : mapCollection.allMapNames()) {
    if (mapCollection.getDefaulted(iName)) {
      cerr << "Logic error: after all intializations, still have map "
        << iName
        << " as defaulted."
        << endl;
      defaultProblem = true;
    }
  }
  if (defaultProblem) exit(1);

  // Now fix all map elements requested to be fixed
  fixMapComponents<Astro>(mapCollection,
                          fixMapList,
                          instruments);
  
  // Recalculate all parameter indices - maps are ready to roll!
  mapCollection.rebuildParameterVector();

  cout << "# Total number of free map elements " << mapCollection.nFreeMaps()
    << " with " << mapCollection.nParams() << " free parameters."
    << endl;
  
}

void FitClass::reprojectWCSs() {
  cerr << "who needs color" << endl;
  // Figure out which extensions' maps require a color entry to function
  whoNeedsColor<Astro>(extensions);
  
  PROGRESS(2,Reprojecting startWcs);

  // Before reading objects, we want to set all starting WCS's to go into
  // field coordinates.
  for (auto extnptr : extensions) {
    if (!extnptr) continue;
    if (extnptr->exposure < 0) continue;
    int ifield = exposures[extnptr->exposure]->field;
    extnptr->startWcs->reprojectTo(*fieldProjections[ifield]);
  }
}

void FitClass::fit(){

  try {
    cerr << "check 0" << endl;

    int cut=0;
    auto im = matches.begin();
    cerr << "Reviewing matches:" << endl;
    while (im != matches.end() ) {
      if ((*im)->fitSize() < minMatches){
        cerr << (*im)->fitSize() << endl;
        ++cut;
      }
      ++im;
    }
    cerr << "Total to cut: "  << cut << endl;
    PROGRESS(2,Purging defective detections and matches);

    // Get rid of Detections with errors too high
    purgeNoisyDetections<Astro>(maxError, matches, exposures, extensions);
    cerr << "Matches size: " << matches.size() << endl;
    cerr << "check 9" << endl;
    PROGRESS(2,Purging sparse matches);
    // Get rid of Matches with too few detections
    purgeSparseMatches<Astro>(minMatches, matches);
    cerr << "Matches size: " << matches.size() << endl;
    cerr << "check 10" << endl;
    PROGRESS(2,Purging out-of-range colors);
    // Get rid of Matches with color out of range (note that default color is 0).
    purgeBadColor<Astro>(minColor, maxColor, matches);
    cerr << "Matches size: " << matches.size() << endl;
    cerr << "check 11" << endl;
    PROGRESS(2,Reserving matches);
    // Reserve desired fraction of matches
    if (reserveFraction>0.) {
      cerr << "reserve info: " << reserveFraction << "  " << randomNumberSeed << endl;
      reserveMatches<Astro>(matches, reserveFraction, randomNumberSeed);
    }
    cerr << "Matches size: " << matches.size() << endl;
    cerr << "check 12" << endl;
    PROGRESS(2,Purging underpopulated exposures);
    // Find exposures whose parameters are free but have too few
    // Detections being fit to the exposure model.
    auto badExposures = findUnderpopulatedExposures<Astro>(minFitExposures,
							   matches,
							   exposures,
							   extensions,
							   mapCollection);
    cerr << "Matches size: " << matches.size() << endl;
    cerr << "check 13" << endl;
    PROGRESS(2,Purging bad exposure parameters and Detections);
    // Freeze parameters of an exposure model and clip all
    // Detections that were going to use it.
    for (auto i : badExposures) {
      cout << "WARNING: Shutting down exposure map " << i.first
	   << " with only " << i.second
	   << " fitted detections "
	   << endl;
      freezeMap<Astro>(i.first, matches, extensions, mapCollection);
    }
    
    if (purgeOutput) {
      PROGRESS(2,Purging unfittable maps);
      mapCollection.purgeInvalid();
    }
    
    cerr << "check 14" << endl;
    PROGRESS(2,Match census);
    matchCensus<Astro>(matches, cout);

    ///////////////////////////////////////////////////////////
    // Now do the re-fitting 
    ///////////////////////////////////////////////////////////

    PROGRESS(1,Begin fitting process);
    
    // make CoordAlign class
    CoordAlign ca(mapCollection, matches);

    int nclip;
    double oldthresh=0.;

    // Start off in a "coarse" mode so we are not fine-tuning the solution
    // until most of the outliers have been rejected:
    bool coarsePasses = true;
    cerr << "relative Tolerance: " << chisqTolerance << endl;
    ca.setRelTolerance(10.*chisqTolerance);
    // Here is the actual fitting loop 
    do {
      // Report number of active Matches / Detections in each iteration:
      {
        cerr << "Matches size: " << matches.size() << endl;
	long int mcount=0;
	long int dcount=0;
	ca.count(mcount, dcount, false, 2);
	double maxdev=0.;
	int dof=0;
	double chi= ca.chisqDOF(dof, maxdev, false);
	if (verbose>=1)
	  cerr << "Fitting " << mcount << " matches with " << dcount << " detections "
	       << " chisq " << chi << " / " << dof << " dof,  maxdev " << maxdev 
	       << " sigma" << endl;
      }

      // Do the fit here!!
      double chisq = ca.fitOnce(verbose>=1,divideInPlace);  // save space if selected
      // Note that fitOnce() remaps *all* the matches, including reserved ones.
      double max;
      int dof;
      ca.chisqDOF(dof, max, false);	// Exclude reserved Matches
      double thresh = sqrt(chisq/dof) * clipThresh; // ??? change dof to expectedChisq?
      if (verbose>=1)
        cerr << "After iteration: chisq " << chisq 
              << " / " << dof << " dof, max deviation " << max
              << "  new clip threshold at: " << thresh << " sigma"
              << endl;
      if (thresh >= max || (oldthresh>0. && (1-thresh/oldthresh)<minimumImprovement)) {
        cerr << "minImp check: " << thresh << " " << max << " " << oldthresh << " " << (1-thresh/oldthresh) << " " << minimumImprovement << endl;
        // Sigma clipping is no longer doing much.  Quit if we are at full precision,
        // else require full target precision and initiate another pass.
        if (coarsePasses) {
          coarsePasses = false;
          ca.setRelTolerance(chisqTolerance);
          PROGRESS(1,Starting strict tolerance passes);
          if (clipEntireMatch && verbose>=1) cerr << "-->clipping full matches" << endl;
          oldthresh = thresh;
          nclip = ca.sigmaClip(thresh, false, clipEntireMatch && !coarsePasses,
                    verbose>=1);
          if (verbose>=0)
            cerr << "# Clipped " << nclip << " matches " << endl;
          continue;
        } else {
          // Done!
          break;
        }
      }
      oldthresh = thresh;
      // Clip entire matches on final passes if clipEntireMatch=true
      nclip = ca.sigmaClip(thresh, false, clipEntireMatch && !coarsePasses,
			   verbose>=1);
      if (nclip==0 && coarsePasses) {
        // Nothing being clipped; tighten tolerances and re-fit
        coarsePasses = false;
        ca.setRelTolerance(chisqTolerance);
        PROGRESS(1,Starting strict tolerance passes);
        if (clipEntireMatch && verbose>=1) cerr << "-->clipping full matches" << endl;
        continue;
      }
      if (verbose>=0)
	cerr << "# Clipped " << nclip << " matches " << endl;
      
    } while (coarsePasses || nclip>0);
  
    
    // If there are reserved Matches, run sigma-clipping on them now.
    if (reserveFraction > 0.) {
      PROGRESS(1,Clipping reserved matches);
      clipReserved<Astro>(ca, clipThresh, minimumImprovement,
			  false, verbose>=1);  
    }
	
    //////////////////////////////////////
    // Output data and calculate some statistics
    //////////////////////////////////////

    PROGRESS(1,Saving astrometric residuals);
    vector<SphericalCoords*> extensionProjections(extensions.size());
    // Save the pointwise fitting results
    {
      // This routine needs an array of field projections for each extension
      //vector<SphericalCoords*> extensionProjections(extensions.size(),
	  //					nullptr);
      for (int i=0; i<extensions.size(); i++) {
	if (!extensions[i])
	  continue;
	int iExposure = extensions[i]->exposure;
	if (iExposure<0 || !exposures[iExposure])
	  continue;
	int iField = exposures[iExposure]->field;
	extensionProjections[i] = fieldProjections[iField];
      }
      PROGRESS(2, extensionProjections completed);
      
    }
    
  } catch (std::runtime_error& m) {
    quit(m,1);
  }
}

void FitClass::cleanup() {
  //////////////////////////////////////
  // Cleanup:
  //////////////////////////////////////

  PROGRESS(2,Cleaning up);
  // Get rid of matches:
  for (auto im = matches.begin(); im!=matches.end(); ) {
    (*im)->clear(true);  // deletes detections
    // And get rid of match itself.
    im = matches.erase(im);
  }
  // Get rid of the coordinate systems for each field:
  for (int i=0; i<fieldProjections.size(); i++)
    delete fieldProjections[i];
  // Get rid of extensions
  for (int i=0; i<extensions.size(); i++)
    if (extensions[i]) delete extensions[i];
  // Get rid of exposures
  //for (int i=0; i<exposures.size(); i++)
  //  if (exposures[i]) delete exposures[i];
  // Get rid of instruments
  for (int i=0; i<instruments.size(); i++)
    if (instruments[i]) delete instruments[i];
}
