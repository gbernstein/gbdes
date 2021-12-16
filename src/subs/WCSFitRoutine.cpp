// Program to fit astrometric solutions to detection catalogs already matched by WCSFoF.

#include "WCSFitRoutine.h"

#define PROGRESS(val, msg) if (verbose>=val) cerr << "-->" <<  #msg << endl

WCSFit::WCSFit(Fields & fields_,
                   vector<shared_ptr<Instrument>> instruments_,
                   ExposuresHelper exposures_,
                   vector<int> extensionExposureNumbers,
                   vector<int> extensionDevices,
                   YAMLCollector inputYAML,
                   vector<shared_ptr<astrometry::Wcs>> wcss,
                   vector<int> sequence,
                   vector<LONGLONG> extns,
                   vector<LONGLONG> objects,
                   double sysErr,
                   double refSysErr,
                   int minMatches,
                   string skipObjectsFile,
                   string fixMaps,
                   bool usePM,
                   int verbose
                   ) : minMatches(minMatches), verbose(verbose), fields(std::move(fields_)) {

  // Set Instruments:
  instruments.reserve(instruments_.size());
  for (int i=0; i < instruments_.size(); i++) {
    shared_ptr<Instrument> inst_init = instruments_[i];
    unique_ptr<Instrument> inst(new Instrument(inst_init->name));
    inst->band = inst_init->band;
    for (int d=0; d < inst_init->nDevices; d++) {
      inst->addDevice(inst_init->deviceNames.nameOf(d), inst_init->domains[d]);
    }
    instruments.emplace_back(std::move(inst));
  }
  
  // Set Exposures:
  vector<unique_ptr<Exposure>> tmpExposures;
  tmpExposures.reserve(exposures_.expNames.size());
  for (int i=0; i < exposures_.expNames.size(); i++) {
    astrometry::Gnomonic gn(astrometry::Orientation(astrometry::SphericalICRS(exposures_.ras[i],
                                                                              exposures_.decs[i])));
    unique_ptr<Exposure> expo(new Exposure(exposures_.expNames[i], gn));
    expo->field = exposures_.fieldNumbers[i];
    expo->instrument = exposures_.instrumentNumbers[i];
    expo->airmass = exposures_.airmasses[i];
    expo->exptime = exposures_.exposureTimes[i];
    expo->mjd = exposures_.mjds[i];
    tmpExposures.emplace_back(std::move(expo));
  }
  setExposures(std::move(tmpExposures), sysErr, refSysErr);
  
  // Set Extensions:
  extensions.reserve(extensionExposureNumbers.size());
  for (int i=0; i < extensionExposureNumbers.size(); i++) {
    unique_ptr<Extension> extn(new Extension());
    int iExposure = extensionExposureNumbers[i];
    extn->exposure = iExposure;
    const Exposure& expo = *exposures[iExposure];
    extn->device = extensionDevices[i];
    shared_ptr<astrometry::Wcs> tmpWcs = wcss[i];
    extn->startWcs = unique_ptr<astrometry::Wcs>(tmpWcs.get());
    extn->startWcs->reprojectTo(*expo.projection);
    if (expo.instrument < 0) {
      // This is the reference catalog:
      extn->mapName = astrometry::IdentityMap().getName();
    }
    else {
      astrometry::YAMLCollector::Dictionary d;
      d["INSTRUMENT"] = instruments[expo.instrument]->name;
      d["DEVICE"] = instruments[expo.instrument]->deviceNames.nameOf(extn->device);
      d["EXPOSURE"] = expo.name;
      d["BAND"] = instruments[expo.instrument]->band;
      // Add Epoch to the dictionary if it is not empty
      if (!expo.epoch.empty())
        d["EPOCH"] = expo.epoch;
      // This will be the name of the extension's final WCS and photometry map:
      extn->wcsName = d["EXPOSURE"] + "/" + d["DEVICE"];
      // And this is the name of the PixelMap underlying the WCS, or
      // the name of the photometric map (for which we do not want the "base" part).
      extn->mapName = extn->wcsName + "/base";
      if (!inputYAML.addMap(extn->mapName,d)) {
        cerr << "Input YAML files do not have complete information for map "
             << extn->mapName
             << endl;
      }
    }
    extensions.emplace_back(std::move(extn));
  }
  colorExtensions = vector<unique_ptr<typename Astro::ColorExtension>>(extensions.size());

  loadPixelMapParser();

  setRefWCSNames();

  setupMaps(inputYAML, fixMaps);

  ExtensionObjectSet matchSkipSet(skipObjectsFile);
  readMatches<Astro>(sequence, extns, objects, matches, extensions, colorExtensions,
              matchSkipSet, minMatches, usePM);
  
}

WCSFit::WCSFit(int minMatches, int verbose): minMatches(minMatches), verbose(verbose){

    // Teach PixelMapCollection about new kinds of PixelMaps:
    loadPixelMapParser();
    
}

void WCSFit::setExposures(vector<unique_ptr<Exposure>> expos, double sysErr, double refSysErr) {

  exposures = std::move(expos);

  // Convert error parameters from I/O units to internal
  float referenceSysError = refSysErr * RESIDUAL_UNIT/WCS_UNIT;
  float sysError = sysErr * RESIDUAL_UNIT/WCS_UNIT;

  if (sysError > 0.) {
    Matrix22 astrometricCovariance(0.);
    astrometricCovariance(0,0) = sysError*sysError;
    astrometricCovariance(1,1) = sysError*sysError;
    for (auto const & e : exposures) {
      if (e && e->instrument >= 0) {
      e->astrometricCovariance += astrometricCovariance;
      }
    }
  }
  if (referenceSysError > 0.) {
    Matrix22 astrometricCovariance(0.);
    astrometricCovariance(0,0) = referenceSysError*referenceSysError;
    astrometricCovariance(1,1) = referenceSysError*referenceSysError;
    for (auto const & e : exposures) {
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

void WCSFit::setRefWCSNames() {
  PROGRESS(2,Setting reference wcsNames);

  // A special loop here to set the wcsname of reference extensions to the
  // name of their field.
  for (auto const & extnptr : extensions) {
    if (!extnptr) continue;
    const Exposure& expo = *exposures[extnptr->exposure];
    if ( expo.instrument >= 0) continue;
    int ifield = expo.field;
    extnptr->wcsName = fields.names().nameOf(ifield);
  }
}

void WCSFit::addMap(YAMLCollector& inputYAML, string mapName, vector<string> mapParams) {

  astrometry::YAMLCollector::Dictionary d;
  d["INSTRUMENT"] = mapParams[0];
  d["DEVICE"] = mapParams[1];
  d["EXPOSURE"] = mapParams[2];
  d["BAND"] = mapParams[3];
  inputYAML.addMap(mapName,d);
}

void WCSFit::setupMaps(YAMLCollector& inputYAML, string fixMaps) {

  /////////////////////////////////////////////////////
  //  Create and initialize all maps
  /////////////////////////////////////////////////////

  PROGRESS(1,Building initial PixelMapCollection);

  // Now build a preliminary set of pixel maps from the configured YAML files
  PixelMapCollection* pmcInit = new PixelMapCollection;
  pmcInit->learnMap(IdentityMap(), false, false);

  {
    istringstream iss(inputYAML.dump());
    if (!pmcInit->read(iss)) {
      cerr << "Failure parsing the initial YAML map specs" << endl;
      exit(1);
    }
  }

  // Check every map name against the list of those to fix.
  // Includes all devices of any instruments on the fixMapList.
  list<string> fixMapList = splitArgument(fixMaps);
  fixMapComponents<Astro>(*pmcInit,
                          fixMapList,
                          instruments);

  /////////////////////////////////////////////////////
  // First sanity check: Check for maps that are defaulted but fixed
  /////////////////////////////////////////////////////
  PROGRESS(2,Checking for fixed and defaulted maps);
  {
    bool done = false;
    for (auto iName : pmcInit->allMapNames()) {
      if (pmcInit->getFixed(iName) && pmcInit->getDefaulted(iName)) {
        cerr << "ERROR: Map element " << iName
          << " is frozen at defaulted parameters"
          << endl;
        done = true;
      }
    }
    if (done) exit(1);
  }

  /////////////////////////////////////////////////////
  // Second sanity check: in each field that has *any* free maps,
  // do we have at least one map that is *fixed*?
  /////////////////////////////////////////////////////

  PROGRESS(2,Checking for field degeneracy);

  {
    vector<bool> fieldHasFree(fields.names().size(), false);
    vector<bool> fieldHasFixed(fields.names().size(), false);
    for (auto const & extnptr : extensions) {
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
        << fields.names().nameOf(i)
        << " have fixed maps to break shift degeneracy"
        << endl;
      done = true;
      }
    }
    if (done) exit(1);
  }

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
    for (auto const & expoPtr : exposures) {
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

  // Do not need the preliminary PMC any more.
  delete pmcInit;
  pmcInit = 0;
  // And clean out any old maps stored in the YAMLCollector
  inputYAML.clearMaps();

  // Make final map collection from the YAML
  //PixelMapCollection mapCollection;
  createMapCollection<Astro>(instruments,
                              exposures,
                              extensions,
                              inputYAML,
                              mapCollection);

  // Add WCS for every extension, and reproject into field coordinates
  PROGRESS(2,Defining all WCSs);
  setupWCS(fields.projections(), instruments, exposures, extensions, mapCollection);

  /////////////////////////////////////////////////////
  //  Initialize map components that were created with default
  //  parameters by fitting to fake data created from the
  //  starting WCS.  
  /////////////////////////////////////////////////////

  PROGRESS(2,Initializing defaulted maps);
  
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
    // Fit set of extensions to initialize defaulted map(s)
    set<Extension*> defaultedExtensions;
    for (auto iextn : extnSet) {
      defaultedExtensions.insert(extensions[iextn].get());
      initializedExtensions.insert(iextn);
    }
    fitDefaulted(mapCollection,
                  defaultedExtensions,
                  instruments,
                  exposures,
                  verbose>=2);
  }

  // Try to fit on every extension not already initialized just to
  // make sure that we didn't miss any non-Poly map elements.
  // The fitDefaulted routine will just return if there are no
  // defaulted parameters for the extension.
  for (int iextn=0; iextn<extensions.size(); iextn++) {
    // Skip extensions that don't exist or are already initialized
    if (!extensions[iextn] || initializedExtensions.count(iextn))
      continue;
    set<Extension*> defaultedExtensions = {extensions[iextn].get()};
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

  cout << "# Total number of free map elements " << to_string(mapCollection.nFreeMaps())
    << " with " << to_string(mapCollection.nParams()) << " free parameters."
    << endl;
  
  // Figure out which extensions' maps require a color entry to function
  whoNeedsColor<Astro>(extensions);
  
  PROGRESS(2,Reprojecting startWcs);

  // Before reading objects, we want to set all starting WCS's to go into
  // field coordinates.
  for (auto const & extnptr : extensions) {
    if (!extnptr) continue;
    if (extnptr->exposure < 0) continue;
    int ifield = exposures[extnptr->exposure]->field;
    extnptr->startWcs->reprojectTo(*fields.projections()[ifield]);
  }
}

void WCSFit::setMatches(vector<int> sequence, vector<LONGLONG> extns, vector<LONGLONG> objects,
                          ExtensionObjectSet skipSet, bool usePM) {
  readMatches<Astro>(sequence, extns, objects, matches, extensions, colorExtensions,
              skipSet, minMatches, usePM);
}

void WCSFit::setObjects(int i, map<string, vector<double>> tableMap, 
                          string xKey, string yKey,
                          vector<string> xyErrKeys, string idKey, string pmCovKey, string magKey, int magKeyElement, string magErrKey,
                          int magErrKeyElement, string pmRaKey, string pmDecKey, string parallaxKey) {

        img::FTable ff;
        for (auto x : tableMap) {
          ff.addColumn<double>(x.second, x.first);
        }
        readObjects_oneExtension<Astro>(exposures, i, ff, xKey, yKey, idKey, pmCovKey, xyErrKeys,
                                 magKey, magKeyElement, magErrKey, magErrKeyElement, pmRaKey, pmDecKey,
                                 parallaxKey, extensions, fields.projections(), verbose, true);
}

void WCSFit::fit(double maxError, int minFitExposures, double reserveFraction, int randomNumberSeed,
                   double minimumImprovement, double clipThresh, double chisqTolerance, bool clipEntireMatch,
                   bool divideInPlace, bool purgeOutput, double minColor, double maxColor){
  try {
    
    PROGRESS(2,Purging defective detections and matches);

    // Get rid of Detections with errors too high
    maxError *= RESIDUAL_UNIT/WCS_UNIT;
    purgeNoisyDetections<Astro>(maxError, matches, exposures, extensions);
    cerr << "Matches size after purgeNoisyDetections: " << to_string(matches.size()) << endl;

    PROGRESS(2,Purging sparse matches);
    // Get rid of Matches with too few detections
    purgeSparseMatches<Astro>(minMatches, matches);
    cerr << "Matches size after purgeSparseMatches: " << to_string(matches.size()) << endl;
    
    PROGRESS(2,Purging out-of-range colors);
    // Get rid of Matches with color out of range (note that default color is 0).
    purgeBadColor<Astro>(minColor, maxColor, matches);
    cerr << "Matches size after purgeBadColor: " << to_string(matches.size()) << endl;
    
    PROGRESS(2,Reserving matches);
    // Reserve desired fraction of matches
    if (reserveFraction>0.) { 
      reserveMatches<Astro>(matches, reserveFraction, randomNumberSeed);
    }

    PROGRESS(2,Purging underpopulated exposures);
    // Find exposures whose parameters are free but have too few
    // Detections being fit to the exposure model.
    auto badExposures = findUnderpopulatedExposures<Astro>(minFitExposures,
							   matches,
							   exposures,
							   extensions,
							   mapCollection);

    PROGRESS(2,Purging bad exposure parameters and Detections);
    // Freeze parameters of an exposure model and clip all
    // Detections that were going to use it.
    for (auto i : badExposures) {
      cout << "WARNING: Shutting down exposure map " << i.first
	   << " with only " << to_string(i.second)
	   << " fitted detections "
	   << endl;
      freezeMap<Astro>(i.first, matches, extensions, mapCollection);
    }
    
    if (purgeOutput) {
      PROGRESS(2,Purging unfittable maps);
      mapCollection.purgeInvalid();
    }

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

    ca.setRelTolerance(10.*chisqTolerance);
    // Here is the actual fitting loop 
    do {
      // Report number of active Matches / Detections in each iteration:
      {
        cerr << "Matches size: " << to_string(matches.size()) << endl;
        long int mcount=0;
        long int dcount=0;
        ca.count(mcount, dcount, false, 2);
        double maxdev=0.;
        int dof=0;
        double chi= ca.chisqDOF(dof, maxdev, false);
        if (verbose>=1)
          cerr << "Fitting " << to_string(mcount) << " matches with " << to_string(dcount) << " detections "
              << " chisq " << to_string(chi) << " / " << to_string(dof) << " dof,  maxdev " << to_string(maxdev)
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
        cerr << "After iteration: chisq " << to_string(chisq)
              << " / " << to_string(dof) << " dof, max deviation " << to_string(max)
              << "  new clip threshold at: " << to_string(thresh) << " sigma"
              << endl;
      if (thresh >= max || (oldthresh>0. && (1-thresh/oldthresh)<minimumImprovement)) {
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
            cerr << "# Clipped " << to_string(nclip) << " matches " << endl;
          continue;
        } else {
          // Done!
          break;
        }
      }
      oldthresh = thresh;
      // Clip entire matches on final passes if clipEntireMatch=true
      nclip = ca.sigmaClip(thresh, false, clipEntireMatch && !coarsePasses, verbose>=1);
      if (nclip==0 && coarsePasses) {
        // Nothing being clipped; tighten tolerances and re-fit
        coarsePasses = false;
        ca.setRelTolerance(chisqTolerance);
        PROGRESS(1,Starting strict tolerance passes);
        if (clipEntireMatch && verbose>=1) cerr << "-->clipping full matches" << endl;
        continue;
      }
      if (verbose>=0)
        cerr << "# Clipped " << to_string(nclip) << " matches " << endl;
      
    } while (coarsePasses || nclip>0);
  
    
    // If there are reserved Matches, run sigma-clipping on them now.
    if (reserveFraction > 0.) {
      PROGRESS(1,Clipping reserved matches);
      clipReserved<Astro>(ca, clipThresh, minimumImprovement, false, verbose>=1);  
    }

    //////////////////////////////////////
    // Output data and calculate some statistics
    //////////////////////////////////////

    // Report summary of residuals to stdout
    Astro::reportStatistics(matches, exposures, extensions, cout);
  } catch (std::runtime_error& m) {
    quit(m,1);
  }
}

void WCSFit::saveResults(string outWcs, string outCatalog, string starCatalog) {

  // The re-fitting is now complete.  Serialize all the fitted coordinate systems
  PROGRESS(2,Saving astrometric parameters);
  // Save the pointwise fitting results
  {
    ofstream ofs(outWcs.c_str());
    if (!ofs) {
      cerr << "Error trying to open output file for fitted Wcs: "
        << outWcs << endl;
      // *** will not quit before dumping output ***
    } else {
      mapCollection.write(ofs);
    }
  }

  PROGRESS(1,Saving astrometric residuals);
  // Save the pointwise fitting results
  {
    // This routine needs an array of field projections for each extension
    vector<SphericalCoords*> extensionProjections(extensions.size());
    for (int i=0; i<extensions.size(); i++) {
      if (!extensions[i])
        continue;
      int iExposure = extensions[i]->exposure;
      if (iExposure<0 || !exposures[iExposure])
        continue;
      int iField = exposures[iExposure]->field;
      extensionProjections[i] = fields.projections()[iField].get();
    }
    PROGRESS(2, extensionProjections completed);
    Astro::saveResults(matches, outCatalog, starCatalog, extensionProjections);
  }

}