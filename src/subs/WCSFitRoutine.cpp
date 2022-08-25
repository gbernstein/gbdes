// Program to fit astrometric solutions to detection catalogs already matched by WCSFoF.

#include "WCSFitRoutine.h"

#define PROGRESS(val, msg) \
    if (verbose >= val) std::cerr << "-->" << #msg << std::endl

WCSFit::WCSFit(Fields &fields_, std::vector<std::shared_ptr<Instrument>> instruments_,
               ExposuresHelper exposures_, std::vector<int> extensionExposureNumbers,
               std::vector<int> extensionDevices, astrometry::YAMLCollector inputYAML,
               std::vector<std::shared_ptr<astrometry::Wcs>> wcss, const std::vector<int> &sequence,
               const std::vector<LONGLONG> &extns, const std::vector<LONGLONG> &objects,
               const std::vector<int> &exposureColorPriorities, double sysErr, double refSysErr,
               int minMatches, std::string skipObjectsFile, std::string fixMaps, bool usePM, double pmPrior,
               double parallaxPrior, int verbose)
        : minMatches(minMatches), verbose(verbose), fields(std::move(fields_)) {

    if (usePM) astrometry::PMMatch::setPrior(pmPrior, parallaxPrior);

    // Set Instruments:
    instruments.reserve(instruments_.size());
    for (int i = 0; i < instruments_.size(); i++) {
        std::unique_ptr<Instrument> inst(new Instrument(*instruments_[i]));
        instruments.emplace_back(std::move(inst));
    }
    std::cerr << "made instruments" << std::endl;
    setExposures(exposures_.getExposuresVector(), sysErr, refSysErr);
    std::cerr << "made expos" << std::endl;
    // Set Extensions:
    extensions.reserve(extensionExposureNumbers.size());
    colorExtensions = vector<unique_ptr<typename Astro::ColorExtension>>(extensionExposureNumbers.size());
    for (int i = 0; i < extensionExposureNumbers.size(); i++) {
        std::unique_ptr<Extension> extn(new Extension());
        int iExposure = extensionExposureNumbers[i];
        extn->exposure = iExposure;
        const Exposure &expo = *exposures[iExposure];
        extn->device = extensionDevices[i];
        wcss[i]->reprojectTo(*expo.projection);
        extn->startWcs = std::unique_ptr<astrometry::Wcs>(wcss[i]->duplicate());
        if (expo.instrument < 0) {
            // This is the reference catalog:
            extn->mapName = astrometry::IdentityMap().getName();
        } else {
            astrometry::YAMLCollector::Dictionary d;
            d["INSTRUMENT"] = instruments[expo.instrument]->name;
            d["DEVICE"] = instruments[expo.instrument]->deviceNames.nameOf(extn->device);
            d["EXPOSURE"] = expo.name;
            d["BAND"] = instruments[expo.instrument]->band;
            // Add Epoch to the dictionary if it is not empty
            if (!expo.epoch.empty()) d["EPOCH"] = expo.epoch;
            // This will be the name of the extension's final WCS and photometry map:
            extn->wcsName = d["EXPOSURE"] + "/" + d["DEVICE"];
            // And this is the name of the PixelMap underlying the WCS, or
            // the name of the photometric map (for which we do not want the "base" part).
            extn->mapName = extn->wcsName + "/base";
            if (!inputYAML.addMap(extn->mapName, d)) {
                throw std::runtime_error("Input YAML files do not have complete information for map " +
                                         extn->mapName);
            }
        }
        extensions.emplace_back(std::move(extn));

        if (exposureColorPriorities.size() > 0) {
            int colorPriority = exposureColorPriorities[iExposure];
            if (colorPriority >= 0) {
                std::unique_ptr<Astro::ColorExtension> ce(new typename Astro::ColorExtension);
                ce->priority = colorPriority;
                colorExtensions[i] = std::move(ce);
            }
        }
    }
    std::cerr << "made extens" << std::endl;
    loadPixelMapParser();

    setRefWCSNames();
    std::cerr << "start setupMaps" << std::endl;
    setupMaps(inputYAML, fixMaps);
    std::cerr << "setup maps" << std::endl;
    ExtensionObjectSet matchSkipSet(skipObjectsFile);
    readMatches<Astro>(sequence, extns, objects, matches, extensions, colorExtensions, matchSkipSet,
                       minMatches, usePM);
}

WCSFit::WCSFit(int minMatches, int verbose) : minMatches(minMatches), verbose(verbose) {
    // Teach PixelMapCollection about new kinds of PixelMaps:
    loadPixelMapParser();
}

void WCSFit::setExposures(std::vector<std::unique_ptr<Exposure>> expos, double sysErr, double refSysErr) {
    exposures = std::move(expos);

    // Convert error parameters from I/O units to internal
    float referenceSysError = refSysErr * astrometry::RESIDUAL_UNIT / astrometry::WCS_UNIT;
    float sysError = sysErr * astrometry::RESIDUAL_UNIT / astrometry::WCS_UNIT;

    if (sysError > 0.) {
        astrometry::Matrix22 astrometricCovariance(0.);
        astrometricCovariance(0, 0) = sysError * sysError;
        astrometricCovariance(1, 1) = sysError * sysError;
        for (auto const &e : exposures) {
            if (e && e->instrument >= 0) {
                e->astrometricCovariance += astrometricCovariance;
            }
        }
    }
    if (referenceSysError > 0.) {
        astrometry::Matrix22 astrometricCovariance(0.);
        astrometricCovariance(0, 0) = referenceSysError * referenceSysError;
        astrometricCovariance(1, 1) = referenceSysError * referenceSysError;
        for (auto const &e : exposures) {
            if (e && (e->instrument == REF_INSTRUMENT || e->instrument == PM_INSTRUMENT))
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
    PROGRESS(2, Setting reference wcsNames);

    // A special loop here to set the wcsname of reference extensions to the
    // name of their field.
    for (auto const &extnptr : extensions) {
        if (!extnptr) continue;
        const Exposure &expo = *exposures[extnptr->exposure];
        if (expo.instrument >= 0) continue;
        int ifield = expo.field;
        extnptr->wcsName = fields.names().nameOf(ifield);
    }
}

void WCSFit::addMap(astrometry::YAMLCollector &inputYAML, std::string mapName,
                    std::vector<std::string> mapParams) {
    astrometry::YAMLCollector::Dictionary d;
    d["INSTRUMENT"] = mapParams[0];
    d["DEVICE"] = mapParams[1];
    d["EXPOSURE"] = mapParams[2];
    d["BAND"] = mapParams[3];
    inputYAML.addMap(mapName, d);
}

void WCSFit::setupMaps(astrometry::YAMLCollector &inputYAML, std::string fixMaps) {
    /////////////////////////////////////////////////////
    //  Create and initialize all maps
    /////////////////////////////////////////////////////

    PROGRESS(1, Building initial PixelMapCollection);

    // Now build a preliminary set of pixel maps from the configured YAML files
    astrometry::PixelMapCollection *pmcInit = new astrometry::PixelMapCollection;
    pmcInit->learnMap(astrometry::IdentityMap(), false, false);

    {
        std::istringstream iss(inputYAML.dump());
        if (!pmcInit->read(iss)) {
            throw std::runtime_error("Failure parsing the initial YAML map specs");
        }
    }

    // Check every map name against the list of those to fix.
    // Includes all devices of any instruments on the fixMapList.
    std::list<std::string> fixMapList = splitArgument(fixMaps);
    fixMapComponents<Astro>(*pmcInit, fixMapList, instruments);

    /////////////////////////////////////////////////////
    // First sanity check: Check for maps that are defaulted but fixed
    /////////////////////////////////////////////////////
  PROGRESS(2,Checking for fixed and defaulted maps);
  {
      bool done = false;
      for (auto iName : pmcInit->allMapNames()) {
          if (pmcInit->getFixed(iName) && pmcInit->getDefaulted(iName)) {
              std::cerr << "ERROR: Map element " << iName << " is frozen at defaulted parameters"
                        << std::endl;
              done = true;
          }
      }
      if (done) {
          throw std::runtime_error("Some map elements are frozen at defaulted parameters");
      }
  }

  /////////////////////////////////////////////////////
  // Second sanity check: in each field that has *any* free maps,
  // do we have at least one map that is *fixed*?
  /////////////////////////////////////////////////////

  PROGRESS(2,Checking for field degeneracy);

  {
      std::vector<bool> fieldHasFree(fields.names().size(), false);
      std::vector<bool> fieldHasFixed(fields.names().size(), false);
      for (auto const &extnptr : extensions) {
          if (!extnptr) continue;  // Not in use
          int field = exposures[extnptr->exposure]->field;
          if (pmcInit->getFixed(extnptr->mapName)) {
              fieldHasFixed[field] = true;
          } else {
              fieldHasFree[field] = true;
          }
      }
      bool done = false;
      for (int i = 0; i < fieldHasFree.size(); ++i) {
          if (fieldHasFree[i] && !fieldHasFixed[i]) {
              std::cerr << "ERROR: No data in field " << fields.names().nameOf(i)
                        << " have fixed maps to break shift degeneracy" << std::endl;
              done = true;
          }
      }
      if (done) {
          throw std::runtime_error("Cannot break degeneracies");
      }
  }

  /////////////////////////////////////////////////////
  // Next degeneracy: if linear/poly maps are compounded there
  // must be data for which one of them is fixed.  See if it
  // is needed to set some exposure maps to Identity to break
  // such degeneracies.
  /////////////////////////////////////////////////////

  PROGRESS(2,Checking for polynomial degeneracies);
  std::set<std::string> degenerateTypes = {"Poly", "Linear", "Constant"};
  {
      MapDegeneracies<Astro> degen(extensions, *pmcInit, degenerateTypes,
                                   false);  // Any non-fixed maps are examined
      // All exposure maps are candidates for setting to Identity
      // (the code will ignore those which already are Identity)
      std::set<std::string> exposureMapNames;
      for (auto const &expoPtr : exposures) {
          if (expoPtr && !expoPtr->name.empty()) {
              // Check if the exposure maps are composite, and if so add the component elements
              if (pmcInit->mapExists(expoPtr->name)) {
                  std::set<string> dependentElements = pmcInit->dependencies(expoPtr->name);
                  for (std::string depElem : dependentElements) {
                      exposureMapNames.insert(depElem);
                  }
              }
              else {
                exposureMapNames.insert(expoPtr->name);
              }
          }
      }

      auto replaceThese = degen.replaceWithIdentity(exposureMapNames);

      // Supersede their maps if there are any
      if (!replaceThese.empty()) {
          astrometry::PixelMapCollection pmcAltered;
          for (auto mapname : replaceThese) {
              std::cerr << "..Setting map <" << mapname << "> to Identity" << std::endl;
              pmcAltered.learnMap(astrometry::IdentityMap(mapname));
          }
          // Add the altered maps specs to the input YAML specifications
          std::ostringstream oss;
          pmcAltered.write(oss);
          std::istringstream iss(oss.str());
          inputYAML.addInput(iss, "", true);  // Prepend these specs to others
      }
  }  // End of poly-degeneracy check/correction

  PROGRESS(1, Making final mapCollection);

  // Do not need the preliminary PMC any more.
  delete pmcInit;
  pmcInit = 0;
  // And clean out any old maps stored in the YAMLCollector
  inputYAML.clearMaps();

  // Make final map collection from the YAML
  createMapCollection<Astro>(instruments, exposures, extensions, inputYAML, mapCollection);

  // Add WCS for every extension, and reproject into field coordinates
  PROGRESS(2, Defining all WCSs);
  setupWCS(fields.projections(), instruments, exposures, extensions, mapCollection);

  /////////////////////////////////////////////////////
  //  Initialize map components that were created with default
  //  parameters by fitting to fake data created from the
  //  starting WCS.
  /////////////////////////////////////////////////////

  PROGRESS(2, Initializing defaulted maps);

  // This routine figures out an order in which defaulted maps can
  // be initialized without degeneracies
  std::list<std::set<int>> initializeOrder;
  std::set<int> initializedExtensions;
  {
      MapDegeneracies<Astro> degen(extensions, mapCollection, degenerateTypes,
                                   true);  // Only defaulted maps are used
      // Get the recommended initialization order:
      initializeOrder = degen.initializationOrder();
  }

  for (auto extnSet : initializeOrder) {
      // Fit set of extensions to initialize defaulted map(s)
      std::set<Extension *> defaultedExtensions;
      for (auto iextn : extnSet) {
          defaultedExtensions.insert(extensions[iextn].get());
          initializedExtensions.insert(iextn);
      }
      fitDefaulted(mapCollection, defaultedExtensions, instruments, exposures, verbose >= 2);
  }

  // Try to fit on every extension not already initialized just to
  // make sure that we didn't miss any non-Poly map elements.
  // The fitDefaulted routine will just return if there are no
  // defaulted parameters for the extension.
  for (int iextn = 0; iextn < extensions.size(); iextn++) {
      // Skip extensions that don't exist or are already initialized
      if (!extensions[iextn] || initializedExtensions.count(iextn)) continue;
      std::set<Extension *> defaultedExtensions = {extensions[iextn].get()};
      fitDefaulted(mapCollection, defaultedExtensions, instruments, exposures, verbose >= 2);
      initializedExtensions.insert(iextn);
  }

  // As a check, there should be no more defaulted maps
  bool defaultProblem = false;
  for (auto iName : mapCollection.allMapNames()) {
      if (mapCollection.getDefaulted(iName)) {
          std::cerr << "Logic error: after all intializations, still have map " << iName << " as defaulted."
                    << std::endl;
          defaultProblem = true;
      }
  }
  if (defaultProblem) {
      throw std::runtime_error("Logic error: still have some maps as defaulted");
  }

  // Now fix all map elements requested to be fixed
  fixMapComponents<Astro>(mapCollection, fixMapList, instruments);

  // Recalculate all parameter indices - maps are ready to roll!
  mapCollection.rebuildParameterVector();

  cout << "# Total number of free map elements " << mapCollection.nFreeMaps() << " with "
       << mapCollection.nParams() << " free parameters." << std::endl;

  // Figure out which extensions' maps require a color entry to function
  whoNeedsColor<Astro>(extensions);

  PROGRESS(2, Reprojecting startWcs);

  // Before reading objects, we want to set all starting WCS's to go into
  // field coordinates.
  for (auto const &extnptr : extensions) {
      if (!extnptr) continue;
      if (extnptr->exposure < 0) continue;
      int ifield = exposures[extnptr->exposure]->field;
      extnptr->startWcs->reprojectTo(*fields.projections()[ifield]);
  }
}

void WCSFit::setMatches(const std::vector<int> &sequence, const std::vector<LONGLONG> &extns,
                        const std::vector<LONGLONG> &objects, ExtensionObjectSet skipSet, bool usePM) {
    readMatches<Astro>(sequence, extns, objects, matches, extensions, colorExtensions, skipSet, minMatches,
                       usePM);
}

void WCSFit::setObjects(int i, const map<std::string, std::vector<double>> &tableMap, const std::string &xKey,
                        const std::string &yKey, const std::vector<std::string> &xyErrKeys,
                        const std::string &idKey, const std::string &pmCovKey, const std::string &magKey,
                        const int &magKeyElement, const std::string &magErrKey, const int &magErrKeyElement,
                        const std::string &pmRaKey, const std::string &pmDecKey,
                        const std::string &parallaxKey, const std::vector<std::vector<double>> &pmCov) {
    img::FTable ff;
    for (auto x : tableMap) {
        ff.addColumn<double>(x.second, x.first);
    }
    if (pmCovKey != "") {
        ff.addColumn<vector<double>>(pmCov, pmCovKey);
    }

    readObjects_oneExtension<Astro>(exposures, i, ff, xKey, yKey, idKey, pmCovKey, xyErrKeys, magKey,
                                    magKeyElement, magErrKey, magErrKeyElement, pmRaKey, pmDecKey,
                                    parallaxKey, extensions, fields.projections(), verbose, true);
}

void WCSFit::fit(double maxError, int minFitExposures, double reserveFraction, int randomNumberSeed,
                 double minimumImprovement, double clipThresh, double chisqTolerance, bool clipEntireMatch,
                 bool divideInPlace, bool purgeOutput, double minColor, double maxColor) {
    try {
        PROGRESS(2, Purging defective detections and matches);

        // Get rid of Detections with errors too high
        maxError *= astrometry::RESIDUAL_UNIT / astrometry::WCS_UNIT;
        purgeNoisyDetections<Astro>(maxError, matches, exposures, extensions);
        std::cerr << "Matches size after purgeNoisyDetections: " << matches.size()
                  << std::endl;

        PROGRESS(2, Purging sparse matches);
        // Get rid of Matches with too few detections
        purgeSparseMatches<Astro>(minMatches, matches);
        std::cerr << "Matches size after purgeSparseMatches: " << matches.size() << std::endl;

        PROGRESS(2, Purging out - of - range colors);
        // Get rid of Matches with color out of range (note that default color is 0).
        purgeBadColor<Astro>(minColor, maxColor, matches);
        std::cerr << "Matches size after purgeBadColor: " << matches.size() << std::endl;

        PROGRESS(2, Reserving matches);
        // Reserve desired fraction of matches
        if (reserveFraction > 0.) {
            reserveMatches<Astro>(matches, reserveFraction, randomNumberSeed);
        }

        PROGRESS(2, Purging underpopulated exposures);
        // Find exposures whose parameters are free but have too few
        // Detections being fit to the exposure model.
        auto badExposures = findUnderpopulatedExposures<Astro>(minFitExposures, matches, exposures,
                                                               extensions, mapCollection);

        PROGRESS(2, Purging bad exposure parameters and Detections);
        // Freeze parameters of an exposure model and clip all
        // Detections that were going to use it.
        for (auto i : badExposures) {
            cout << "WARNING: Shutting down exposure map " << i.first << " with only "
                 << i.second << " fitted detections " << std::endl;
            freezeMap<Astro>(i.first, matches, extensions, mapCollection);
        }

        if (purgeOutput) {
            PROGRESS(2, Purging unfittable maps);
            mapCollection.purgeInvalid();
        }

        PROGRESS(2, Match census);
        matchCensus<Astro>(matches, cout);

        ///////////////////////////////////////////////////////////
        // Now do the re-fitting
        ///////////////////////////////////////////////////////////

        PROGRESS(1, Begin fitting process);

        // make CoordAlign class
        astrometry::CoordAlign ca(mapCollection, matches);

        int nclip;
        double oldthresh = 0.;

        // Start off in a "coarse" mode so we are not fine-tuning the solution
        // until most of the outliers have been rejected:
        bool coarsePasses = true;

        ca.setRelTolerance(10. * chisqTolerance);
        // Here is the actual fitting loop
        do {
            // Report number of active Matches / Detections in each iteration:
            {
                std::cerr << "Matches size: " << matches.size() << std::endl;
                long int mcount = 0;
                long int dcount = 0;
                ca.count(mcount, dcount, false, 2);
                double maxdev = 0.;
                int dof = 0;
                double chi = ca.chisqDOF(dof, maxdev, false);
                if (verbose >= 1)
                    std::cerr << "Fitting " << mcount << " matches with "
                              << dcount << " detections "
                              << " chisq " << chi << " / " << dof
                              << " dof,  maxdev " << maxdev << " sigma" << std::endl;
            }

            // Do the fit here!!
            double chisq = ca.fitOnce(verbose >= 1, divideInPlace);  // save space if selected
            // Note that fitOnce() remaps *all* the matches, including reserved ones.
            double max;
            int dof;
            ca.chisqDOF(dof, max, false);                    // Exclude reserved Matches
            double thresh = sqrt(chisq / dof) * clipThresh;  // ??? change dof to expectedChisq?
            if (verbose >= 1)
                std::cerr << "After iteration: chisq " << chisq << " / "
                          << dof << " dof, max deviation " << max
                          << "  new clip threshold at: " << thresh << " sigma" << std::endl;
            if (thresh >= max || (oldthresh > 0. && (1 - thresh / oldthresh) < minimumImprovement)) {
                // Sigma clipping is no longer doing much.  Quit if we are at full precision,
                // else require full target precision and initiate another pass.
                if (coarsePasses) {
                    coarsePasses = false;
                    ca.setRelTolerance(chisqTolerance);
                    PROGRESS(1, Starting strict tolerance passes);
                    if (clipEntireMatch && verbose >= 1) std::cerr << "-->clipping full matches" << std::endl;
                    oldthresh = thresh;
                    nclip = ca.sigmaClip(thresh, false, clipEntireMatch && !coarsePasses, verbose >= 1);
                    if (verbose >= 0)
                        std::cerr << "# Clipped " << nclip << " matches " << std::endl;
                    continue;
                } else {
                    // Done!
                    break;
                }
            }
            oldthresh = thresh;
            // Clip entire matches on final passes if clipEntireMatch=true
            nclip = ca.sigmaClip(thresh, false, clipEntireMatch && !coarsePasses, verbose >= 1);
            if (nclip == 0 && coarsePasses) {
                // Nothing being clipped; tighten tolerances and re-fit
                coarsePasses = false;
                ca.setRelTolerance(chisqTolerance);
                PROGRESS(1, Starting strict tolerance passes);
                if (clipEntireMatch && verbose >= 1) std::cerr << "-->clipping full matches" << std::endl;
                continue;
            }
            if (verbose >= 0) std::cerr << "# Clipped " << nclip << " matches " << std::endl;

        } while (coarsePasses || nclip > 0);

        // If there are reserved Matches, run sigma-clipping on them now.
        if (reserveFraction > 0.) {
            PROGRESS(1, Clipping reserved matches);
            clipReserved<Astro>(ca, clipThresh, minimumImprovement, false, verbose >= 1);
        }

        //////////////////////////////////////
        // Output data and calculate some statistics
        //////////////////////////////////////

        // Report summary of residuals to stdout
        Astro::reportStatistics(matches, exposures, extensions, cout);
    } catch (std::runtime_error &m) {
        quit(m, 1);
    }
}

std::map<std::string, vector<float>> WCSFit::getOutputCatalog() {
    std::map<std::string, vector<float>> outputDict;
    // Save the pointwise fitting results
    outputDict = Astro::getOutputCatalog(matches);
    return outputDict;
}

std::map<std::string, vector<float>> WCSFit::getPMCatalog(vector<vector<float>> &PMMean,
                                                          vector<vector<float>> &PMInvCov) {
    std::map<std::string, vector<float>> outputPM;
    // Save the pointwise fitting results
    outputPM = Astro::getPMCatalog(matches, PMMean, PMInvCov); 
    return outputPM;
}

std::map<std::string, vector<float>> WCSFit::getStarCatalog(vector<vector<float>> &starInvCov) {
    std::map<std::string, vector<float>> starCatalog;
    // Save the pointwise fitting results
    {
        // This routine needs an array of field projections for each extension
        std::vector<astrometry::SphericalCoords *> extensionProjections(extensions.size());
        for (int i = 0; i < extensions.size(); i++) {
            if (!extensions[i]) continue;
            int iExposure = extensions[i]->exposure;
            if (iExposure < 0 || !exposures[iExposure]) continue;
            int iField = exposures[iExposure]->field;
            extensionProjections[i] = fields.projections()[iField].get();
        }
        PROGRESS(2, extensionProjections completed);
        
        starCatalog = Astro::getStarCatalog(matches, extensionProjections, starInvCov); 
    }
    
    return starCatalog;
}

void WCSFit::saveResults(std::string outWcs, std::string outCatalog, std::string starCatalog) {
    // The re-fitting is now complete.  Serialize all the fitted coordinate systems
    PROGRESS(2, Saving astrometric parameters);
    // Save the pointwise fitting results
    {
        ofstream ofs(outWcs.c_str());
        if (!ofs) {
            std::cerr << "Error trying to open output file for fitted Wcs: " << outWcs << std::endl;
            // *** will not quit before dumping output ***
        } else {
            mapCollection.write(ofs);
        }
    }

    PROGRESS(1, Saving astrometric residuals);
    // Save the pointwise fitting results
    {
        // This routine needs an array of field projections for each extension
        std::vector<astrometry::SphericalCoords *> extensionProjections(extensions.size());
        for (int i = 0; i < extensions.size(); i++) {
            if (!extensions[i]) continue;
            int iExposure = extensions[i]->exposure;
            if (iExposure < 0 || !exposures[iExposure]) continue;
            int iField = exposures[iExposure]->field;
            extensionProjections[i] = fields.projections()[iField].get();
        }
        PROGRESS(2, extensionProjections completed);
        Astro::saveResults(matches, outCatalog, starCatalog, extensionProjections);
    }
}
