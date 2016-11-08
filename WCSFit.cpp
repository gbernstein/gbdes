// Program to fit astrometric solutions to detection catalogs already matched by WCSFoF.
#include <fstream>
#include <sstream>
#include <map>
#include <algorithm>

#include "Std.h"
#include "Astrometry.h"
#include "FitsTable.h"
#include "StringStuff.h"
#include "Pset.h"
#include "PixelMapCollection.h"
#include "YAMLCollector.h"
#include "Match.h"
#include "Instrument.h"

#include "FitSubroutines.h"
#include "WcsSubs.h"
#include "MapDegeneracies.h"

#ifdef _OPENMP
#include <omp.h>
#endif

using namespace std;
using namespace astrometry;
using namespace stringstuff;
using img::FTable;

string usage=
  "WCSFit: Refine coordinate solutions for a matched set of catalogs.\n"
  "usage: WCSFit <match file> [parameter file] [parameter file...]\n"
  "   [-parameter[=]value...]\n"
  "      <match file>:  FITS file with binary tables produced by WCSFoF\n"
  "      Program parameters specified as command-line options or read from\n"
  "          parameter file(s) specified on cmd line";

// Temporary documentation:
// Note that this is assuming that regexes do not include = or , characters.
// whitespace will be stripped from edges of each name.

// parameter inputMaps is a string with
// [<mapName>@]<filename>, ...
// which says that maps matching the regex mapName should be deserialized from the YAML
// file filename.  If no mapName is given, anything matches.  The inputMaps are searched
// in order given.  The input maps may be uninitialized (i.e. no parameters given) in which
// case an initial fit based on the starting WCS will be done.  The inputMaps files
// will specify the functional forms used for the coordinate maps.  They may
// contain strings like INSTRUMENT, EXPOSURE, BAND, DEVICE which will be replaced
// from a dictionary.
// Same caveats: no @ or commas in regexes, whitespace stripped.

// parameter fixMaps is a string with
// <mapName>, ...
// where any given mapName should have its parameters fixed at initial values during
// the fitting.  Regexes allowed (no commas!).

// canonicalExposures ???
// are exposures that will be given an identity exposure map in order to break
// the usual degeneracy between exposure and instrument maps.  There must be
// 0 or 1 of these specified for any instrument that has Instrument Map with free parameters
// but no exposures in which the either the instrument map or exposure map is fixed.
// Default is to find an exposure that has data in all devices and use it.
// Will have an error if there is more than one constraint on any Instrument.

// Note that pixel maps for devices within instrument will get names <instrument>/<device>.
// And Wcs's for individual exposures' extension will get names <exposure>/<device>.

int
main(int argc, char *argv[])
{
  double reserveFraction;
  int randomNumberSeed;
  string skipFile;

  double clipThresh;
  double maxPixError;
  double pixSysError;
  double referenceSysError;

  int minMatches;
  int minFitExposures;
  bool clipEntireMatch;
  double chisqTolerance;

  string inputMaps;
  string fixMaps;
  string useInstruments;
  string skipExposures;

  string outCatalog;
  string outWcs;

  string colorExposures;
  double minColor;
  double maxColor;
  
  Pset parameters;
  {
    const int def=PsetMember::hasDefault;
    const int low=PsetMember::hasLowerBound;
    const int up=PsetMember::hasUpperBound;
    const int lowopen = low | PsetMember::openLowerBound;
    const int upopen = up | PsetMember::openUpperBound;

    parameters.addMemberNoValue("INPUTS");
    parameters.addMember("maxPixError",&maxPixError, def | lowopen,
			 "Cut objects with pixel posn uncertainty above this", 0.1, 0.);
    parameters.addMember("pixSysError",&pixSysError, def | low,
			 "Additional systematic error for detections (pixels)", 0.01, 0.);
    parameters.addMember("referenceSysError",&referenceSysError, def | low,
			 "Reference object additional error (arcsec)", 0.003, 0.);
    parameters.addMember("minMatch",&minMatches, def | low,
			 "Minimum number of detections for usable match", 2, 2);
    parameters.addMember("minFitExposures",&minFitExposures, def | low,
			 "Minimum number of detections to fit exposure map", 200, 2);
    parameters.addMember("useInstruments",&useInstruments, def,
			 "the instruments to include in fit",".*");
    parameters.addMember("skipExposures",&skipExposures, def,
			 "exposures to ignore during fitting","");
    parameters.addMemberNoValue("CLIPPING");
    parameters.addMember("clipThresh",&clipThresh, def | low,
			 "Clipping threshold (sigma)", 5., 2.);
    parameters.addMember("clipEntireMatch",&clipEntireMatch, def,
			 "Discard entire object if one outlier on later passes", false);
    parameters.addMember("skipFile",&skipFile, def,
			 "optional file holding extension/object of detections to ignore","");
    parameters.addMemberNoValue("FITTING");
    parameters.addMember("reserveFraction",&reserveFraction, def | low,
			 "Fraction of matches reserved from re-fit", 0., 0.);
    parameters.addMember("seed",&randomNumberSeed, def,
			 "seed for reserving randomizer, <=0 to seed with time", 0);
    parameters.addMember("chisqTolerance",&chisqTolerance, def | lowopen,
			 "Fractional change in chisq for convergence", 0.001, 0.);
    parameters.addMember("inputMaps",&inputMaps, def,
			 "list of YAML files specifying maps","");
    parameters.addMember("fixMaps",&fixMaps, def,
			 "list of map components or instruments to hold fixed","");

    parameters.addMemberNoValue("COLORS");
    parameters.addMember("colorExposures",&colorExposures, def,
			 "exposures holding valid colors for stars","");
    parameters.addMember("minColor",&minColor, def,
			 "minimum value of color to be used",-10.);
    parameters.addMember("maxColor",&maxColor, def,
			 "maximum value of color to be used",+10.);

    parameters.addMemberNoValue("OUTPUTS");
    parameters.addMember("outCatalog",&outCatalog, def,
			 "Output FITS binary catalog", "wcscat.fits");
    parameters.addMember("outWcs",&outWcs, def,
			 "Output serialized Wcs systems", "wcsfit.wcs");
  }

  // Positional accuracy (in degrees) demanded of numerical solutions for inversion of 
  // pixel maps: 
  const double worldTolerance = 0.001*ARCSEC/DEGREE;
  // Fractional reduction in RMS required to continue sigma-clipping:
  const double minimumImprovement=0.02;

  try {
    
    // Read all the command-line and parameter-file program parameters
    processParameters(parameters, usage, 1, argc, argv);
    string inputTables = argv[1];

    referenceSysError *= ARCSEC/DEGREE;

    /////////////////////////////////////////////////////
    // Parse all the parameters
    /////////////////////////////////////////////////////

    // Teach PixelMapCollection about new kinds of PixelMaps:
    loadPixelMapParser();

    // This is list of regexes of PixelMap names (or instrument names) that should
    // have their parameters held fixed.
    list<string> fixMapList = splitArgument(fixMaps);

    // The list of instruments that we will be matching together in this run:
    // have their parameters held fixed.
    list<string> useInstrumentList = splitArgument(useInstruments);

    // The list of exposures that are considered valid sources of color information:
    list<string> useColorList = splitArgument(colorExposures);
    
    // Objects to ignore on input:
    ExtensionObjectSet skipSet(skipFile);

    // Exposures to ignore:
    list<string> skipExposureList = splitArgument(skipExposures);
    
    // Class that will build a starting YAML config for all extensions
    YAMLCollector inputYAML(inputMaps, PixelMapCollection::magicKey);
    // Make sure inputYAML knows about the Identity transformation:
    {
      istringstream iss("Identity:\n  Type:  Identity\n");
      inputYAML.addInput(iss);
    }
    
    /////////////////////////////////////////////////////
    //  Read in properties of all Fields, Instruments, Devices, Exposures
    /////////////////////////////////////////////////////

    // All the names will be stripped of leading/trailing white space, and internal
    // white space replaced with single underscore - this keeps PixelMap parsing functional.

    /**/cerr << "Reading fields" << endl;

    // All we care about fields are names and orientations:
    NameIndex fieldNames;
    vector<SphericalCoords*> fieldProjections;
    // Read the Fields table from input, copy to a new output FITS file, extract needed info
    readFields(inputTables, outCatalog, fieldNames, fieldProjections);


    /**/cerr << "Reading instruments" << endl;

    // Let's figure out which of our FITS extensions are Instrument or MatchCatalog
    vector<int> instrumentHDUs;
    vector<int> catalogHDUs;
    inventoryFitsTables(inputTables, instrumentHDUs, catalogHDUs);
    
    
    // This flag is set since we have already opened (and overwritten) the
    // output FITS catalog.
    bool outputCatalogAlreadyOpen = true;

    // Read in all the instrument extensions and their device info from input
    // FITS file, save useful ones and write to output FITS file.
    vector<Instrument*> instruments =
      readInstruments(instrumentHDUs, useInstrumentList, inputTables, outCatalog,
		      outputCatalogAlreadyOpen);


    /**/cerr << "Reading exposures" << endl;

    // This vector will hold the color-priority value of each exposure.  -1 means an exposure
    // that does not hold color info.
    vector<int> exposureColorPriorities;
    // Read in the table of exposures
    vector<Exposure*> exposures =
      readExposures(instruments,
		    exposureColorPriorities,
		    useColorList,
		    inputTables,
		    outCatalog,
		    skipExposureList,
		    true, // Use reference exposures for astrometry
		    outputCatalogAlreadyOpen);

		    
    /**/cerr << "Reading extensions" << endl;

    // Read info about all Extensions - we will keep the Table around.
    FTable extensionTable;
    {
      FITS::FitsTable ft(inputTables, FITS::ReadOnly, "Extensions");
      extensionTable = ft.extract();
      FITS::FitsTable out(outCatalog, FITS::ReadWrite+FITS::Create, "Extensions");
      out.copy(extensionTable);
    }
    vector<ColorExtension*> colorExtensions;
    vector<Extension*> extensions =
      readExtensions<Astro>(extensionTable,
			    instruments,
			    exposures,
			    exposureColorPriorities,
			    colorExtensions,
			    inputYAML);

		    
    /**/cerr << "Setting reference wcsNames" << endl;

    // A special loop here to set the wcsname of reference extensions to the
    // name of their field.
    for (auto extnptr : extensions) {
      if (!extnptr) continue;
      const Exposure& expo = *exposures[extnptr->exposure];
      if ( expo.instrument >= 0) continue;
      int ifield = expo.field;
      extnptr->wcsName = fieldNames.nameOf(ifield); // ??? mapName instead???
    }

    
    /////////////////////////////////////////////////////
    //  Create and initialize all maps
    /////////////////////////////////////////////////////

    /**/cerr << "Building initial PixelMapCollection" << endl;

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
    fixMapComponents<Astro>(*pmcInit,
			    fixMapList,
			    instruments);
    
    /////////////////////////////////////////////////////
    // First sanity check: Check for maps that are defaulted but fixed
    /////////////////////////////////////////////////////
    /**/cerr << "Checking for fixed+default" << endl;
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

    /**/cerr << "Checking for field degeneracy" << endl;

    {
      vector<bool> fieldHasFree(fieldNames.size(), false);
      vector<bool> fieldHasFixed(fieldNames.size(), false);
      for (auto extnptr : extensions) {
	if (!extnptr) continue; // Not in use
	int field = exposures[extnptr->exposure]->field;
	if (pmcInit->getFixed(extnptr->mapName))
	  fieldHasFixed[field] = true;
	else
	  fieldHasFree[field] = true;
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
    
    /////////////////////////////////////////////////////
    // Next degeneracy: if linear/poly maps are compounded there
    // must be data for which one of them is fixed.  See if it
    // is needed to set some exposure maps to Identity to break
    // such degeneracies.
    /////////////////////////////////////////////////////

    /**/cerr << "Checking for polynomial degeneracies" << endl;
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

    /**/cerr << "Making final mapCollection" << endl;

    // Do not need the preliminary PMC any more.
    delete pmcInit;
    pmcInit = 0;
    // And clean out any old maps stored in the YAMLCollector
    inputYAML.clearMaps();
    
    // Make final map collection from the YAML
    PixelMapCollection mapCollection;
    createMapCollection<Astro>(instruments,
			       exposures,
			       extensions,
			       inputYAML,
			       mapCollection);


    // Add WCS for every extension, and reproject into field coordinates
    /**/cerr << "Defining all WCS's" << endl;
    setupWCS(fieldProjections, instruments, exposures, extensions, mapCollection);

    /////////////////////////////////////////////////////
    //  Initialize map components that were created with default
    //  parameters by fitting to fake data created from the
    //  starting WCS.  
    /////////////////////////////////////////////////////

    /**/cerr << "Initializing defaulted maps" << endl;

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
	defaultedExtensions.insert(extensions[iextn]);
	initializedExtensions.insert(iextn);
      }
      fitDefaulted(mapCollection,
		   defaultedExtensions,
		   instruments,
		   exposures);
    }

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
		   exposures);
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


    //////////////////////////////////////////////////////////
    // Read in all the data
    //////////////////////////////////////////////////////////

    // List of all Matches - they will hold pointers to all Detections too.
    list<Match*> matches;
    
    // Figure out which extensions' maps require a color entry to function
    whoNeedsColor<Astro>(extensions);
    
    /**/cerr << "Reprojecting startWcs" << endl;
    
    // Before reading objects, we want to set all starting WCS's to go into
    // field coordinates.
    for (auto extnptr : extensions) {
      if (!extnptr) continue;
      if (extnptr->exposure < 0) continue;
      int ifield = exposures[extnptr->exposure]->field;
      extnptr->startWcs->reprojectTo(*fieldProjections[ifield]);
    }

    // Start by reading all matched catalogs, creating Detection and Match arrays, and 
    // telling each Extension which objects it should retrieve from its catalog

    for (int icat = 0; icat < catalogHDUs.size(); icat++) {
      FITS::FitsTable ft(inputTables, FITS::ReadOnly, catalogHDUs[icat]);
      FTable ff = ft.use();
      string dummy1, dummy2;
      ff.getHdrValue("Field", dummy1);
      ff.getHdrValue("Affinity", dummy2);
      /**/cerr << "Parsing catalog field " << dummy1 << " Affinity " << dummy2 << endl;
      
      readMatches<Astro>(ff, matches, extensions, colorExtensions, skipSet, minMatches);
      
    } // End loop over input matched catalogs

    /**/cerr << "Total match count: " << matches.size() << endl;

    // Now loop over all original catalog bintables, reading the desired rows
    // and collecting needed information into the Detection structures
    /**/cerr << "Reading catalogs." << endl;
    readObjects<Astro>(extensionTable, exposures, extensions, pixSysError, referenceSysError);


    // Now loop again over all catalogs being used to supply colors,
    // and insert colors into all the Detections they match
    /**/cerr << "Reading colors" << endl;
    readColors<Astro>(extensionTable, colorExtensions);

    /**/cerr << "Purging defective detections and matches" << endl;

    // Get rid of Detections with errors too high
    purgeNoisyDetections<Astro>(maxPixError, matches, exposures, extensions);
			 
    // Get rid of Matches with too few detections
    purgeSparseMatches<Astro>(minMatches, matches);

    // Get rid of Matches with color out of range (note that default color is 0).
    purgeBadColor<Astro>(minColor, maxColor, matches);
    
    // Reserve desired fraction of matches
    if (reserveFraction>0.) 
      reserveMatches<Astro>(matches, reserveFraction, randomNumberSeed);

    // Find exposures whose parameters are free but have too few
    // Detections being fit to the exposure model.
    auto badExposures = findUnderpopulatedExposures<Astro>(minFitExposures,
							   matches,
							   exposures,
							   extensions,
							   mapCollection);

    // Freeze parameters of an exposure model and clip all
    // Detections that were going to use it.
    for (auto i : badExposures) {
      cout << "WARNING: Shutting down exposure map " << i.first
	   << " with only " << i.second
	   << " fitted detections "
	   << endl;
      freezeMap<Astro>(i.first, matches, extensions, mapCollection);
    } 

    matchCensus<Astro>(matches, cout);

    ///////////////////////////////////////////////////////////
    // Now do the re-fitting 
    ///////////////////////////////////////////////////////////

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
	long int mcount=0;
	long int dcount=0;
	ca.count(mcount, dcount, false, 2);
	double maxdev=0.;
	int dof=0;
	double chi= ca.chisqDOF(dof, maxdev, false);
	/**/cerr << "Fitting " << mcount << " matches with " << dcount << " detections "
		 << " chisq " << chi << " / " << dof << " dof,  maxdev " << maxdev 
		 << " sigma" << endl;
      }

      // Do the fit here!!
      double chisq = ca.fitOnce();
      // Note that fitOnce() remaps *all* the matches, including reserved ones.
      double max;
      int dof;
      ca.chisqDOF(dof, max, false);	// Exclude reserved Matches
      double thresh = sqrt(chisq/dof) * clipThresh;
      /**/cerr << "After iteration: chisq " << chisq 
	       << " / " << dof << " dof, max deviation " << max
	       << "  new clip threshold at: " << thresh << " sigma"
	       << endl;
      if (thresh >= max || (oldthresh>0. && (1-thresh/oldthresh)<minimumImprovement)) {
	// Sigma clipping is no longer doing much.  Quit if we are at full precision,
	// else require full target precision and initiate another pass.
	if (coarsePasses) {
	  coarsePasses = false;
	  ca.setRelTolerance(chisqTolerance);
	  /**/cerr << "--Starting strict tolerance passes";
	  /**/if (clipEntireMatch) cerr << "; clipping full matches";
	  /**/cerr << endl;
	  oldthresh = thresh;
	  nclip = ca.sigmaClip(thresh, false, clipEntireMatch && !coarsePasses);
	  /**/cerr << "Clipped " << nclip
		   << " matches " << endl;
	  continue;
	} else {
	  // Done!
	  break;
	}
      }
      oldthresh = thresh;
      // Clip entire matches on final passes if clipEntireMatch=true
      nclip = ca.sigmaClip(thresh, false, clipEntireMatch && !coarsePasses);
      if (nclip==0 && coarsePasses) {
	// Nothing being clipped; tighten tolerances and re-fit
	coarsePasses = false;
	ca.setRelTolerance(chisqTolerance);
	/**/cerr << "--Starting strict tolerance passes";
	/**/if (clipEntireMatch) cerr << "; clipping full matches";
	/**/cerr << endl;
	continue;
      }
      /**/cerr << "Clipped " << nclip
	       << " matches " << endl;
      
    } while (coarsePasses || nclip>0);
  
    // The re-fitting is now complete.  Serialize all the fitted coordinate systems
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

    // If there are reserved Matches, run sigma-clipping on them now.
    if (reserveFraction > 0.) {
      /**/cerr << "** Clipping reserved matches: " << endl;
      clipReserved<Astro>(ca, clipThresh, minimumImprovement,
			  false, true);  //turn on cerr logging
      //**clipEntireMatch, true);  //turn on cerr logging
    }

    //////////////////////////////////////
    // Output data and calculate some statistics
    //////////////////////////////////////

    // Save the pointwise fitting results
    saveResults<Astro>(matches, outCatalog);
    
    /**/cerr << "Saved results to FITS table" << endl;
    
    // Report summary of residuals to stdout
    reportStatistics<Astro>(matches, exposures, extensions, cout);

    //////////////////////////////////////
    // Cleanup:
    //////////////////////////////////////

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
    for (int i=0; i<exposures.size(); i++)
      if (exposures[i]) delete exposures[i];
    // Get rid of instruments
    for (int i=0; i<instruments.size(); i++)
      if (instruments[i]) delete instruments[i];

  } catch (std::runtime_error& m) {
    quit(m,1);
  }
}


