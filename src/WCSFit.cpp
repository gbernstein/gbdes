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

#define PROGRESS(val, msg) if (verbose>=val) cerr << "-->" <<  #msg << endl

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
  double maxError;
  double sysError;
  double referenceSysError;
  bool   freePM; 
  double pmEpoch;
  double parallaxPrior;
  double pmPrior;

  int minMatches;
  int minFitExposures;
  bool clipEntireMatch;
  double chisqTolerance;
  bool divideInPlace;
  bool purgeOutput;

  string inputMaps;
  string fixMaps;
  string useInstruments;
  string skipExposures;

  string outCatalog;
  string outWcs;
  string starCatalog;

  string colorExposures;
  double minColor;
  double maxColor;

  int verbose;
  
  Pset parameters;
  {
    const int def=PsetMember::hasDefault;
    const int low=PsetMember::hasLowerBound;
    const int up=PsetMember::hasUpperBound;
    const int lowopen = low | PsetMember::openLowerBound;
    const int upopen = up | PsetMember::openUpperBound;

    parameters.addMemberNoValue("INPUTS");

    // These 3 are assumed to be in RESIDUAL_UNIT from Units.h:
    parameters.addMember("maxError",&maxError, def | lowopen,
			 "Cut objects with posn uncertainty above this (mas)", 100., 0.);
    parameters.addMember("sysError",&sysError, def | low,
			 "Additional systematic error for detections (mas)", 2., 0.);
    parameters.addMember("referenceSysError",&referenceSysError, def | low,
			 "Additional systematic error for non-PM reference objects (mas)", 2., 0.);

    parameters.addMember("freePM",&freePM, def,
			 "Allow free proper motion and parallax?", true);
    parameters.addMember("pmEpoch",&pmEpoch, def,
			 "Time origin for proper motion (2015.5)", 2015.5);
    parameters.addMember("parallaxPrior",&parallaxPrior, def | low,
			 "Prior on parallax for each star (mas)", 10., 0.);
    parameters.addMember("pmPrior",&pmPrior, def | low,
			 "Prior on proper motion per axis for each star (mas/yr)", 100., 0.);
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
    parameters.addMember("divideInPlace",&divideInPlace, def,
			 "use in-place Cholesky to save memory but lose debug of degeneracies",false);
    parameters.addMemberNoValue("FITTING");
    parameters.addMember("reserveFraction",&reserveFraction, def | low,
			 "Fraction of matches reserved from fit", 0., 0.);
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
    parameters.addMember("purgeOutput",&purgeOutput, def,
			 "Purge un-fittable maps from output", false);
    parameters.addMember("outWcs",&outWcs, def,
			 "Output serialized Wcs systems", "wcsfit.wcs");
    parameters.addMember("outCatalog",&outCatalog, def,
			 "Output FITS binary catalog", "wcscat.fits");
    parameters.addMember("starCatalog",&starCatalog, def,
			 "Output stellar PM catalog", "starcat.fits");
    parameters.addMember("verbose", &verbose, def,
			 "stderr detail level", 1);
  }

  // Positional accuracy demanded of numerical solutions for inversion of 
  // pixel maps: 
  const double worldTolerance = 0.1*MILLIARCSEC/WCS_UNIT;
  // Fractional reduction in RMS required to continue sigma-clipping:
  const double minimumImprovement=0.02;

  try {
    
    // Read all the command-line and parameter-file program parameters
    processParameters(parameters, usage, 1, argc, argv);
    string inputTables = argv[1];

    // Convert error parameters from I/O units to internal
    referenceSysError *= RESIDUAL_UNIT/WCS_UNIT;
    sysError *= RESIDUAL_UNIT/WCS_UNIT;
    maxError *= RESIDUAL_UNIT/WCS_UNIT;

    PMMatch::setPrior(pmPrior, parallaxPrior);


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

    PROGRESS(1,Reading fields);
      
    // All we care about fields are names and orientations:
    NameIndex fieldNames;
    vector<SphericalCoords*> fieldProjections;
    vector<double> fieldEpochs;
    // Read the Fields table from input, copy to a new output FITS file, extract needed info
    readFields(inputTables, outCatalog, fieldNames, fieldProjections,
	       fieldEpochs, pmEpoch);

    PROGRESS(1,Reading instruments);

    // Let's figure out which of our FITS extensions are Instrument or MatchCatalog
    vector<int> instrumentHDUs;
    vector<int> catalogHDUs;
    inventoryFitsTables(inputTables, instrumentHDUs, catalogHDUs);
    
    
    // This flag is set since we have already opened (and overwritten) the
    // output FITS catalog.
    bool outputCatalogAlreadyOpen = true;

    // Read in all the instrument extensions and their device info from input
    // FITS file, save useful ones and write to output FITS file.
    vector<unique_ptr<Instrument>> instruments =
      readInstruments(instrumentHDUs, useInstrumentList, inputTables, outCatalog,
		      outputCatalogAlreadyOpen);


    PROGRESS(1,Reading exposures);

    // This vector will hold the color-priority value of each exposure.
    // -1 means an exposure that does not hold color info.
    vector<int> exposureColorPriorities;
    // Read in the table of exposures
    vector<Exposure*> exposures =
      readExposures(instruments,
		    fieldEpochs,
		    exposureColorPriorities,
		    useColorList,
		    inputTables,
		    outCatalog,
		    skipExposureList,
		    true, // Use reference exposures for astrometry
		    outputCatalogAlreadyOpen);

    if (sysError > 0.) {
      Matrix22 astrometricCovariance(0.);
      astrometricCovariance(0,0) = sysError*sysError;
      astrometricCovariance(1,1) = sysError*sysError;
      for (auto e : exposures) {
	if (e && e->instrument >= 0)
	  e->astrometricCovariance += astrometricCovariance;
      }
    }
    if (referenceSysError > 0.) {
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
		    
    PROGRESS(1,Reading extensions);

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
			    inputYAML,
			    verbose>=1);  // Print reading progress?

		    
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
    PROGRESS(2,Defining all WCSs);

    setupWCS(fieldProjections, instruments, exposures, extensions, mapCollection);

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
	defaultedExtensions.insert(extensions[iextn]);
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


    //////////////////////////////////////////////////////////
    // Read in all the data
    //////////////////////////////////////////////////////////

    // List of all Matches - they will hold pointers to all Detections too.
    MCat matches;
    
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

    // Start by reading all matched catalogs, creating Detection and Match arrays, and 
    // telling each Extension which objects it should retrieve from its catalog

    for (int icat = 0; icat < catalogHDUs.size(); icat++) {
      FITS::FitsTable ft(inputTables, FITS::ReadOnly, catalogHDUs[icat]);
      FTable ff = ft.use();
      string dummy1, affinity;
      ff.getHdrValue("Field", dummy1);
      ff.getHdrValue("Affinity", affinity);
      stringstuff::stripWhite(affinity);
      
      // Only use STELLAR affinity for astrometry
      if (!stringstuff::nocaseEqual(affinity,stellarAffinity))
	continue;
      if (verbose>=2)
	cerr << "-->Parsing catalog field " << dummy1 << " Affinity " << affinity << endl;

      // Set this true if we are going to want to create PMMatches from
      // this extension's matches:
      bool usePM = freePM;
      
      readMatches<Astro>(ff, matches, extensions, colorExtensions, skipSet, minMatches,
			 usePM);
      
    } // End loop over input matched catalogs

    if (verbose>=0) cout << "# Total match count: " << matches.size() << endl;

    // Now loop over all original catalog bintables, reading the desired rows
    // and collecting needed information into the Detection structures
    PROGRESS(1,Reading catalogs);
    readObjects<Astro>(extensionTable, exposures, extensions,fieldProjections);

    // Now loop again over all catalogs being used to supply colors,
    // and insert colors into all the Detections they match
    PROGRESS(1,Reading colors);
    readColors<Astro>(extensionTable, colorExtensions);

    PROGRESS(2,Purging defective detections and matches);

    // Get rid of Detections with errors too high
    purgeNoisyDetections<Astro>(maxError,
				matches, exposures, extensions);
			 
    PROGRESS(2,Purging sparse matches);
    // Get rid of Matches with too few detections
    purgeSparseMatches<Astro>(minMatches, matches);

    PROGRESS(2,Purging out-of-range colors);
    // Get rid of Matches with color out of range (note that default color is 0).
    purgeBadColor<Astro>(minColor, maxColor, matches);
    
    PROGRESS(2,Reserving matches);
    // Reserve desired fraction of matches
    if (reserveFraction>0.) 
      reserveMatches<Astro>(matches, reserveFraction, randomNumberSeed);

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
	   << " with only " << i.second
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
    // Save the pointwise fitting results
    {
      // This routine needs an array of field projections for each extension
      vector<SphericalCoords*> extensionProjections(extensions.size(),
							nullptr);
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
      Astro::saveResults(matches, outCatalog, starCatalog,extensionProjections);
    }
    
    PROGRESS(2,Saving FITS tables);
    // Report summary of residuals to stdout
    Astro::reportStatistics(matches, exposures, extensions, cout);

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
    for (int i=0; i<exposures.size(); i++)
      if (exposures[i]) delete exposures[i];

  } catch (std::runtime_error& m) {
    quit(m,1);
  }
}


