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

// parameter canonicalExposures is a string with
// <exposureID>, ....
// which are are exposures that will be given an identity exposure map in order to break
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
			 "Discard entire object if one outlier", false);
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

    // All we care about fields are names and orientations:
    NameIndex fieldNames;
    vector<SphericalCoords*> fieldProjections;
    // Read the Fields table from input, copy to a new output FITS file, extract needed info
    readFields(inputTables, outCatalog, fieldNames, fieldProjections);

    /**/cerr << "Done readFields" << endl;

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

    /**/cerr << "Done readInstruments" << endl;

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

    /**/cerr << "Done readExposures" << endl;
		    
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

    /**/cerr << "Done readExtensions" << endl;
		    
    // A special loop here to set the wcsname of reference extensions to the
    // name of their field.
    for (auto extnptr : extensions) {
      if (!extnptr) continue;
      const Exposure& expo = *exposures[extnptr->exposure];
      if ( expo.instrument >= 0) continue;
      int ifield = expo.field;
      extnptr->wcsName = fieldNames.nameOf(ifield); // ??? mapName instead???
    }

    /**/cerr << "Done setting reference wcsNames" << endl;
    
    /////////////////////////////////////////////////////
    //  Create and initialize all magnitude maps
    /////////////////////////////////////////////////////

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

    /**/cerr << "Successfully built initial PixelMapCollection" << endl;
    
    // Check every map name against the list of those to fix.
    // Includes all devices of any instruments on the fixMapList.
    fixMapComponents<Astro>(*pmcInit,
			    fixMapList,
			    instruments);
    
    /////////////////////////////////////////////////////
    // First sanity check: Check for maps that are defaulted but fixed
    /////////////////////////////////////////////////////
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
    
    /**/cerr << "Done fixed+default check" << endl;

    /////////////////////////////////////////////////////
    // Second sanity check: in each field that has *any* free maps,
    // do we have at least one map that is *fixed*?
    /////////////////////////////////////////////////////
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
    
    /**/cerr << "Done field degeneracy check" << endl;

    /////////////////////////////////////////////////////
    // Now for each instrument, determine if we need a canonical
    // exposure to break an exposure/instrument degeneracy.  If
    // so, choose one, and replace its exposure map with Identity.
    /////////////////////////////////////////////////////

    // Here's where we'll collect the map definitions that override the
    // inputs in order to eliminate degeneracies.
    PixelMapCollection pmcAltered;  

    for (int iInst=0; iInst < instruments.size(); iInst++) {
      if (!instruments[iInst]) continue;  // Not using instrument
      auto& instr = *instruments[iInst];

      int canonicalExposure =
	findCanonical<Astro>(instr, iInst, exposures, extensions, *pmcInit);

      if (canonicalExposure >= 0)
	// Make a new map spec for the canonical exposure
	pmcAltered.learnMap(IdentityMap(exposures[canonicalExposure]->name));
    } // End instrument loop

    // Re-read all of the maps, assigning WCS to each extension this time
    // and placing into the final PixelMapCollection.

    // Add the altered maps specs to the input YAML specifications
    {
      ostringstream oss;
      pmcAltered.write(oss);
      istringstream iss(oss.str());
      inputYAML.addInput(iss, "", true); // Prepend these specs to others
    }
    
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

    /**/cerr << "Done making final mapCollection!" << endl;

    // Add WCS for every extension, and reproject into field coordinates
    setupWCS(fieldProjections, instruments, exposures, extensions, mapCollection);

    /**/cerr << "Done defining all WCS's" << endl;

    /////////////////////////////////////////////////////
    //  Initialize map components that were created with default
    //  parameters by fitting to fake data created from the
    //  starting WCS.  
    /////////////////////////////////////////////////////

    // Will do this first by finding any device maps that are defaulted,
    // and choosing an exposure with which to initialize them.
    // The chosen exposure should have non-defaulted exposure solution.

    list<int> exposuresToInitialize =
      pickExposuresToInitialize(instruments,
				exposures,
				extensions,
				mapCollection);
    
    // Initialize defaulted maps, one exposure at a time, starting
    // with the ones that solve for devices, then all the others.
    // The defaulted status
    // in the mapCollection should be getting updated as we go.
    for (int i=0; i<exposures.size(); i++)
      if (exposures[i])
	exposuresToInitialize.push_back(i);

    for (auto iExpo : exposuresToInitialize) {
      //   Find all defaulted extensions using this exposure
      set<Extension*> defaultedExtensions;
      for (auto extnptr : extensions) {
	if (!extnptr) continue; // Not in use
	if (extnptr->exposure != iExpo)
	  continue;
	if (mapCollection.getDefaulted(extnptr->mapName))
	  defaultedExtensions.insert(extnptr);
      }
      if (!defaultedExtensions.empty()) {
	fitDefaulted(mapCollection,
		     defaultedExtensions,
		     instruments,
		     exposures);
      }
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

    cout << "Total number of free map elements " << mapCollection.nFreeMaps()
	 << " with " << mapCollection.nParams() << " free parameters."
	 << endl;


    //////////////////////////////////////////////////////////
    // Read in all the data
    //////////////////////////////////////////////////////////

    // List of all Matches - they will hold pointers to all Detections too.
    list<Match*> matches;
    
    // Figure out which extensions' maps require a color entry to function
    whoNeedsColor<Astro>(extensions);
    
    // Before reading objects, we want to set all starting WCS's to go into
    // field coordinates.
    for (auto extnptr : extensions) {
      if (!extnptr) continue;
      if (extnptr->exposure < 0) continue;
      int ifield = exposures[extnptr->exposure]->field;
      extnptr->startWcs->reprojectTo(*fieldProjections[ifield]);
    }

    /**/cerr << "Done reprojecting startWcs" << endl;
    
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
    readObjects<Astro>(extensionTable, exposures, extensions, pixSysError, referenceSysError);

    /**/cerr << "Done reading catalogs." << endl;

    // Now loop again over all catalogs being used to supply colors,
    // and insert colors into all the PhotoArguments for Detections they match
    //** ???    readColors<Astro>(extensionTable, colorExtensions);

    /**/cerr << "Done reading colors" << endl;

    // Get rid of Detections with errors too high
    purgeNoisyDetections<Astro>(maxPixError, matches, exposures, extensions);
			 
    // Get rid of Matches with too few detections
    purgeSparseMatches<Astro>(minMatches, matches);

    // Get rid of Matches with color out of range (note that default color is 0).
    purgeBadColor<Astro>(minColor, maxColor, matches); // ??? Nop right now
    
    /**/cerr << "Done purging defective detections and matches" << endl;

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

    matchCensus<Astro>(matches);

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
	  /**/cerr << "--Starting strict tolerance passes, clipping full matches" << endl;
	  oldthresh = thresh;
	  nclip = ca.sigmaClip(thresh, false, true);
	  /**/cerr << "Clipped " << nclip
		   << " matches " << endl;
	  continue;
	} else {
	  // Done!
	  break;
	}
      }
      oldthresh = thresh;
      nclip = ca.sigmaClip(thresh, false, clipEntireMatch || !coarsePasses);
      if (nclip==0 && coarsePasses) {
	// Nothing being clipped; tighten tolerances and re-fit
	coarsePasses = false;
	ca.setRelTolerance(chisqTolerance);
	/**/cerr << "No clipping--Starting strict tolerance passes, clipping full matches" << endl;
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
      clipReserved(ca, clipThresh, minimumImprovement,
		   clipEntireMatch, true);  //turn on cerr logging
    }

    //////////////////////////////////////
    // Output data and calculate some statistics
    //////////////////////////////////////

    // Create Accum instances for fitted and reserved Detections on every
    // exposure, plus total accumulator for all reference
    // and all non-reference objects.
    vector<Accum> vaccFit(exposures.size());
    vector<Accum> vaccReserve(exposures.size());
    Accum refAccFit;
    Accum refAccReserve;
    Accum accFit;
    Accum accReserve;

    // Open the output bintable
    FITS::FitsTable ft(outCatalog, FITS::ReadWrite + FITS::Create, "WCSOut");
    FTable outTable = ft.use();;

    // Create vectors that will be put into output table
    vector<int> sequence;
    outTable.addColumn(sequence, "SequenceNumber");
    vector<long> catalogNumber;
    outTable.addColumn(catalogNumber, "Extension");
    vector<long> objectNumber;
    outTable.addColumn(objectNumber, "Object");
    vector<bool> clip;
    outTable.addColumn(clip, "Clip");
    vector<bool> reserve;
    outTable.addColumn(reserve, "Reserve");
    vector<float> xpix;
    outTable.addColumn(xpix, "xPix");
    vector<float> ypix;
    outTable.addColumn(ypix, "yPix");
    vector<float> sigpix;
    outTable.addColumn(sigpix, "sigPix");
    vector<float> xrespix;
    outTable.addColumn(xrespix, "xresPix");
    vector<float> yrespix;
    outTable.addColumn(yrespix, "yresPix");

    vector<float> xw;
    outTable.addColumn(xw, "xW");
    vector<float> yw;
    outTable.addColumn(yw, "yW");
    vector<float> sigw;
    outTable.addColumn(sigw, "sigW");
    vector<float> xresw;
    outTable.addColumn(xresw, "xresW");
    vector<float> yresw;
    outTable.addColumn(yresw, "yresW");

    vector<float> wtFrac;
    outTable.addColumn(wtFrac, "wtFrac");

    // Cumulative counter for rows written to table:
    long pointCount = 0;
    // Write vectors to table when this many rows accumulate:
    const long WriteChunk = 100000;

    // Write all matches to output catalog, deleting them along the way
    // and accumulating statistics of each exposure.
    // 
    list<Match*>::iterator im = matches.begin();
    while ( im != matches.end()) {
      // First, write current vectors to table if they've gotten big
      if ( sequence.size() > WriteChunk) {
	outTable.writeCells(sequence, "SequenceNumber", pointCount);
	outTable.writeCells(catalogNumber, "Extension", pointCount);
	outTable.writeCells(objectNumber, "Object", pointCount);
	outTable.writeCells(clip, "Clip", pointCount);
	outTable.writeCells(reserve, "Reserve", pointCount);
	outTable.writeCells(xpix, "xPix", pointCount);
	outTable.writeCells(ypix, "yPix", pointCount);
	outTable.writeCells(xrespix, "xresPix", pointCount);
	outTable.writeCells(yrespix, "yresPix", pointCount);
	outTable.writeCells(xw, "xW", pointCount);
	outTable.writeCells(yw, "yW", pointCount);
	outTable.writeCells(xresw, "xresW", pointCount);
	outTable.writeCells(yresw, "yresW", pointCount);
	outTable.writeCells(sigpix, "sigPix", pointCount);
	outTable.writeCells(sigw, "sigW", pointCount);
	outTable.writeCells(wtFrac, "wtFrac", pointCount);

	pointCount += sequence.size();

	sequence.clear();
	catalogNumber.clear();
	objectNumber.clear();
	clip.clear();
	reserve.clear();
	xpix.clear();
	ypix.clear();
	xrespix.clear();
	yrespix.clear();
	xw.clear();
	yw.clear();
	xresw.clear();
	yresw.clear();
	sigpix.clear();
	sigw.clear();
	wtFrac.clear();
      }	// Done flushing the vectors to Table

      Match* m = *im;
      double xc, yc;
      double wtotx, wtoty;
      m->centroid(xc, yc, wtotx, wtoty);
      // Calculate number of DOF per Detection coordinate after
      // allowing for fit to centroid:
      int nFit = m->fitSize();
      double dofPerPt = (nFit > 1) ? 1. - 1./nFit : 0.;
	    
      int detcount=0;
      for (Match::const_iterator idet=m->begin();
	   idet != m->end();
	   ++idet, ++detcount) {
	const Detection* d = *idet;

	// Prepare output quantities
	sequence.push_back(detcount);
	catalogNumber.push_back(d->catalogNumber);
	objectNumber.push_back(d->objectNumber);
	clip.push_back(d->isClipped);
	reserve.push_back(m->getReserved());
	wtFrac.push_back( d->isClipped ? 0. : (d->wtx + d->wty) / (wtotx+wtoty));  //Just take mean of x & y

	xpix.push_back(d->xpix);
	ypix.push_back(d->ypix);
	sigpix.push_back(d->sigma);

	xw.push_back(d->xw);
	yw.push_back(d->yw);
	double sigmaWorld=0.5*(d->clipsqx + d->clipsqy);
	sigmaWorld = (sigmaWorld > 0.) ? 1./sqrt(sigmaWorld) : 0.;
	sigw.push_back(sigmaWorld*1000.*DEGREE/ARCSEC);

	// Calculate residuals if we have a centroid for the match:
	double xcpix=0., ycpix=0.;
	double xerrw=0., yerrw=0.;
	double xerrpix=0., yerrpix=0.;

	if (xc!=0. || yc!=0.) {
	  xcpix = d->xpix;
	  ycpix = d->ypix;
	  Assert (d->map);
	  d->map->toPix(xc, yc, xcpix, ycpix);

	  xerrw = d->xw - xc;
	  yerrw = d->yw - yc;
	  xerrpix = d->xpix - xcpix;
	  yerrpix = d->ypix - ycpix;

	  // Accumulate statistics for meaningful residuals
	  if (dofPerPt >= 0. && !d->isClipped) {
	    int exposureNumber = extensions[d->catalogNumber]->exposure;
	    Exposure* expo = exposures[exposureNumber];
	    if (m->getReserved()) {
	      if (expo->instrument==REF_INSTRUMENT) {
		refAccReserve.add(d, xc, yc, dofPerPt);
		vaccReserve[exposureNumber].add(d, xc, yc, dofPerPt);
	      } else if (expo->instrument==TAG_INSTRUMENT) {
		// do nothing
	      } else {
		accReserve.add(d, xc, yc, dofPerPt);
		vaccReserve[exposureNumber].add(d, xc, yc, dofPerPt);
	      }
	    } else {
	      // Not a reserved object:
	      if (expo->instrument==REF_INSTRUMENT) {
		refAccFit.add(d, xc, yc, dofPerPt);
		vaccFit[exposureNumber].add(d, xc, yc, dofPerPt);
	      } else if (expo->instrument==TAG_INSTRUMENT) {
		// do nothing
	      } else {
		accFit.add(d, xc, yc, dofPerPt);
		vaccFit[exposureNumber].add(d, xc, yc, dofPerPt);
	      }
	    }
	  } // end statistics accumulation

	} // End residuals calculation

	// Put world residuals into milliarcsec
	xrespix.push_back(xerrpix);
	yrespix.push_back(yerrpix);
	xresw.push_back(xerrw*1000.*DEGREE/ARCSEC);
	yresw.push_back(yerrw*1000.*DEGREE/ARCSEC);

      } // End detection loop

      // Done with this match, delete it with its Detections
      m->clear(true);
      // And get rid of match itself.
      im = matches.erase(im);
    } // end match loop

    // Write remaining results to output table:
    outTable.writeCells(sequence, "SequenceNumber", pointCount);
    outTable.writeCells(catalogNumber, "Extension", pointCount);
    outTable.writeCells(objectNumber, "Object", pointCount);
    outTable.writeCells(clip, "Clip", pointCount);
    outTable.writeCells(reserve, "Reserve", pointCount);
    outTable.writeCells(xpix, "xPix", pointCount);
    outTable.writeCells(ypix, "yPix", pointCount);
    outTable.writeCells(xrespix, "xresPix", pointCount);
    outTable.writeCells(yrespix, "yresPix", pointCount);
    outTable.writeCells(xw, "xW", pointCount);
    outTable.writeCells(yw, "yW", pointCount);
    outTable.writeCells(xresw, "xresW", pointCount);
    outTable.writeCells(yresw, "yresW", pointCount);
    outTable.writeCells(sigpix, "sigPix", pointCount);
    outTable.writeCells(sigw, "sigW", pointCount);
    outTable.writeCells(wtFrac, "wtFrac", pointCount);

    // Print summary statistics for each exposure
    cout << "#    Exp    N    DOF    dx    +-    dy    +-   RMS chi_red  "
      "xpix  ypix      xw       yw \n"
	 << "#                     |.....milliarcsec...........|                     |....degrees....|"
	 << endl;
    for (int iexp=0; iexp<exposures.size(); iexp++) {
      // ??? Sort here, omit unused exposures ???
      cout << "Fit     " << setw(3) << iexp
	   << " " << vaccFit[iexp].summary()
	   << endl;
      if (reserveFraction > 0. && vaccReserve[iexp].n > 0) 
	  cout << "Reserve " << setw(3) << iexp
	       << " " << vaccReserve[iexp].summary()
	       << endl;
    } // exposure summary loop
    
    // Output summary data for reference catalog and detections
    cout << "# " << endl;
    cout << "#                   N    DOF    dm    +-    RMS chi_red  \n"
	 << "#                             |.....millimag......|         "
	 << endl;
    cout << "Reference fit     " << refAccFit.summary() << endl;
    if (reserveFraction>0. && refAccReserve.n>0)
      cout << "Reference reserve " << refAccReserve.summary() << endl;

    cout << "Detection fit     " << accFit.summary() << endl;
    if (reserveFraction>0. && accReserve.n>0)
      cout << "Detection reserve " << accReserve.summary() << endl;

    // Cleanup:

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


