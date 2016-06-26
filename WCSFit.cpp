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
#include "Random.h"
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
  bool clipEntireMatch;
  double chisqTolerance;

  string inputMaps;
  string fixMaps;
  string useInstruments;

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
    parameters.addMember("useInstruments",&useInstruments, def,
			 "the instruments to include in fit","");
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

#ifdef COLORS
    parameters.addMemberNoValue("COLORS");
    parameters.addMember("colorExposures",&colorExposures, def,
			 "exposures holding valid colors for stars","");
    parameters.addMember("minColor",&minColor, def,
			 "minimum value of color to be used",-10.);
    parameters.addMember("maxColor",&maxColor, def,
			 "maximum value of color to be used",+10.);
#endif
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
		    true, // Use reference exposures for astrometry
		    outputCatalogAlreadyOpen);
		    
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

    // A special loop here to set the wcsname of reference extensions to the
    // name of their field.
    for (auto extnptr : extensions) {
      if (!extnptr) continue;
      const Exposure& expo = *exposures[extnptr->exposure];
      if ( expo.instrument >= 0) continue;
      int ifield = expo.field;
      extnptr->wcsName = fieldNames.nameOf(ifield);
    }

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
    for (auto extnptr : extensions) {
      if (!extnptr) continue; // Not in use
      Exposure& expo = *exposures[extnptr->exposure];
      int ifield = expo.field;
      if ( expo.instrument < 0) {
	// Tag & reference exposures have no instruments and no fitting
	// being done.  Coordinates are fixed to xpix = xw.
	// Build a simple Wcs that takes its name from the field
	SphericalICRS icrs;
	if (!mapCollection.wcsExists(extnptr->wcsName)) {
	  // If we are the first reference/tag exposure in this field:
	  mapCollection.defineWcs(extnptr->wcsName, icrs, 
				  extnptr->mapName,
				  DEGREE);
	  auto wcs = mapCollection.issueWcs(extnptr->wcsName);
	  // And have this Wcs reproject into field coordinates, learn as map
	  wcs->reprojectTo(*fieldProjections[ifield]);
	  mapCollection.learnMap(*wcs, false, false);
	}
	extnptr->wcs = mapCollection.issueWcs(extnptr->wcsName);
	extnptr->map = mapCollection.issueMap(extnptr->wcsName);
      } else {
	// Real instrument, WCS goes into its exposure coordinates
	mapCollection.defineWcs(extnptr->wcsName, *expo.projection,
				extnptr->mapName,
				DEGREE);
	extnptr->wcs = mapCollection.issueWcs(extnptr->wcsName);

	// Reproject this Wcs into the field system and get a SubMap including reprojection:
	extnptr->wcs->reprojectTo(*fieldProjections[ifield]);
	mapCollection.learnMap(*extnptr->wcs, false, false);
	extnptr->map = mapCollection.issueMap(extnptr->wcsName);
      }
    } // end extension loop

    /**/cerr << "Done defining all WCS's" << endl;

    /////////////////////////////////////////////////////
    //  Initialize map components that were created with default
    //  parameters by fitting to fake data created from the
    //  starting WCS.  
    /////////////////////////////////////////////////////

    // Will do this first by finding any device maps that are defaulted,
    // and choosing an exposure with which to initialize them.
    // The chosen exposure should have non-defaulted exposure solution.

    // We'll put these into a list of exposures to de-default first.
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
	if (!mapCollection.mapExists(mapName)) {
	  unusedDevices.insert(iDev);
	  continue;
	}
	if (mapCollection.getDefaulted(mapName)) {
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
	  if (!mapCollection.mapExists(mapName)) {
	    unusedExposures.insert(iExpo);
	  } else if (mapCollection.getDefaulted(mapName)) {
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
	  if (mapCollection.dependsOn(extnptr->mapName, exposures[iExpo]->name)
	      && mapCollection.dependsOn(extnptr->mapName, instr.mapNames[iDev])) {
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

    // Now fix all map elements requested to be fixed, checking for
    // any remaining defaulted maps at this time.
    {
      set<string> fixTheseMaps;
      bool defaultProblem=false;
      for (auto iName : mapCollection.allMapNames()) {
	if (regexMatchAny(fixMapList, iName))
	  fixTheseMaps.insert(iName);
	if (mapCollection.getDefaulted(iName)) {
	  cerr << "Logic error: after all intializations, still have map "
	       << iName
	       << " as defaulted."
	       << endl;
	  defaultProblem = true;
	}
      }
      if (defaultProblem) exit(1);
      
      // Add names of all devices of instruments on the fixMapList
      for (auto instptr : instruments) {
	if (!instptr) continue; // Not in use
	if (regexMatchAny(fixMapList, instptr->name)) {
	  // Loop through all devices
	  for (int i=0; i<instptr->nDevices; i++) {
	    string devMap = instptr->name + "/" + instptr->deviceNames.nameOf(i);
	    // Freeze the device's map if it's in use
	    if (mapCollection.mapExists(devMap))
	      fixTheseMaps.insert(devMap);
	  }
	}
      }
      mapCollection.setFixed(fixTheseMaps);
    }
    
    // Recalculate all parameter indices
    mapCollection.rebuildParameterVector();

    /**/cerr << "Total number of free map elements " << mapCollection.nFreeMaps()
	     << " with " << mapCollection.nParams() << " free parameters."
	     << endl;


    //////////////////////////////////////////////////////////
    // Read in all the data
    //////////////////////////////////////////////////////////

    // List of all Matches - they will hold pointers to all Detections too.
    list<Match*> matches;    // ??? Store objects, not pointers??

    // Start by reading all matched catalogs, creating Detection and Match arrays, and 
    // telling each Extension which objects it should retrieve from its catalog

    for (int icat = 0; icat < catalogHDUs.size(); icat++) {
      /**/cerr << "# Reading catalog extension " << catalogHDUs[icat] << endl;
      FITS::FitsTable ft(inputTables, FITS::ReadOnly, catalogHDUs[icat]);
      FTable ff = ft.use();
      vector<int> seq;
      vector<LONGLONG> extn;
      vector<LONGLONG> obj;
      ff.readCells(seq, "SequenceNumber");
      ff.readCells(extn, "Extension");
      ff.readCells(obj, "Object");
      // Smaller collections for each match
      vector<long> extns;
      vector<long> objs;
      for (int i=0; i<=seq.size(); i++) {
	if ( (i==seq.size() || seq[i]==0) && !extns.empty()) {
	  // Make a match from previous few entries
	  Detection* d = new Detection;
	  d->catalogNumber = extns[0];
	  d->objectNumber = objs[0];
	  matches.push_back(new Match(d));
	  extensions[extns[0]]->keepers.insert(std::pair<long, Detection*>(objs[0], d));
	  for (int j=1; j<extns.size(); j++) {
	    d = new Detection;
	    d->catalogNumber = extns[j];
	    d->objectNumber = objs[j];
	    matches.back()->add(d);
	    extensions[extns[j]]->keepers.insert(std::pair<long, Detection*>(objs[j], d));
	  }
	  extns.clear();
	  objs.clear();
	} // Finished processing previous match
	if (i>=seq.size()) break;
	// Should we ignore this object?
	if (skipSet(extn[i], obj[i])) continue;

	// Add to things being used:
	extns.push_back(extn[i]);
	objs.push_back(obj[i]);
      } // End loop of catalog entries
      
    } // End loop over input matched catalogs

    // Now loop over all original catalog bintables, reading the desired rows
    // and collecting needed information into the Detection structures

    // Should be safe to multithread this loop as different threads write
    // only to distinct parts of memory.  Perhaps protect the FITS table read?

#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic,4)
#endif
    for (int iext = 0; iext < extensions.size(); iext++) {
      if (!extensions[iext]) continue; // Skip unused 
      string filename;
      extensionTable.readCell(filename, "FILENAME", iext);
      /**/if (iext%10==0) cerr << "# Reading object catalog " << iext
			       << "/" << extensions.size()
			       << " from " << filename << endl;
      int hduNumber;
      extensionTable.readCell(hduNumber, "EXTENSION", iext);
      string xKey;
      extensionTable.readCell(xKey, "xKey", iext);
      string yKey;
      extensionTable.readCell(yKey, "yKey", iext);
      string idKey;
      extensionTable.readCell(idKey, "idKey", iext);
      string errKey;
      extensionTable.readCell(errKey, "errKey", iext);
      double weight;
      extensionTable.readCell(weight, "Weight", iext);

      // Relevant structures for this extension
      Extension& extn = *extensions[iext];
      Exposure& expo = *exposures[extn.exposure];

      const SubMap* pixmap=extn.map;
      bool isReference = (expo.instrument == REF_INSTRUMENT);
      bool isTag = (expo.instrument == TAG_INSTRUMENT);

      Wcs* startWcs = extn.startWcs;
      startWcs->reprojectTo(*fieldProjections[expo.field]);

      if (!startWcs) {
	cerr << "Failed to find initial Wcs for device " 
	     << instruments[expo.instrument]->deviceNames.nameOf(extn.device)
	     << " of exposure " << expo.name
	     << endl;
	exit(1);
      }

      bool useRows = stringstuff::nocaseEqual(idKey, "_ROW");
      // What we need to read from the FitsTable:
      vector<string> neededColumns;
      if (!useRows)
	neededColumns.push_back(idKey);
      neededColumns.push_back(xKey);
      neededColumns.push_back(yKey);
      if (!isTag)
	neededColumns.push_back(errKey);

      FTable ff;
#ifdef _OPENMP
#pragma omp critical(fitsio)
#endif
      {
	FITS::FitsTable ft(filename, FITS::ReadOnly, hduNumber);
	ff = ft.extract(0, -1, neededColumns);
      }
      vector<long> id;
      if (useRows) {
	id.resize(ff.nrows());
	for (long i=0; i<id.size(); i++)
	  id[i] = i;
      } else {
	ff.readCells(id, idKey);
      }
      Assert(id.size() == ff.nrows());

      bool errorColumnIsDouble = true;
      try {
	double d;
	if (!isTag)
	  ff.readCell(d, errKey, 0);
      } catch (img::FTableError& e) {
	errorColumnIsDouble = false;
      }

      double sysError = isReference ? referenceSysError : pixSysError;

      for (long irow = 0; irow < ff.nrows(); irow++) {
	map<long, Detection*>::iterator pr = extn.keepers.find(id[irow]);
	if (pr == extn.keepers.end()) continue; // Not a desired object

	// Have a desired object now.  Fill its Detection structure
	Detection* d = pr->second;
	extn.keepers.erase(pr);

	d->map = pixmap;
	ff.readCell(d->xpix, xKey, irow);
	ff.readCell(d->ypix, yKey, irow);
	double sigma;
	if (isTag) {
	  sigma = 0.;  // Don't need sigma for tags
	} else if (errorColumnIsDouble) {
	  ff.readCell(sigma, errKey, irow);
	} else {
	  float f;
	  ff.readCell(f, errKey, irow);
	  sigma = f;
	}

	sigma = std::sqrt(sysError*sysError + sigma*sigma);
	d->sigmaPix = sigma;
	startWcs->toWorld(d->xpix, d->ypix, d->xw, d->yw);
	Matrix22 dwdp = startWcs->dWorlddPix(d->xpix, d->ypix);
	// no clips on tags
	double wt = isTag ? 0. : pow(sigma,-2.);
	d->clipsqx = wt / (dwdp(0,0)*dwdp(0,0)+dwdp(0,1)*dwdp(0,1));
	d->clipsqy = wt / (dwdp(1,0)*dwdp(1,0)+dwdp(1,1)*dwdp(1,1));
	d->wtx = d->clipsqx * weight;
	d->wty = d->clipsqy * weight;
	    
      } // End loop over catalog objects

      if (!extn.keepers.empty()) {
	cerr << "Did not find all desired objects in catalog " << filename
	     << " extension " << hduNumber
	     << endl;
	exit(1);
      }
    } // end loop over catalogs to read
	
    cerr << "Done reading catalogs." << endl;

    // Make a pass through all matches to reserve as needed and purge 
    // those not wanted.  

    {
      ran::UniformDeviate u;
      if (randomNumberSeed > 0)  u.Seed(randomNumberSeed);
      long dcount=0;
      int dof=0;
      double chi=0.;
      double maxdev=0.;

      list<Match*>::iterator im = matches.begin();
      while (im != matches.end() ) {
	// Remove Detections with too-large errors:
	Match::iterator j=(*im)->begin();
	while (j != (*im)->end()) {
	  Detection* d = *j;
	  // Keep it if pixel error is small or if it's from a tag or reference "instrument"
	  if ( d->sigmaPix > maxPixError
	       && exposures[ extensions[d->catalogNumber]->exposure]->instrument >= 0) {
	    j = (*im)->erase(j, true);  // will delete the detection too.
	  } else {
	    ++j;
	  }
	}
	(*im)->countFit(); // Recount the fittable objects after changing wts
	int nFit = (*im)->fitSize();
	if ( nFit < minMatches) {
	  // Remove entire match if it's too small, and kill its Detections too
	  (*im)->clear(true);
	  im = matches.erase(im);
	} else {
	  // Still a good match.

	  dcount += nFit;	// Count total good detections
	  chi += (*im)->chisq(dof, maxdev);

	  // Decide whether to reserve this match
	  if (reserveFraction > 0.)
	    (*im)->setReserved( u < reserveFraction );

	  ++im;
	}
      } // End loop over all matches

      cout << "Using " << matches.size() 
	   << " matches with " << dcount << " total detections." << endl;
      cout << " chisq " << chi << " / " << dof << " dof maxdev " << maxdev << endl;

    } // End Match-culling section.


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
	cout << "Fitting " << mcount << " matches with " << dcount << " detections "
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
      cout << "Final chisq: " << chisq 
	   << " / " << dof << " dof, max deviation " << max
	   << "  new clip threshold at: " << thresh << " sigma"
	   << endl;
      if (thresh >= max || (oldthresh>0. && (1-thresh/oldthresh)<minimumImprovement)) {
	// Sigma clipping is no longer doing much.  Quit if we are at full precision,
	// else require full target precision and initiate another pass.
	if (coarsePasses) {
	  coarsePasses = false;
	  ca.setRelTolerance(chisqTolerance);
	  cout << "--Starting strict tolerance passes, clipping full matches" << endl;
	  oldthresh = thresh;
	  nclip = ca.sigmaClip(thresh, false, true);
	  cout << "Clipped " << nclip
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
	cout << "No clipping--Starting strict tolerance passes, clipping full matches" << endl;
	continue;
      }
      cout << "Clipped " << nclip
	   << " matches " << endl;
      
    } while (coarsePasses || nclip>0);
  
    // The re-fitting is now complete.  Serialize all the fitted coordinate systems
    {
      ofstream ofs(outWcs.c_str());
      if (!ofs) {
	cerr << "Error trying to open output file for fitted Wcs: " << outWcs << endl;
	// *** will not quit before dumping output ***
      } else {
	mapCollection.write(ofs);
      }
    }

    // If there are reserved Matches, run sigma-clipping on them now.
    if (reserveFraction > 0.) {
      cout << "** Clipping reserved matches: " << endl;
      oldthresh = 0.;
      do {
	// Report number of active Matches / Detections in each iteration:
	long int mcount=0;
	long int dcount=0;
	ca.count(mcount, dcount, true, 2);
	double max;
	int dof=0;
	double chisq= ca.chisqDOF(dof, max, true);
	cout << "Clipping " << mcount << " matches with " << dcount << " detections "
	     << " chisq " << chisq << " / " << dof << " dof,  maxdev " << max 
	     << " sigma" << endl;
      
	double thresh = sqrt(chisq/dof) * clipThresh;
	cout << "  new clip threshold: " << thresh << " sigma"
	     << endl;
	if (thresh >= max) break;
	if (oldthresh>0. && (1-thresh/oldthresh)<minimumImprovement) break;
	oldthresh = thresh;
	nclip = ca.sigmaClip(thresh, true, clipEntireMatch);
	cout << "Clipped " << nclip
	     << " matches " << endl;
      } while (nclip>0);
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
	sigpix.push_back(d->sigmaPix);

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


