// Program to fit astrometric solutions to detection catalogs already matched by WCSFoF.
#include <fstream>
#include <sstream>

#include "WCSFit_fit.h"

using namespace std;

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

  // Fractional reduction in RMS required to continue sigma-clipping:
  double minimumImprovement=0.02;
  try {
  // Read all the command-line and parameter-file program parameters
  processParameters(parameters, usage, 1, argc, argv);
  string inputTables = argv[1];

  PMMatch::setPrior(pmPrior, parallaxPrior);

  FitClass fitclass(minMatches=minMatches, verbose=verbose);

  list<string> useInstrumentList = splitArgument(useInstruments);

  // Objects to ignore on input:
  ExtensionObjectSet skipSet(skipFile);
  
  // The list of exposures that are considered valid sources of color information:
  list<string> useColorList = splitArgument(colorExposures);

  // Exposures to ignore:
  list<string> skipExposureList = splitArgument(skipExposures);

  //fitclass.addInputYAML(inputMaps);
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
  // Read the Fields table from input, copy to a new output FITS file, extract needed info  
  fitclass.fields = Fields::read(inputTables, outCatalog, pmEpoch);
    
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
  fitclass.instruments = readInstruments(instrumentHDUs, useInstrumentList, inputTables, outCatalog,
                                         outputCatalogAlreadyOpen);
  
  PROGRESS(1,Reading exposures);

  // This vector will hold the color-priority value of each exposure.
  // -1 means an exposure that does not hold color info.
  vector<int> exposureColorPriorities;
  vector<unique_ptr<Exposure>> exposures = readExposures(fitclass.instruments,
                                     fitclass.fields.epochs(),
                                     exposureColorPriorities,
                                     useColorList,
                                     inputTables, outCatalog, skipExposureList, 
                                     true, // Use reference exposures for astrometry
                                     outputCatalogAlreadyOpen);

  fitclass.setExposures(std::move(exposures), sysError, referenceSysError);

  PROGRESS(1,Reading extensions);

  // Read info about all Extensions - we will keep the Table around.
  FTable extensionTable;
  {
    FITS::FitsTable ft(inputTables, FITS::ReadOnly, "Extensions");
    extensionTable = ft.extract();
    FITS::FitsTable out(outCatalog, FITS::ReadWrite+FITS::Create, "Extensions");
    out.copy(extensionTable);
  }
  
  vector<unique_ptr<ColorExtension>> colorExtensions;
  fitclass.extensions = readExtensions<Astro>(extensionTable,
                                              fitclass.instruments,
                                              fitclass.exposures,
                                              exposureColorPriorities,
                                              colorExtensions,
                                              inputYAML,
                                              //fitclass.inputYAML,
                                              verbose>=1);  // Print reading progress?
  
  
  fitclass.setRefWCSNames();
  
  fitclass.setupMaps(inputYAML, fixMaps=fixMaps);
  
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

    readMatches<Astro>(ff, fitclass.matches, fitclass.extensions, colorExtensions, skipSet,
                       fitclass.minMatches, usePM);
  } // End loop over input matched catalogs

  if (verbose>=0) cout << "# Total match count: " << fitclass.matches.size() << endl;

  // Now loop over all original catalog bintables, reading the desired rows
  // and collecting needed information into the Detection structures
  PROGRESS(1,Reading catalogs);
  readObjects<Astro>(extensionTable, fitclass.exposures, fitclass.extensions, fitclass.fields.projections());

  // Now loop again over all catalogs being used to supply colors,
  // and insert colors into all the Detections they match
  PROGRESS(1,Reading colors);
  readColors<Astro>(extensionTable, colorExtensions);

  fitclass.fit(maxError=maxError, minFitExposures=minFitExposures, reserveFraction=reserveFraction,
               randomNumberSeed=randomNumberSeed, minimumImprovement=minimumImprovement,
               clipThresh=clipThresh, chisqTolerance=chisqTolerance, clipEntireMatch=clipEntireMatch,
               divideInPlace=divideInPlace, purgeOutput=purgeOutput, minColor=minColor, maxColor=maxColor);

  fitclass.saveResults(outWcs, outCatalog, starCatalog);

  }
  catch (std::runtime_error& m) {
    quit(m,1);
  }
}
