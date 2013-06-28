// Program to fit photometric solutions to detection catalogs already matched by WCSFoF.
#include <fstream>
#include <sstream>
#include <map>

#include "Std.h"
#include "Astrometry.h"
#include "FitsTable.h"
#include "StringStuff.h"
#include "Pset.h"
#include "PixelMapCollection.h"
#include "TemplateMap.h"
#include "TPVMap.h"
#include "Random.h"
#include "PhotoMatch.h"
#include "PhotoMapCollection.h"

using namespace std;
using namespace stringstuff;
using namespace photometry;
using img::FTable;

string usage=
  "PhotoFit: Refine photometric solutions for a matched set of catalogs.\n"
  "usage: WCSFit <match file> [parameter file] [parameter file...]\n"
  "      <match file>:  FITS file with binary tables produced by WCSFoF\n"
  "     stdin is read as parameters.";

// Temporary documentation:
// parameter renameInstruments is a string that has format like
// <instrument name>=<map name>, <instrument name> = <map name>, ...
// where the Instrument's PixelMap will be given the map name.  Regex allowed.
// Note that this is assuming that regexes do not include = or , characters.
// whitespace will be stripped from edges of each name.

// parameter existingMaps is a string with
// <mapName>@<file name>, ...
// which says that maps with the mapName should be initially be deserialized from file.
// Regexes ok; same caveats: no @ or , in names, whitespace stripped.

// parameter fixMaps is a string with
// <mapName>, ...
// where any given mapName should have its parameters fixed at initial values during
// the fitting.  Regexes allowed (no commas!).

// parameter canonicalExposures is a string with
// <exposureID>, ....
// where ares are exposures that will be given an identity exposure map.  There must be
// 0 or 1 of these specified for any instrument that has Instrument Map with free parameters.
// Default is to find an exposure that has data in all devices and use it.
// Will have an error if there is more than one constraint on any Instrument.

// Note that PhotoMaps for devices within instrument will get names <instrument>/<device>.
// And PhotoMaps for individual exposures will get names <exposure>/<device>.

#include "PhotoInstrument.h"

// A helper function that strips white space from front/back of a string and replaces
// internal white space with underscores:
void
spaceReplace(string& s) {
  stripWhite(s);
  s = regexReplace("[[:space:]]+","_",s);
}

// A function to parse strings describing PixelMap models, setting tolerance for inversion of world coords:
PhotoMap* photoMapDecode(string code, string name, PolyMap::ArgumentType arg=PolyMap::Device);

class Accum {
public:
  double sum_m;
  double sum_mw;
  double sum_mm;
  double sum_mmw;
  double sum_w;
  int n;
  double sum_x;
  double sum_y;
  double sum_dof;
  Accum(): sum_m(0.), sum_mw(0.),
	   sum_mm(0.), sum_mmw(0.),
	   sum_w(0.),
	   sum_dof(0.), 
	   sum_x(0.), sum_y(0.),
	   n(0) {}
  void add(const Detection* d, double magoff=0., double dof=1.) {
    double dm = d->magOut - magoff;
    sum_m += dm;
    sum_mw += dm * d->wt;
    sum_mm += dm * dm;
    sum_mmw += dm * dm *d->wt;
    sum_w += d->wt;
    sum_dof += dof;
    sum_x += d->args.xExposure;
    sum_y += d->args.yExposure;
    ++n;
    }
  double rms() const {
    return n > 0 ? sqrt( sum_mm/n ) : 0.;
  }
  double reducedChisq() const {
    return sum_dof>0 ? sum_mmw/sum_dof : 0.;
  }
  string summary() const {
    ostringstream oss;
    double dm = 0., sigm=0.;
    double x=0., y=0.;
    if (n>0) {
      dm = sum_mw / sum_w;
      sigm = 1./sqrt(sum_w);
      x = sum_x / n;
      y = sum_y / n;
    }
    oss << setw(4) << n 
	<< fixed << setprecision(1)
	<< " " << setw(6) << sum_dof
	<< " " << setw(5) << dm*1000.
	<< " " << setw(5) << sigm*1000.
	<< " " << setw(5) << rms()*1000.
	<< setprecision(2) 
	<< " " << setw(5) << reducedChisq()
	<< setprecision(5) << showpoint << showpos
	<< " " << setw(9) << x 
	<< " " << setw(9) << y ;
    return oss.str();
  }
};

int
main(int argc, char *argv[])
{
  double reserveFraction;
  double clipThresh;
  double priorClipThresh;
  bool clipEntirePrior;
  double maxMagError;
  double magSysError;
  double referenceSysError;
  int minMatches;
  bool clipEntireMatch;
  bool exposureModel;
  string deviceModel;
  string outCatalog;
  string catalogSuffix;
  string newHeadSuffix;
  double chisqTolerance;
  string renameInstruments;
  string useInstruments;
  string colorExposures;
  string existingMaps;
  string fixMaps;
  string canonicalExposures;
  string outPhotFile;
  string priorFile;
  string outPriorFile;
  Pset parameters;
  {
    const int def=PsetMember::hasDefault;
    const int low=PsetMember::hasLowerBound;
    const int up=PsetMember::hasUpperBound;
    const int lowopen = low | PsetMember::openLowerBound;
    const int upopen = up | PsetMember::openUpperBound;

    parameters.addMemberNoValue("INPUTS");
    parameters.addMember("maxMagError",&maxMagError, def | lowopen,
			 "cut objects with magnitude uncertainty above this", 0.05, 0.);
    parameters.addMember("magSysError",&magSysError, def | low,
			 "additional systematic error for all detected mags", 0.002, 0.);
    parameters.addMember("referenceSysError",&referenceSysError, def | low,
			 "Reference object systematic error (mag)", 0.05, 0.);
    parameters.addMember("minMatch",&minMatches, def | low,
			 "minimum number of detections for usable match", 2, 2);
    parameters.addMember("useInstruments",&useInstruments, def,
			 "the instruments that are assumed to produce matching mags","");
    parameters.addMember("colorExposures",&colorExposures, def,
			 "exposures holding valid colors for stars","");
    parameters.addMemberNoValue("CLIPPING");
    parameters.addMember("clipThresh",&clipThresh, def | low,
			 "Clipping threshold (sigma)", 5., 2.);
    parameters.addMember("clipEntireMatch",&clipEntireMatch, def,
			 "Discard entire object if one outlier", false);
    parameters.addMember("priorClipThresh",&priorClipThresh, def | low,
			 "Clipping threshold (sigma)", 5., 2.);
    parameters.addMember("clipEntirePrior",&clipEntirePrior, def,
			 "Discard entire night's prior if one outlier", false);
    parameters.addMemberNoValue("FITTING");
    parameters.addMember("priorFile", &priorFile, def,
			 "File specifying any priors to apply to zeropoints and colors", "");
    parameters.addMember("reserveFraction",&reserveFraction, def | low,
			 "Fraction of matches reserved from re-fit", 0., 0.);
    parameters.addMember("exposureModel",&exposureModel, def,
			 "Form of per-exposure map", "Constant");
    parameters.addMember("deviceModel",&deviceModel, def,
			 "Form of device map", "Poly 2");
    parameters.addMember("chisqTolerance",&chisqTolerance, def | lowopen,
			 "Fractional change in chisq for convergence", 0.001, 0.);
    parameters.addMember("renameInstruments",&renameInstruments, def,
			 "list of new names to give to instruments for maps","");
    parameters.addMember("existingMaps",&existingMaps, def,
			 "list of photometric maps to draw from existing files","");
    parameters.addMember("fixMaps",&fixMaps, def,
			 "list of map components or instruments to hold fixed","");
    parameters.addMember("canonicalExposures",&canonicalExposures, def,
			 "list of exposures that will have identity exposure maps","");
    parameters.addMemberNoValue("OUTPUTS");
    parameters.addMember("outCatalog",&outCatalog, def,
			 "Output FITS binary catalog", "wcscat.fits");
    parameters.addMember("outPhotFile",&outPhotFile, def,
			 "Output serialized photometric solutions", "photfit.phot");
    parameters.addMember("outPriorFile",&outPriorFile, def,
			 "Output listing of zeropoints etc of exposures tied by priors","");
  }

  // Fractional reduction in RMS required to continue sigma-clipping:
  const double minimumImprovement=0.02;
  const int REF_INSTRUMENT=-1;	// Instrument for reference objects (no fitting) ??? color terms???
  const int TAG_INSTRUMENT=-2;	// Exposure number for tag objects (no fitting nor contrib to stats)

  const string stellarAffinity="STELLAR";

  try {
    
    // Read parameters
    if (argc<2) {
      cerr << usage << endl;
      exit(1);
    }
    string inputTables = argv[1];

    parameters.setDefault();
    for (int i=2; i<argc; i++) {
      // Open & read all specified input files
      ifstream ifs(argv[i]);
      if (!ifs) {
	cerr << "Can't open parameter file " << argv[i] << endl;
	exit(1);
      }
      try {
	parameters.setStream(ifs);
      } catch (std::runtime_error &m) {
	cerr << "In file " << argv[i] << ":" << endl;
	quit(m,1);
      }
    }
    // and stdin:
    try {
      parameters.setStream(cin);
    } catch (std::runtime_error &m) {
      cerr << "In stdin:" << endl;
      quit(m,1);
    }

    // List parameters in use
    parameters.dump(cout);

    /////////////////////////////////////////////////////
    // Parse all the parameters describing maps etc. 
    /////////////////////////////////////////////////////

    // Teach PhotoMapCollection about new kinds of PhotoMaps:
    // Currently all map types are automatically loaded in PhotoMapCollection constructor,
    // but if new ones are invented, this is the format:
    // PhotoMapCollection::registerMapType<PolyMap>();

    const char listSeperator=',';

    // First is a regex map from instrument names to the names of their PhotoMaps
    RegexReplacements instrumentTranslator;
    {
      list<string> ls = split(renameInstruments, listSeperator);
      for (list<string>::const_iterator i = ls.begin();
	   i != ls.end();
	   ++i) {
	if (i->empty()) continue;
	list<string> ls2 = split(*i, '=');
	if (ls2.size() != 2) {
	  cerr << "renameInstruments has bad translation spec: " << *i << endl;
	  exit(1);
	}
	string regex = ls2.front();
	string replacement = ls2.back();
	stripWhite(regex);
	stripWhite(replacement);
	instrumentTranslator.addRegex(regex, replacement);
      }
    }

    // Next is regex map from PhotoMap or instrument map names to files 
    // from which we should de-serialize them
    RegexReplacements existingMapFinder;
    {
      list<string> ls = split(existingMaps, listSeperator);
      for (list<string>::const_iterator i = ls.begin();
	   i != ls.end();
	   ++i) {
	if (i->empty()) continue;
	list<string> ls2 = split(*i, '@');
	if (ls2.size() != 2) {
	  cerr << "existingMaps has bad spec: " << *i << endl;
	  exit(1);
	}
	string regex = ls2.front();
	string replacement = ls2.back();
	stripWhite(regex);
	stripWhite(replacement);
	existingMapFinder.addRegex(regex, replacement);
      }
    }

    // This is list of regexes of PhotoMap names (or instrument names) that should
    // have their parameters held fixed.
    list<string> fixMapList = split(fixMaps, listSeperator);
    for (list<string>::iterator i = fixMapList.begin();
	 i != fixMapList.end(); ) 
      if (i->empty()) {
	i = fixMapList.erase(i);
      } else {
	stripWhite(*i);
	++i;
      }

    list<string> canonicalExposureList = split(canonicalExposures, listSeperator);
    for (list<string>::iterator i = canonicalExposureList.begin();
	 i != canonicalExposureList.end(); ) {
      stripWhite(*i);
      if (i->empty()) {
	i = canonicalMapList.erase(i);
      } else {
	++i;
      }
    }

    // The list of instruments that we will be matching together in this run:
    // have their parameters held fixed.
    list<string> useInstrumentList = split(useInstruments, listSeperator);
    for (list<string>::iterator i = useInstrumentList.begin();
	 i != useInstrumentList.end(); )  {
      stripWhite(*i);
      if (i->empty()) {
	i = useInstrumentList.erase(i);
      } else {
	++i;
      }
    }

    // The list of exposures that are considered valid sources of color information:
    list<string> useColorList = split(colorExposures, listSeperator);
    for (list<string>::iterator i = useColorList.begin();
	 i != useColorList.end(); ) {
      stripWhite(*i);
      if (i->empty()) {
	i = useColorList.erase(i);
      } else {
	++i;
      }
    }

    /////////////////////////////////////////////////////
    //  Read in properties of all Instruments, Devices, Exposures
    /////////////////////////////////////////////////////

    // All the names will be stripped of leading/trailing white space, and internal
    // white space replaced with single underscore - this keeps PhotoMap parsing functional.

    // Let's figure out which of our FITS extensions are Instrument or MatchCatalog
    vector<int> instrumentHDUs;
    vector<int> catalogHDUs;
    {
      // Find which extensions are instrument tables
      FITS::FitsFile ff(inputTables);
      for (int i=1; i<ff.HDUCount(); i++) {
	FITS::Hdu h(inputTables, FITS::HDUAny, i);
	if (stringstuff::nocaseEqual(h.getName(), "Instrument")) {
	  instrumentHDUs.push_back(i);
	} else if (stringstuff::nocaseEqual(h.getName(), "MatchCatalog")) {
	  catalogHDUs.push_back(i);
	}
      }
    }
    
    // Read in all the instrument extensions and their device info.
    vector<Instrument*> instruments(instrumentHDUs.size(),0);
    for (int i=0; i<instrumentHDUs.size(); i++) {
      FITS::FitsTable ft(inputTables, FITS::ReadOnly, instrumentHDUs[i]);
      // Append new table to output:
      Assert( stringstuff::nocaseEqual(ft.getName(), "Instrument"));
      FTable ff=ft.extract();

      string instrumentName;
      int instrumentNumber;
      if (!ff.header()->getValue("Name", instrumentName)
	  || !ff.header()->getValue("Number", instrumentNumber)) {
	cerr << "Could not read name and/or number of instrument at extension "
	     << instrumentHDUs[i] << endl;
      }
      spaceReplace(instrumentName);
      // Is this an instrument we are going to care about?
      string name = instrumentName;
      instrumentTranslator(name);  // Apply any name remapping specified.
      if (!regexMatchAny(useInstrumentList, name)) 
	continue;	// Skip this instrument if we aren't using its mags

      FITS::FitsTable out(outCatalog, FITS::ReadWrite+FITS::Create, -1);
      out.setName("Instrument");
      out.adopt(ff);
      Assert(instrumentNumber < instruments.size());
      Instrument* instptr = new Instrument(instrumentName);
      instruments[instrumentNumber] = instptr;
      
      vector<string> devnames;
      vector<double> vxmin;
      vector<double> vxmax;
      vector<double> vymin;
      vector<double> vymax;
      ff.readCells(devnames, "Name");
      ff.readCells(vxmin, "XMin");
      ff.readCells(vxmax, "XMax");
      ff.readCells(vymin, "YMin");
      ff.readCells(vymax, "YMax");
      for (int j=0; j<devnames.size(); j++) {
	spaceReplace(devnames[j]);
	instptr->addDevice( devnames[j],
			    Bounds<double>( vxmin[j], vxmax[j], vymin[j], vymax[j]));
      }
    }

    // Read in the table of exposures
    vector<Exposure*> exposures;
    // This vector will hold the color-priority value of each exposure.  -1 means an exposure
    // that does not hold color info.
    vector<int> exposureColorPriorities;
    {
      FITS::FitsTable ft(inputTables, FITS::ReadOnly, "Exposures");
      FITS::FitsTable out(outCatalog, FITS::ReadWrite+FITS::Create, "Exposures");
      FTable ff = ft.extract();
      out.adopt(ff);
      vector<string> names;
      vector<double> ra;
      vector<double> dec;
      vector<int> fieldNumber;
      vector<int> instrumentNumber;
      ff.readCells(names, "Name");
      ff.readCells(ra, "RA");
      ff.readCells(dec, "Dec");
      ff.readCells(fieldNumber, "fieldNumber");
      ff.readCells(instrumentNumber, "InstrumentNumber");
      exposureColorPriorities = vector<int>(names.size(), -1);
      for (int i=0; i<names.size(); i++) {
	// The projection we will use for this exposure:
	astrometry::Gnomonic gn(astrometry::Orientation(astrometry::SphericalICRS(ra[i]*DEGREE,
										  dec[i]*DEGREE)));
	spaceReplace(names[i]);
	Exposure* expo = new Exposure(names[i],gn);
	expo->field = fieldNumber[i];
	expo->instrument = instrumentNumber[i];

	// See if this exposure name matches any of the color exposures
	int priority = 0;
	for (list<string>::const_iterator j = useColorList.begin();
	     j != useColorList.end();
	     ++j, ++priority) {
	  if (regexMatch(*j, names[i])) {
	    // Yes: give this exposure a priority according to order of the exposure it first matches.
	    exposureColorPriorities[i] = priority;
	    continue;
	  }
	}
	// Is this exposure in an instrument of interest, or a reference or tag?
	if (expo->instrument>=0 && !instruments[expo->instrument]) {
	  // No - don't save info, insert null into exposure list
	  // ??? Need a way to specify whether to use reference exposures, tags ???
	  delete expo;
	  expo = 0;
	}
	exposures.push_back(expo);
      }
    }

    // Read info about all Extensions - we will keep the Table around.
    FTable extensionTable;
    vector<Extension*> extensions;
    vector<ColorExtension*> colorExtensions;
    {
      FITS::FitsTable ft(inputTables, FITS::ReadOnly, "Extensions");
      extensionTable = ft.extract();
      FITS::FitsTable out(outCatalog, FITS::ReadWrite+FITS::Create, "Extensions");
      out.copy(extensionTable);
    }
    for (int i=0; i<extensionTable.nrows(); i++) {
      Extension* extn = new Extension;
      int j;
      extensionTable.readCell(j, "ExposureNumber", i);

      // Determine whether this extension might be used to provide colors
      int colorPriority = exposureColorPriorities[j];
      if (colorPriority < 0) {
	colorExtensions.push_back(0);
      } else {
	ColorExtension* ce = new ColorExtension;
	ce->priority = colorPriority;
	colorExtensions.push_back(ce);
      }
	
      if (!exposures[j]) {
	// This extension is not in an exposure of interest.  Skip it.
	delete extn;
	extn = 0;
	extensions.push_back(extn);
	continue;
      }
      extensions.push_back(extn);
      extn->exposure = j;
      extensionTable.readCell(j, "DeviceNumber", i);
      extn->device = j;
      string s;
      extensionTable.readCell(s, "WCS", i);
      if (stringstuff::nocaseEqual(s, "ICRS")) {
	// Create a Wcs that just takes input as RA and Dec in degrees;
	astrometry::IdentityMap identity;
	astrometry::SphericalICRS icrs;
	extn->startWcs = new astrometry::Wcs(&identity, icrs, "ICRS_degrees", DEGREE);
      } else {
	istringstream iss(s);
	astrometry::PixelMapCollection pmcTemp;
	if (!pmcTemp.read(iss)) {
	  cerr << "Did not find Wcs format for ???" << endl;
	  exit(1);
	}
	string wcsName = pmcTemp.allWcsNames().front();
	extn->startWcs = pmcTemp.cloneWcs(wcsName);
      }
      // destination projection is the Exposure projection, whose coords used for any
      // exposure-level magnitude corrections
      extn->startWcs->reprojectTo(*exposures[extn->exposure]->projection);
      // ??? or use the field center for reference/tag exposures ???
    }

    /////////////////////////////////////////////////////
    //  Create and initialize all magnitude maps
    /////////////////////////////////////////////////////

    // Logic now for initializing instrument maps:
    // 1) Assign a name to each instrument map; either it's name of instrument,
    //    or we have a mapping from input parameters
    // 2) If the instrument name is something we were told to reuse:
    //        * get its form and values from the input map libraries for each device
    // 3) If not reusing, make new maps for each device.
    //        * Select canonical exposure either from cmd line or finding one
    //        * set reference exposure's map to identity
    //        * Initialize new map to identity.
    // 4) Set initial exposure maps for each exposure that is not already canonical.
    // 5) Freeze parameters of all map elements specified in input.

    // Create a PhotoMapCollection to hold all the components of the PhotoMaps.
    PhotoMapCollection mapCollection;
    // Create an identity map used by all reference catalogs
    mapCollection.learnMap(IdentityMap());

    // This is list of PhotoMaps whose parameters will be held fixed during fitting:
    list<string> fixedMapNames;

    // Set up coordinate maps for each instrument & exposure
    for (int iinst=0; iinst<instruments.size(); iinst++) {
      Instrument* inst = instruments[iinst];
      if (!inst) continue; // Only using some instruments

      // (1) choose a name to be used for this Instrument's PhotoMap:
      string mapName = inst->name;
      instrumentTranslator(mapName);  // Apply any remapping specified.

      // (2) See if this instrument has a map that we are re-using, and if params are fixed
      string loadFile = mapName;
      bool loadExistingMaps = existingMapFinder(loadFile);
      bool fixInstrumentParameters = regexMatchAny(fixMapList, mapName);

      // (3) decide for all devices whether they are looked up, or need creating,
      //     and whether they are fixed
      vector<bool> deviceMapsExist(inst->nDevices, loadExistingMaps);
      vector<bool> fixDeviceMaps(inst->nDevices, fixInstrumentParameters);

      PhotoMapCollection pmcInstrument;
      if (loadExistingMaps) {
	ifstream ifs(loadFile.c_str());
	if (!ifs) {
	  cerr << "Could not open file " + loadFile + " holding existing PhotoMaps" << endl;
	  exit(1);
	}
	pmcInstrument.read(ifs);
      }

      for (int idev=0; idev<inst->nDevices; idev++) {
	string devMapName = mapName + '/' + inst->deviceNames.nameOf(idev);
	inst->mapNames[idev] = devMapName;
	// A PhotoMap file for Device overrides one for the Instrument
	string devLoadFile = devMapName;
	if (existingMapFinder(devLoadFile)) {
	  deviceMapsExist[idev] = true;
	} else if (deviceMapsExist[idev]) {
	  // Already have a map file from the Instrument
	  devLoadFile = loadFile;
	}
	if (regexMatchAny(fixMapList, devMapName))
	  fixDeviceMaps[idev] = true;

	if (fixDeviceMaps[idev] && !deviceMapsExist[idev]) {
	  cerr << "WARNING: Requested fixed parameters for Device map "
	       << devMapName
	       << " that does not have initial values.  Ignoring request to fix params."
	       << endl;
	}

	// Done with this device for now if we have no existing map to read.
	if (!deviceMapsExist[idev]) continue;
	  
	// Load the PhotoMap for this device if it is to be extracted from a file.
	// We will not be throwing exceptions for duplicate names
	PhotoMap* deviceMap;
	if (devLoadFile==loadFile) {
	  // Get map from the Instrument's collection:
	  deviceMap = pmcInstrument.issueMap(devMapName);
	  mapCollection.learnMap(*deviceMap);
	} else {
	  // Device has its own serialized file:
	  PhotoMapCollection pmcDev;
	  ifstream ifs(devLoadFile.c_str());
	  if (!ifs) {
	    cerr << "Could not open file " + devLoadFile + " holding existing PhotoMaps" << endl;
	    exit(1);
	  }
	  pmcDev.read(ifs);
	  deviceMap = pmcDev.issueMap(devMapName);
	  mapCollection.learnMap(*deviceMap);
	}

	if (fixDeviceMaps[idev])
	  fixedMapNames.push_back(devMapName);
      } // end of device loop

      // Canonical exposure will have its Exposure map will be fixed at identity 
      // if none of the Device's maps are
      // fixed at initial values;  we assume that any single device being fixed will
      // break degeneracy between Exposure and Instrument models.

      // First find all exposures from this instrument,
      // and find if one is canonical;
      long canonicalExposure = -1;
      list<long> exposuresWithInstrument;
      for (long iexp = 0; iexp < exposures.size(); iexp++) {
	if (!exposures[iexp]) continue;  // Not using this exposure or its instrument
	if (exposures[iexp]->instrument == iinst) {
	  exposuresWithInstrument.push_back(iexp);
	  if (regexMatchAny(canonicalExposureList, exposures[iexp]->name)) {
	    if (canonicalExposure < 0) {
	      // This is our canonical exposure
	      canonicalExposure = iexp;
	    } else {
	      // Duplicate canonical exposures is an error
	      cerr << "More than one canonical exposure for instrument " 
		   << inst->name << ": " << endl;
	      cerr << exposures[canonicalExposure]->name
		   << " and " << exposures[iexp]->name
		   << endl;
	      exit(1);
	    }
	  }
	}
      } // end exposure loop

      // Do we need a canonical exposure? Only if we have an an uninitialized device map
      bool needCanonical = false;
      for (int idev=0; idev<deviceMapsExist.size(); idev++)
	if (!deviceMapsExist[idev]) {
	  needCanonical = true;
	  break;
	}

      // Or if none of the devices have their parameters fixed
      bool noDevicesFixed = true;
      for (int idev=0; idev<deviceMapsExist.size(); idev++)
	if (deviceMapsExist[idev] && fixDeviceMaps[idev]) {
	  noDevicesFixed = false;
	  break;
	}

      if (noDevicesFixed) needCanonical = true;

      if (needCanonical && canonicalExposure > 0) {
	// Make sure our canonical exposure has all devices:
	set<long> vexp;
	vexp.insert(canonicalExposure);
	FTable exts = extensionTable.extractRows("ExposureNumber", vexp);
	if (exts.nrows() != inst->nDevices) {
	  cerr << "Canonical exposure " << exposures[canonicalExposure]->name
	       << " for Instrument " << inst->name
	       << " only has " << exts.nrows()
	       << " devices out of " << inst->nDevices
	       << endl;
	  exit(1);
	}
      }

      if (needCanonical && canonicalExposure < 0) {
	// Need canonical but don't have one.  
	// Find an exposure that has all devices for this Instrument.
	for (list<long>::const_iterator i = exposuresWithInstrument.begin();
	     i != exposuresWithInstrument.end();
	     ++i) {
	  // Filter the extension table for this exposure:
	  set<long> vexp;
	  vexp.insert(*i);
	  FTable exts = extensionTable.extractRows("ExposureNumber", vexp);
	  // Stop here if this exposure has right number of extensions (exactly 1 per device):
	  if (exts.nrows() == inst->nDevices) {
	    canonicalExposure = *i;
	    break; 
	  }
	}

	if (canonicalExposure < 0) {
	  cerr << "Could not find exposure with all devices for intrument "
	       << inst->name << endl;
	  exit(1);
	}

      } // end finding a canonical exposure for the Instrument

      // Now we create new PhotoMaps for each Device that does not already have one.
      for (int idev=0; idev < inst->nDevices; idev++) {
	if (deviceMapsExist[idev]) continue;

	// Create a new PhotoMap for this device. 
	PhotoMap* pm=0;
	list<string> pmCodes = stringstuff::split(deviceModel,'+');
	if (pmCodes.size() == 1) {
	  pm = photoMapDecode(pmCodes.front(),inst->mapNames[idev]);
	} else if (pmCodes.size() > 1) {
	  // Build a compound SubMap with desired maps in order:
	  int index = 0;
	  list<PhotoMap*> pmList;
	  for (list<string>::const_iterator ipm = pmCodes.begin();
	       ipm != pmCodes.end();
	       ++ipm, ++index) {
	    ostringstream oss;
	    // Components will have index number appended to device map name:
	    oss << inst->mapNames[idev] << "_" << index;
	    pmList.push_back(photoMapDecode(*ipm, oss.str()));
	  }
	  pm = new SubMap(pmList, inst->mapNames[idev]);
	  // Delete the original components, which have now been duplicated in SubMap:
	  for (list<PhotoMap*>::iterator ipm=pmList.begin();
	       ipm != pmList.end();
	       ++ipm)
	    delete *ipm;
	} else {
	  cerr << "Empty deviceModel specification" << endl;
	  exit(1);
	}

	// Save this map into the collection
	mapCollection.learnMap(*pm);
	delete pm;

      } // end loop creating PhotoMaps for all Devices.

      // Now create an exposure map for all exposures with this instrument,
      // and set initial parameters to be identity ???

      for (list<long>::const_iterator i= exposuresWithInstrument.begin();
	   i != exposuresWithInstrument.end();
	   ++i) {
	int iexp = *i;
	Exposure& expo = *exposures[iexp];
	PhotoMap* expoMap=0;
	// First we check to see if we are going to use an existing exposure map:
	string loadFile = expo.name;
	if (existingMapFinder(loadFile)) {
	  // Adopt the existing map for this exposure:
	  PhotoMapCollection pmcExpo;
	  ifstream ifs(loadFile.c_str());
	  if (!ifs) {
	    cerr << "Could not open serialized map file " 
		 << loadFile
		 << endl;
	    exit(1);
	  }
	  pmcExpo.read(ifs);
	  expoMap = pmcExpo.issueMap(expo.name);
	  mapCollection.learnMap(*expoMap);
	  expo.mapName = expoMap->getName();
	  // If this imported map is to be held fixed, add it to the list:
	  if (regexMatchAny(fixMapList, expo.mapName))
	    fixedMapNames.push_back(expo.mapName);

	  // We're done with this exposure if we have pre-set map for it.
	  continue;
	}

	if (iexp == canonicalExposure && noDevicesFixed) {
	  /**/cerr << "Giving identity map to canonical exposure " << expo.name <<endl;
	  // Give this canonical exposure identity map to avoid degeneracy with Instrument
	  expo.mapName = IdentityMap().getName();
	  continue;
	}
	
	// We will create a new exposure map 
	// ??? more general exposure map now ???  Initialize to identity ???
	expoMap = new PolyMap(exposureOrder, 
			      PolyMap::Exposure,	// Use exposure coords as polynomial args
			      expo.name);
	expo.mapName = expo.name;

	// Insert the Exposure map into the collection
	mapCollection.learnMap(*expoMap);

	// Check for fixed maps
	if (regexMatchAny(fixMapList, expoMap->getName())) {
	  cerr << "WARNING: Requested fixed parameters for Exposure map "
	       << expoMap->getName()
	       << " that does not have initial values.  Ignoring request to fix params."
	       << endl;
	}

	delete expoMap;
      } // end creation of this exposure map

    } // End instrument loop

    // Freeze the parameters of all the fixed maps
    mapCollection.setFixed(fixedMapNames, true);

    // Now construct a SubMap for every extension
    for (int iext=0; iext<extensions.size(); iext++) {
      Extension& extn = *extensions[iext];
      Exposure& expo = *exposures[extn.exposure];
      int ifield = expo.field;
      if ( expo.instrument < 0) {
	// Tag & reference exposures have no instruments and no fitting
	// being done. 
	// ??? color term in reference exposures ???
	extn.map = mapCollection.issueMap(IdentityMap().getName());
      } else {
	// Real instrument, make a map combining its exposure with its Device map:
	string mapName = expo.name + "/" 
	  + instruments[expo.instrument]->deviceNames.nameOf(extn.device);
	list<string> mapElements;
	// The extension map is the device map:
	mapElements.push_back(instruments[expo.instrument]->mapNames[extn.device]);
	// Followed by the exposure map:
	mapElements.push_back(expo.mapName);
	mapCollection.defineChain( mapName, mapElements);

	extn.map = mapCollection.issueMap(mapName);
      }
    } // Extension loop

    /**/cerr << "Total number of free map elements " << mapCollection.nFreeMaps()
	     << " with " << mapCollection.nParams() << " free parameters."
	     << endl;


    //////////////////////////////////////////////////////////
    // Read in all the data
    //////////////////////////////////////////////////////////

    // List of all Matches - they will hold pointers to all Detections too.
    list<Match*> matches;    

    // Start by reading all matched catalogs, creating Detection and Match arrays, and 
    // telling each Extension which objects it should retrieve from its catalog

    for (int icat = 0; icat < catalogHDUs.size(); icat++) {
      /**/cerr << "# Reading catalog extension " << catalogHDUs[icat] << endl;
      FITS::FitsTable ft(inputTables, FITS::ReadOnly, catalogHDUs[icat]);
      FTable ff = ft.use();
      {
	// Only do photometric fitting on the stellar objects:
	string affinity;
	if (!ff.getHdrValue("Affinity", affinity)) {
	  cerr << "Could not find affinity keyword in header of extension " 
	       << catalogHDUs[icat] 
	       << endl;
	  exit(1);
	}
	if (!stringstuff::nocaseEqual(affinity, stellarAffinity))
	  continue;
      }
      vector<int> seq;
      vector<LONGLONG> extn;
      vector<LONGLONG> obj;
      ff.readCells(seq, "SequenceNumber");
      ff.readCells(extn, "Extension");
      ff.readCells(obj, "Object");
      // Smaller collections for each match
      vector<long> extns;
      vector<long> objs;
      // These variables determine what the highest-priority color information available
      // in the match is so far.
      long matchColorExtension = -1;
      long matchColorObject = 0;
      int colorPriority = -1;
      for (int i=0; i<=seq.size(); i++) {
	if (i>=seq.size() || seq[i]==0) {
	  // Processing the previous (or final) match.
	  int nValid = extns.size();
	  if (matchColorExtension < 0) {
	    // There is no color information for this match.  So discard any detection which requires
	    // a color to produce its map:
	    nValid = 0;
	    for (int j=0; j<extns.size(); j++) {
	      Assert(extensions[extns[j]]);
	      if (extensions[extns[j]]->map->needsColor()) {
		// Mark this detection as useless
		extns[j] = -1;
	      } else {
		nValid++;
	      }
	    }
	  }
	  if (nValid >= minMatches) {
	    // Make a match from the valid entries, and note need to get data for the detections and color
	    int j=0;
	    while (extns[j]<0 && j<extns.size()) ++j;  // Skip detections starved of their color
	    Assert(j>=extns.size());
	    Detection* d = new Detection;
	    d->catalogNumber = extns[j];
	    d->objectNumber = objs[j];
	    matches.push_back(new Match(d));
	    extensions[extns[j]]->keepers.insert(std::pair<long, Detection*>(objs[j], d));
	    for (++j; j<extns.size(); j++) {
	      if (extns[j]<0) continue;  // Skip detections needing unavailable color
	      d = new Detection;
	      d->catalogNumber = extns[j];
	      d->objectNumber = objs[j];
	      matches.back()->add(d);
	      extensions[extns[j]]->keepers.insert(std::pair<long, Detection*>(objs[j], d));
	    }
	    if (matchColorExtension >=0) {
	      // Tell the color catalog that it needs to look this guy up:
	      Assert(colorExtensions[matchColorExtension]);
	      colorExtensions[matchColorExtension]->keepers.insert(std::pair<long,Match*>(objs[j],
											  matches.back()));
	    }
	  }
	  // Clear out previous Match:
	  extns.clear();
	  objs.clear();
	  matchColorExtension = -1;
	  colorPriority = -1;
	} // Finished processing previous match

	// If we done reading entries, quit this loop
	if (i >= seq.size()) break;

	// Read in next detection in the catalog
	if (extn[i]<0 || extensions[extn[i]]) {
	  // Record a Detection in a useful extension:
	  extns.push_back(extn[i]);
	  objs.push_back(obj[i]);
	}
      } // End loop of catalog entries
      
    } // End loop over input matched catalogs

    // Now loop over all original catalog bintables, reading the desired rows
    // and collecting needed information into the Detection structures
    for (int iext = 0; iext < extensions.size(); iext++) {
      if (!extensions[iext]) continue; // Skip unused 
      string filename;
      extensionTable.readCell(filename, "Filename", iext);
      /**/cerr << "# Reading object catalog " << iext
	       << "/" << extensions.size()
	       << " from " << filename << endl;
      int hduNumber;
      extensionTable.readCell(hduNumber, "HDUNumber", iext);
      string xKey;
      extensionTable.readCell(xKey, "xKey", iext);
      string yKey;
      extensionTable.readCell(yKey, "yKey", iext);
      string idKey;
      extensionTable.readCell(idKey, "idKey", iext);
      string magKey;
      extensionTable.readCell(magKey, "magKey", iext);
      string magErrKey;
      extensionTable.readCell(magErrKey, "magErrKey", iext);
      double weight;
      extensionTable.readCell(weight, "magWeight", iext);
      // Assume an airmass column is in the extension table.
      // Values <1 mean (including a default 0) mean it's not available.
      double airmass;
      extensionTable.readCell(airmass, "Airmass", iext);

      // Relevant structures for this extension
      Extension& extn = *extensions[iext];
      Exposure& expo = *exposures[extn.exposure];

      const SubMap* sm=extn.map;
      bool isReference = (expo.instrument == REF_INSTRUMENT);
      bool isTag = (expo.instrument == TAG_INSTRUMENT);

      astrometry::Wcs* startWcs = extn.startWcs;

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
      neededColumns.push_back(magKey);
      neededColumns.push_back(magErrKey);

      FITS::FitsTable ft(filename, FITS::ReadOnly, hduNumber);
      FTable ff = ft.extract(0, -1, neededColumns);
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
	ff.readCell(d, magErrKey, 0);
      } catch (img::FTableError& e) {
	errorColumnIsDouble = false;
      }

      double sysError = isReference ? referenceSysError : magSysError;

      for (long irow = 0; irow < ff.nrows(); irow++) {
	map<long, Detection*>::iterator pr = extn.keepers.find(id[irow]);
	if (pr == extn.keepers.end()) continue; // Not a desired object

	// Have a desired object now.  Fill its Detection structure
	Detection* d = pr->second;
	extn.keepers.erase(pr);
	d->map = sm;

	// Read the PhotometryArguments for this Detection:
	ff.readCell(d->args.xDevice, xKey, irow);
	ff.readCell(d->args.yDevice, yKey, irow);
	startWcs->toWorld(d->args.xDevice, d->args.yDevice,
			  d->args.xExposure, d->args.yExposure);

	if (useColor) ff.readCell(d->args.color, colorKey, irow);
	d->args.airmass = airmass;
	
	// Get the mag input and its error
	ff.readCell(d->magIn, magKey, irow);
	double sigma;
	if (errorColumnIsDouble) {
	  ff.readCell(sigma, magErrKey, irow);
	} else {
	  float f;
	  ff.readCell(f, magErrKey, irow);
	  sigma = f;
	}

	// Map to output and estimate output error
	// ??? problem for nonlinear mag functions that may not be initialized here ???
	d->magOut = sm->forward(d->magIn, d->args);

	sigma = std::sqrt(sysError*sysError + sigma*sigma);
	d->sigmaMag = sigma;
	sigma *= sm->derivative(d->magIn, d->args);
	// no clips on tags
	double wt = isTag ? 0. : pow(sigma,-2.);
	d->clipsq = wt;
	d->wt = wt * weight;
	    
      } // End loop over catalog objects

      if (!extn.keepers.empty()) {
	cerr << "Did not find all desired objects in catalog " << filename
	     << " extension " << hduNumber
	     << endl;
	exit(1);
      }
    } // end loop over catalogs to read
	
    cerr << "Done reading catalogs for magnitudes." << endl;

    // Now loop again over all catalogs being used to supply colors,
    // and insert colors into all the PhotoArguments for Detections they match
    for (int iext = 0; iext < colorExtensions.size(); iext++) {
      if (!colorExtensions[iext]) continue; // Skip unused 
      Extension& extn = *colorExtensions[iext];
      if (extn.keepers.empty()) continue; // Not using any colors from this catalog
      string filename;
      extensionTable.readCell(filename, "Filename", iext);
      /**/cerr << "# Reading color catalog " << iext
	       << "/" << colorExtensions.size()
	       << " from " << filename << endl;
      int hduNumber;
      extensionTable.readCell(hduNumber, "HDUNumber", iext);
      string idKey;
      extensionTable.readCell(idKey, "idKey", iext);
      string colorExpression;
      extensionTable.readCell(colorExpression, "colorExpression", iext);
      stripWhite(colorExpression);
      if (colorExpression.empty()) {
	cerr << "No colorExpression specified for filename " << filename
	     << " HDU " << hduNumber
	     << endl;
	exit(1);
      }

      // Read the entire catalog for this extension
      FITS::FitsTable ft(filename, FITS::ReadOnly, hduNumber);
      FTable ff = ft.use();

      bool useRows = stringstuff::nocaseEqual(idKey, "_ROW");
      vector<long> id;
      if (useRows) {
	id.resize(ff.nrows());
	for (long i=0; i<id.size(); i++)
	  id[i] = i;
      } else {
	ff.readCells(id, idKey);
      }
      Assert(id.size() == ff.nrows());

      vector<double> color(id.size(), 0.);
      ff.evaluate(color, colorExpression);

      for (long irow = 0; irow < ff.nrows(); irow++) {
	map<long, Match*>::iterator pr = extn.keepers.find(id[irow]);
	if (pr == extn.keepers.end()) continue; // Not a desired object

	// Have a desired object. Put the color into everything it matches
	Match* m = pr->second;
	extn.keepers.erase(pr);

	for ( Match::iterator j = m->begin();
	      j != m->end();
	      ++j) 
	  (*j)->args.color = color[irow];
      } // End loop over catalog objects

      if (!extn.keepers.empty()) {
	cerr << "Did not find all desired objects in catalog " << filename
	     << " extension " << hduNumber
	     << endl;
	exit(1);
      }
    } // end loop over catalogs to read
	
    cerr << "Done reading catalogs for colors." << endl;

    // Make a pass through all matches to reserve as needed and purge 
    // those not wanted.  

    // ??? Count fitted detections per exposure to freeze parameters of those without enough??
    {
      ran::UniformDeviate u;
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
	  // Keep it if mag error is small or if it's from a tag or reference "instrument"
	  if ( d->sigmaMag > maxMagError
	       && exposures[ extensions[d->catalogNumber]->exposure]->instrument >= 0) {
	    j = (*im)->erase(j, true);  // will delete the detection too.
	  } else {
	    ++j;
	  }
	}
	int nFit = (*im)->countFit();
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
      cout << " chisq " << chi << " / " << dof << " dof maxdev " << sqrt(maxdev) << endl;

    } // End Match-culling section.

    ///////////////////////////////////////////////////////////
    // Construct priors
    ///////////////////////////////////////////////////////////

    // ??????????

    ///////////////////////////////////////////////////////////
    // Now do the re-fitting 
    ///////////////////////////////////////////////////////////

    // make CoordAlign class
    PhotoAlign ca(mapCollection, matches);

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
	     << " chisq " << chi << " / " << dof << " dof,  maxdev " << sqrt(maxdev) 
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
      long int dcount=0;
      long int mcount=0;
      
    } while (coarsePasses || nclip>0);
  
    // The re-fitting is now complete.  Serialize all the fitted magnitude solutions
    {
      ofstream ofs(outPhotFile.c_str());
      if (!ofs) {
	cerr << "Error trying to open output file for fitted photometry: " << outPhotFile << endl;
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
	     << " chisq " << chisq << " / " << dof << " dof,  maxdev " << sqrt(max) 
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
    FITS::FitsTable ft(outCatalog, FITS::ReadWrite + FITS::Create, "PhotoOut");
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
    vector<float> xPix;
    outTable.addColumn(xPix, "xPix");
    vector<float> yPix;
    outTable.addColumn(yPix, "yPix");
    vector<float> xExposure;
    outTable.addColumn(xExposure, "xExpo");
    vector<float> yExposure;
    outTable.addColumn(yExposure, "yExpo");
    vector<float> color;
    outTable.addColumn(color, "Color");

    vector<float> magOut;
    outTable.addColumn(magOut, "magOut");
    vector<float> sigMag;
    outTable.addColumn(sigMag, "sigMag");
    vector<float> magRes;
    outTable.addColumn(magRes, "magRes");

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
	outTable.writeCells(xPix, "xPix", pointCount);
	outTable.writeCells(yPix, "yPix", pointCount);
	outTable.writeCells(xExposure, "xExpo", pointCount);
	outTable.writeCells(yExposure, "yExpo", pointCount);
	outTable.writeCells(color, "Color", pointCount);
	outTable.writeCells(magOut, "magOut", pointCount);
	outTable.writeCells(magRes, "magRes", pointCount);
	outTable.writeCells(sigMag, "sigMag", pointCount);
	outTable.writeCells(wtFrac, "wtFrac", pointCount);

	pointCount += sequence.size();

	sequence.clear();
	catalogNumber.clear();
	objectNumber.clear();
	clip.clear();
	reserve.clear();
	xPix.clear();
	yPix.clear();
	xExposure.clear();
	yExposure.clear();
	color.clear();
	magOut.clear();
	magRes.clear();
	sigMag.clear();
	wtFrac.clear();
      }	// Done flushing the vectors to Table

      Match* m = *im;
      double mean;
      double wtot;
      m->getMean(mean, wtot);
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
	wtFrac.push_back( d->isClipped ? 0. : d->wt / wtot);

	xPix.push_back(d->args.xDevice);
	yPix.push_back(d->args.yDevice);
	xExposure.push_back(d->args.xExposure);
	yExposure.push_back(d->args.yExposure);
	color.push_back(d->args.color);

	magOut.push_back(d->magOut);
	// ?it's updated already, right??
	magRes.push_back( mean==0. ? 0. : d->magOut - mean);
	double sig=d->clipsq;
	sig = (sig > 0.) ? 1./sqrt(sig) : 0.;
	sigMag.push_back(sig);

	if (mean!=0.) {
	  // Accumulate statistics for meaningful residuals
	  if (dofPerPt >= 0. && !d->isClipped) {
	    int exposureNumber = extensions[d->catalogNumber]->exposure;
	    Exposure* expo = exposures[exposureNumber];
	    if (m->getReserved()) {
	      if (expo->instrument==REF_INSTRUMENT) {
		refAccReserve.add(d, mean, dofPerPt);
		vaccReserve[exposureNumber].add(d, mean, dofPerPt);
	      } else if (expo->instrument==TAG_INSTRUMENT) {
		// do nothing
	      } else {
		accReserve.add(d, mean, dofPerPt);
		vaccReserve[exposureNumber].add(d, mean, dofPerPt);
	      }
	    } else {
	      // Not a reserved object:
	      if (expo->instrument==REF_INSTRUMENT) {
		refAccFit.add(d, mean, dofPerPt);
		vaccFit[exposureNumber].add(d, mean, dofPerPt);
	      } else if (expo->instrument==TAG_INSTRUMENT) {
		// do nothing
	      } else {
		accFit.add(d, mean, dofPerPt);
		vaccFit[exposureNumber].add(d, mean, dofPerPt);
	      }
	    }
	  } // end statistics accumulation

	} // End residuals calculation

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
    outTable.writeCells(xPix, "xPix", pointCount);
    outTable.writeCells(yPix, "yPix", pointCount);
    outTable.writeCells(xExposure, "xExpo", pointCount);
    outTable.writeCells(yExposure, "yExpo", pointCount);
    outTable.writeCells(color, "Color", pointCount);
    outTable.writeCells(magOut, "magOut", pointCount);
    outTable.writeCells(magRes, "magRes", pointCount);
    outTable.writeCells(sigMag, "sigMag", pointCount);
    outTable.writeCells(wtFrac, "wtFrac", pointCount);

    // Print summary statistics for each exposure
    cout << "#    Exp    N    DOF    dmag    +-    RMS chi_red   xExposure    yExposure \n"
	 << "#                      |.....millimag....|          |.......degrees......|"
	 << endl;
    for (int iexp=0; iexp<exposures.size(); iexp++) {
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
    cout << "#                   N    DOF    dmag    +-    RMS chi_red \n"
	 << "#                             |.....millimag....|         "
	 << endl;
    cout << "Reference fit     " << refAccFit.summary() << endl;
    if (reserveFraction>0. && refAccReserve.n>0)
      cout << "Reference reserve " << refAccReserve.summary() << endl;

    cout << "Detection fit     " << accFit.summary() << endl;
    if (reserveFraction>0. && accReserve.n>0)
      cout << "Detection reserve " << accReserve.summary() << endl;

    // Cleanup:

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

/////////////////////////////////////////////////////////
PhotoMap* photoMapDecode(string code, string name, PolyMap::ArgumentType argType) {
  istringstream iss(code);
  string type;
  iss >> type;
  PhotoMap* pm;
  if (stringstuff::nocaseEqual(type,"Poly")) {
    int xOrder, yOrder;
    if (!(iss >> xOrder)) 
      throw runtime_error("Could not decode PolyMap spec: <" + code + ">");
    if ( (iss >> yOrder)) {
      PolyMap* poly = new PolyMap(xOrder, yOrder,argType, name);
      // Set PolyMap to identity for starters:
      DVector p = pm->getParams();
      p.setZero();
      pm->setParams(p);
      pm = poly;
    } else {
      PolyMap* poly = new PolyMap(xOrder, argType, name);
      // PolyMap has zero coefficients, identity on mags:
      DVector p = pm->getParams();
      p.setZero();
      pm->setParams(p);
      pm = poly;
    }
    /** Not in place yet:
  } else if (stringstuff::nocaseEqual(type,"Template")) {
    string filename;
    iss >> filename;
    ifstream ifs(filename.c_str());
    if (!ifs) 
      throw runtime_error("Could not open file for building TemplateMap1d <" + filename + ">");
    pm = TemplateMap1d::create(ifs, name);
    **/
  } else if (stringstuff::nocaseEqual(type,"Identity")) {
    pm = new IdentityMap;
  } else {
    throw runtime_error("Unknown type of PhotoMap: <" + type + ">");
  }
  return pm;
}
