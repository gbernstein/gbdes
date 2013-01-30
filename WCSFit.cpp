// Program to fit astrometric solutions to detection catalogs already matched by WCSFoF.
#include <fstream>
#include <sstream>
#include <map>

#include "Std.h"
#include "Astrometry.h"
#include "FitsTable.h"
#include "StringStuff.h"
#include "Pset.h"
#include "PixelMapCollection.h"
#include "TPVMap.h"
#include "Random.h"
#include "Match.h"

using namespace std;
using namespace astrometry;
using namespace stringstuff;
using img::FTable;

string usage=
  "WCSFit: Refine coordinate solutions for a matched set of catalogs.\n"
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

// Note that pixel maps for devices within instrument will get names <instrument>/<device>.
// And Wcs's for individual exposures will get names <exposure>/<device>.
// There will be SubMaps by these names too that include the reprojection to field coordinates.

#include "Instrument.h"

class MapMarq {
public:
  MapMarq(const list<Detection*>& dets_, PixelMap* pm_): dets(dets_), pm(pm_) {}
  void operator()(const DVector& a, double& chisq, 
		  DVector& beta, tmv::SymMatrix<double>& alpha);
private:
  const list<Detection*>& dets;
  PixelMap* pm;
};

class Accum {
public:
  double sumx;
  double sumy;
  double sumxx;
  double sumyy;
  double sumxxw;
  double sumyyw;
  double sumxw;
  double sumyw;
  double sumwx;
  double sumwy;
  int n;
  double xpix;
  double ypix;
  double xw;
  double yw;
  double sumdof;
  Accum(): sumxw(0.), sumyw(0.), 
	   sumx(0.), sumy(0.), sumwx(0.), sumwy(0.), 
	   sumxx(0.), sumyy(0.), sumxxw(0.), sumyyw(0.),
	   sumdof(0.), xpix(0.), ypix(0.), xw(0.), yw(0.),
	   n(0) {}
  void add(const Detection* d, double xoff=0., double yoff=0., double dof=1.) {
    sumx += (d->xw-xoff);
    sumy += (d->yw-yoff);
    sumxw += (d->xw-xoff)*d->wtx;
    sumyw += (d->yw-yoff)*d->wty;
    sumxxw += (d->xw-xoff)*(d->xw-xoff)*d->wtx;
    sumyyw += (d->yw-yoff)*(d->yw-yoff)*d->wty;
    sumxx += (d->xw-xoff)*(d->xw-xoff);
    sumyy += (d->yw-yoff)*(d->yw-yoff);
    sumwx += d->wtx;
    sumwy += d->wty;
    sumdof += dof;
    ++n;
  }
  double rms() const {
    return n > 0 ? sqrt( (sumxx+sumyy)/(2.*n)) : 0.;
  }
  double reducedChisq() const {
    return sumdof>0 ? (sumxxw+sumyyw)/(2.*sumdof) : 0.;
  }
  string summary() const {
    ostringstream oss;
    double dx = 0., dy = 0., sigx=0., sigy=0.;
    if (n>0) {
      dx = sumxw / sumwx;
      sigx = 1./sqrt(sumwx);
      dy = sumyw / sumwy;
      sigy = 1./sqrt(sumwy);
    }
    oss << setw(4) << n 
	<< fixed << setprecision(1)
	<< " " << setw(6) << sumdof
	<< " " << setw(5) << dx*1000.*DEGREE/ARCSEC 
	<< " " << setw(5) << sigx*1000.*DEGREE/ARCSEC
	<< " " << setw(5) << dy*1000.*DEGREE/ARCSEC 
	<< " " << setw(5) << sigy*1000.*DEGREE/ARCSEC
	<< " " << setw(5) << rms()*1000.*DEGREE/ARCSEC
	<< setprecision(2) 
	<< " " << setw(5) << reducedChisq()
	<< setprecision(0) << noshowpoint
	<< " " << setw(5) << xpix 
	<< " " << setw(5) << ypix
	<< setprecision(5) << showpoint << showpos
	<< " " << setw(9) << xw 
	<< " " << setw(9) << yw ;
    return oss.str();
  }
};

int
main(int argc, char *argv[])
{
  bool reFit;	// until we separate the statistics part ???
  double reserveFraction;
  double clipThresh;
  double maxPixError;
  double pixSysError;
  double referenceSysError;
  int minMatches;
  bool clipEntireMatch;
  int exposureOrder;
  int deviceOrder;
  string outCatalog;
  string catalogSuffix;
  string newHeadSuffix;
  double chisqTolerance;
  string renameInstruments;
  string existingMaps;
  string fixMaps;
  string canonicalExposures;
  string outWcs;
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
    parameters.addMemberNoValue("CLIPPING");
    parameters.addMember("clipThresh",&clipThresh, def | low,
			 "Clipping threshold (sigma)", 5., 2.);
    parameters.addMember("clipEntireMatch",&clipEntireMatch, def,
			 "Discard entire object if one outlier", false);
    parameters.addMemberNoValue("FITTING");
    parameters.addMember("reserveFraction",&reserveFraction, def | low,
			 "Fraction of matches reserved from re-fit", 0., 0.);
    parameters.addMember("exposureOrder",&exposureOrder, def | low,
			 "Order of per-exposure map", 1, 0);
    parameters.addMember("deviceOrder",&deviceOrder, def | low,
			 "Order of per-extension map", 2, 0);
    parameters.addMember("chisqTolerance",&chisqTolerance, def | lowopen,
			 "Fractional change in chisq for convergence", 0.001, 0.);
    parameters.addMember("renameInstruments",&renameInstruments, def,
			 "list of new names to give to instruments for maps","");
    parameters.addMember("existingMaps",&existingMaps, def,
			 "list of astrometric maps to draw from existing files","");
    parameters.addMember("fixMaps",&fixMaps, def,
			 "list of map components or instruments to hold fixed","");
    parameters.addMember("canonicalExposures",&canonicalExposures, def,
			 "list of exposures that will have identity exposure maps","");

    parameters.addMemberNoValue("FILES");
    parameters.addMember("outCatalog",&outCatalog, def,
			 "Output FITS binary catalog", "wcscat.fits");
    parameters.addMember("outWcs",&outWcs, def,
			 "Output serialized Wcs systems", "wcsfit.wcs");

    // ?? a grammar to define instrument map forms beyond polynomial?
    // ?? distinct outputs for each extension's Wcs and TPV fits?
  }

  // Positional accuracy (in degrees) demanded of numerical solutions for inversion of 
  // pixel maps: 
  const double worldTolerance = 0.001*ARCSEC/DEGREE;
  // Fractional reduction in RMS required to continue sigma-clipping:
  const double minimumImprovement=0.02;
  const int REF_INSTRUMENT=-1;	// Instrument for reference objects (no WCS fitting)
  const int TAG_INSTRUMENT=-2;	// Exposure number for tag objects (no WCS fitting nor contrib to stats)

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

    referenceSysError *= ARCSEC/DEGREE;

    // ??? log parameters to output file somehow?

    /////////////////////////////////////////////////////
    // Parse all the parameters describing maps etc. ????
    /////////////////////////////////////////////////////

    const char listSeperator=',';

    // First is a regex map from instrument names to the names of their PixelMaps
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

    // Next is regex map from PixelMap or instrument map names to files 
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

    // This is list of regexes of PixelMap names (or instrument names) that should
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
	 i != canonicalExposureList.end(); )
      if (i->empty()) {
	i = fixMapList.erase(i);
      } else {
	stripWhite(*i);
	++i;
      }

    /////////////////////////////////////////////////////
    //  Read in properties of all Fields, Instruments, Devices, Exposures
    /////////////////////////////////////////////////////

    // All we care about fields are names and orientations:
    NameIndex fieldNames;
    vector<SphericalCoords*> fieldProjections;

    {
      FITS::FitsTable in(inputTables, FITS::ReadOnly, "Fields");
      FITS::FitsTable out(outCatalog, FITS::ReadWrite + FITS::OverwriteFile + FITS::Create, "Fields");
      FTable ft = in.extract();
      out.adopt(ft);
      vector<double> ra;
      vector<double> dec;
      vector<string> name;
      ft.readCells(name, "Name");
      ft.readCells(ra, "RA");
      ft.readCells(dec, "Dec");
      for (int i=0; i<name.size(); i++) {
	fieldNames.append(name[i]);
	Orientation orient(SphericalICRS(ra[i]*DEGREE, dec[i]*DEGREE));
	fieldProjections.push_back( new Gnomonic(orient));
      }
    }

    // Let's figure out which of our FITS extensions are Instrument or MatchCatalog
    vector<int> instrumentHDUs;
    vector<int> catalogHDUs;
    {
      // Find which extensions are instrument tables
      FITS::FitsFile ff(inputTables);
      for (int i=1; i<ff.HDUCount(); i++) {
	FITS::Hdu h(inputTables, FITS::HDUAny, i);
	if (stringstuff::nocaseEqual(h.getName(), "Instrument"))
	  instrumentHDUs.push_back(i);
	else if (stringstuff::nocaseEqual(h.getName(), "MatchCatalog"))
	  catalogHDUs.push_back(i);
      }
    }
    
    // Read in all the instrument extensions and their device info.
    vector<Instrument*> instruments(instrumentHDUs.size());
    for (int i=0; i<instrumentHDUs.size(); i++) {
      FITS::FitsTable ft(inputTables, FITS::ReadOnly, instrumentHDUs[i]);
      FITS::FitsTable out(outCatalog, FITS::ReadWrite+FITS::Create, "Instrument");
      Assert( stringstuff::nocaseEqual(ft.getName(), "Instrument"));
      FTable ff=ft.extract();
      out.adopt(ff);
      string instrumentName;
      int instrumentNumber;
      if (!ff.header()->getValue("Name", instrumentName)
	  || !ff.header()->getValue("Number", instrumentNumber)) {
	cerr << "Could not read name and/or number of instrument at extension "
	     << instrumentHDUs[i] << endl;
      }
      stripWhite(instrumentName);
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
	instptr->addDevice( devnames[j],
			    Bounds<double>( vxmin[j], vxmax[j], vymin[j], vymax[j]));
      }
    }

    // Read in the table of exposures
    vector<Exposure*> exposures;
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
      for (int i=0; i<names.size(); i++) {
	// The projection we will use for this exposure:
	Gnomonic gn(Orientation(SphericalICRS(ra[i]*DEGREE,dec[i]*DEGREE)));
	Exposure* expo = new Exposure(names[i],gn);
	expo->field = fieldNumber[i];
	expo->instrument = instrumentNumber[i];
	exposures.push_back(expo);
      }
    }

    // Read and keep the table giving instructions for all input files:
    FTable fileTable;
    {
      FITS::FitsTable ft(inputTables, FITS::ReadOnly, "Files");
      FITS::FitsTable out(outCatalog, FITS::ReadWrite+FITS::Create, "Files");
      fileTable = ft.extract();
      out.copy(fileTable);
    }

    // Read info about all Extensions - we will keep the Table around.
    FTable extensionTable;
    vector<Extension*> extensions;
    {
      FITS::FitsTable ft(inputTables, FITS::ReadOnly, "Extensions");
      extensionTable = ft.extract();
      FITS::FitsTable out(outCatalog, FITS::ReadWrite+FITS::Create, "Extensions");
      out.copy(extensionTable);
    }
    for (int i=0; i<extensionTable.nrows(); i++) {
      extensions.push_back(new Extension);
      Extension* extn = extensions.back();
      int j;
      extensionTable.readCell(j, "ExposureNumber", i);
      extn->exposure = j;
      extensionTable.readCell(j, "DeviceNumber", i);
      extn->device = j;
      string s;
      extensionTable.readCell(s, "WCS", i);
      if (stringstuff::nocaseEqual(s, "ICRS")) {
	// Create a Wcs that just takes input as RA and Dec in degrees;
	IdentityMap identity;
	SphericalICRS icrs;
	extn->startWcs = new Wcs(&identity, icrs, "ICRS_degrees", DEGREE);
      } else {
	istringstream iss(s);
	PixelMapCollection pmcTemp;
	if (!pmcTemp.read(iss)) {
	  cerr << "Did not find Wcs format for ???" << endl;
	  exit(1);
	}
	string wcsName = pmcTemp.allWcsNames().front();
	extn->startWcs = pmcTemp.cloneWcs(wcsName);
      }
      // destination projection is field projection:
      extn->startWcs->reprojectTo(*fieldProjections[exposures[extn->exposure]->field]);

      string tpvOutFile;
      extensionTable.readCell(tpvOutFile, "TPVOut", i);
      extn->tpvOutFile = tpvOutFile;
    }

    /////////////////////////////////////////////////////
    //  Create and initialize all coordinate maps
    /////////////////////////////////////////////////////

    // Logic now for initializing instrument maps:
    // 1) Assign a name to each instrument map; either it's name of instrument,
    //    or we have a mapping from input parameters
    // 2) If the instrument name is something we were told to reuse:
    //        * get its form and values from the input map libraries for each device
    // 3) If not reusing, make new maps for each device.
    //        * Select canonical exposure either from cmd line or finding one
    //        * set reference exposure's map to identity
    //        * fit the initial instrument map to starting map.
    // 4) Fit initial exposure maps for each exposure that is not already canonical.
    // 5) Freeze parameters of all map elements specified in input.

    // Create a PixelMapCollection to hold all the components of the PixelMaps.
    PixelMapCollection mapCollection;
    // Create an identity map used by all reference catalogs
    mapCollection.learnMap(IdentityMap());

    // This is list of PixelMaps whose parameters will be held fixed during fitting:
    list<string> fixedMapNames;

    // Set up coordinate maps for each instrument & exposure
    for (int iinst=0; iinst<instruments.size(); iinst++) {
      Instrument* inst = instruments[iinst];

      // (1) choose a name to be used for this Instrument's PixelMap:
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

      PixelMapCollection pmcInstrument;
      if (loadExistingMaps) {
	ifstream ifs(loadFile.c_str());
	if (!ifs) {
	  cerr << "Could not open file " + loadFile + " holding existing PixelMaps" << endl;
	  exit(1);
	}
	pmcInstrument.read(ifs);
      }

      for (int idev=0; idev<inst->nDevices; idev++) {
	string devMapName = mapName + '/' + inst->deviceNames.nameOf(idev);
	inst->mapNames[idev] = devMapName;
	// A PixelMap file for Device overrides one for the Instrument
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
	  
	// Load the PixelMap for this device if it is to be extracted from a file.
	// We will not be throwing exceptions for duplicate names
	PixelMap* deviceMap;
	if (devLoadFile==loadFile) {
	  // Get map from the Instrument's collection:
	  deviceMap = pmcInstrument.issueMap(devMapName);
	  mapCollection.learnMap(*deviceMap);
	} else {
	  // Device has its own serialized file:
	  PixelMapCollection pmcDev;
	  ifstream ifs(devLoadFile.c_str());
	  if (!ifs) {
	    cerr << "Could not open file " + devLoadFile + " holding existing PixelMaps" << endl;
	    exit(1);
	  }
	  pmcDev.read(ifs);
	  deviceMap = pmcDev.issueMap(devMapName);
	  mapCollection.learnMap(*deviceMap);
	}

	if (fixDeviceMaps[idev])
	  fixedMapNames.push_back(devMapName);
      } // end of device loop

      // Canonical exposure will be used to initialize new Device maps.
      // And its Exposure map will be fixed at identity if none of the Device's maps are
      // fixed at initial values;  we assume that any single device being fixed will
      // break degeneracy between Exposure and Instrument models.

      // First find all exposures from this instrument,
      // and find if one is canonical;
      long canonicalExposure = -1;
      list<long> exposuresWithInstrument;
      for (long iexp = 0; iexp < exposures.size(); iexp++) {
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

	/**/cerr << "noDevicesFixed: " << noDevicesFixed
		 << " needCanonical: " << needCanonical
		 << " canonicalExposure: " << canonicalExposure
		 << endl;
	
      } // end finding a canonical exposure for the Instrument

      // Now we create new PixelMaps for each Device that does not already have one.
      for (int idev=0; idev < inst->nDevices; idev++) {
	if (deviceMapsExist[idev]) continue;

	// Create a new PixelMap for this device. ???? More elaborate choices ???
	PixelMap* pm = new PolyMap(deviceOrder, inst->mapNames[idev],
				   worldTolerance);

	// Find extension that has this Device for canonical Exposure:
	long iextn;
	for (iextn=0; iextn<extensions.size(); iextn++)
	  if (extensions[iextn]->exposure == canonicalExposure
	      && extensions[iextn]->device == idev) break;
	if (iextn >= extensions.size()) {
	  cerr << "Did not find extension for canonical exposure " 
	       << exposures[canonicalExposure]->name
	       << " and device " << inst->deviceNames.nameOf(idev)
	       << endl;
	  exit(1);
	}

	Exposure& expo=*exposures[canonicalExposure];
	Extension& extn = *extensions[iextn];

	// Initialize the new PixelMap so that it is best fit to the canonical exposures'
	// map from pixel coordinates to world coords=projection in exposure frame
	extn.startWcs->reprojectTo(*expo.projection);

	// ??? The following should be put into a subroutine and able to do nonlinear fit:
	{
	  const int nGridPoints=400;	// Number of test points for map initialization

	  Bounds<double> b=inst->domains[idev];
	  double step = sqrt(b.area()/nGridPoints);
	  int nx = static_cast<int> (ceil((b.getXMax()-b.getXMin())/step));
	  int ny = static_cast<int> (ceil((b.getYMax()-b.getYMin())/step));
	  double xstep = (b.getXMax()-b.getXMin())/nx;
	  double ystep = (b.getYMax()-b.getYMin())/ny;
	  list<Detection*> testPoints;

	  for (int ix=0; ix<=nx; ix++) {
	    for (int iy=0; iy<=ny; iy++) {
	      double xpix = b.getXMin() + ix*xstep;
	      double ypix = b.getYMin() + iy*ystep;
	      double xw, yw;
	      extn.startWcs->toWorld(xpix, ypix, xw, yw);
	      Detection* d = new Detection;
	      d->xpix = xpix;
	      d->ypix = ypix;
	      d->xw = xw;
	      d->yw = yw;
	      // Note MapMarq does not use positional uncertainties.
	      testPoints.push_back(d);
	    }
	  }
	  MapMarq mm(testPoints, pm);
	  // Since polynomial fit is linear, should just take one iteration:
	  double var=0.;
	  DVector beta(pm->nParams(), 0.);
	  DVector params(pm->nParams(), 0.);
	  tmv::SymMatrix<double> alpha(pm->nParams(), 0.);
	  mm(params, var, beta, alpha);
	  beta /= alpha;
	  params = beta;
	  beta.setZero();
	  alpha.setZero();
	  var = 0.;
	  mm(params, var, beta, alpha);
	  beta /= alpha;
	  params += beta;
	  pm->setParams(params);

	  // Clear out the testPoints:
	  for (list<Detection*>::iterator i = testPoints.begin();
	       i != testPoints.end();
	       ++i)
	    delete *i;
	} // Done initializing parameters of new PixelMap

	// Save this map into the collection
	mapCollection.learnMap(*pm);
      } // end loop creating PixelMaps for all Devices.

      // Now create an exposure map for all exposures with this instrument,
      // and set initial parameters by fitting to input WCS map.

      for (list<long>::const_iterator i= exposuresWithInstrument.begin();
	   i != exposuresWithInstrument.end();
	   ++i) {
	int iexp = *i;
	Exposure& expo = *exposures[iexp];
	PixelMap* warp=0;
	// First we check to see if we are going to use an existing exposure map:
	string loadFile = expo.name;
	if (existingMapFinder(loadFile)) {
	  // Adopt the existing map for this exposure:
	  PixelMapCollection pmcExpo;
	  ifstream ifs(loadFile.c_str());
	  if (!ifs) {
	    cerr << "Could not open serialized map file " 
		 << loadFile
		 << endl;
	    exit(1);
	  }
	  pmcExpo.read(ifs);
	  warp = pmcExpo.issueMap(expo.name);
	  mapCollection.learnMap(*warp);
	  expo.mapName = warp->getName();
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
	
	// We will create a new exposure map and then initialize it
	// ??? Potentially more sophisticated here
	if (exposureOrder==1) 
	  warp = new LinearMap(expo.name);
	else
	  warp = new PolyMap(exposureOrder, expo.name, worldTolerance);
	expo.mapName = expo.name;
	// Set numerical-derivative step to 1 arcsec, polynomial is working in degrees:
	warp->setPixelStep(ARCSEC/DEGREE);

	// Build test points
	list<Detection*> testPoints;

	const int nGridPoints=400;	// Number of test points for map initialization
	int pointsPerDevice = nGridPoints / inst->nDevices + 1;

	// Get points from each extension that is part of this exposure
	for (long iextn=0; iextn<extensions.size(); iextn++) {
	  if (extensions[iextn]->exposure != iexp) continue;
	  Extension& extn = *extensions[iextn];
	  int idev = extn.device;
	  Bounds<double> b=inst->domains[idev];
	  double step = sqrt(b.area()/pointsPerDevice);
	  int nx = static_cast<int> (ceil((b.getXMax()-b.getXMin())/step));
	  int ny = static_cast<int> (ceil((b.getYMax()-b.getYMin())/step));
	  double xstep = (b.getXMax()-b.getXMin())/nx;
	  double ystep = (b.getYMax()-b.getYMin())/ny;
	  // starting pixel map for the first exposure with this instrument
	  if (!extn.startWcs) {
	    cerr << "Failed to get starting Wcs for device "
		 << inst->deviceNames.nameOf(idev)
		 << " in exposure " << expo.name
		 << endl;
	    exit(1);
	  }
	  extn.startWcs->reprojectTo(*expo.projection);
	  PixelMap* devpm = mapCollection.issueMap(inst->mapNames[idev]);
	
	  for (int ix=0; ix<=nx; ix++) {
	    for (int iy=0; iy<=ny; iy++) {
	      double xpix = b.getXMin() + ix*xstep;
	      double ypix = b.getYMin() + iy*ystep;
	      double xw, yw;
	      // Map coordinates into the exposure projection = world:
	      extn.startWcs->toWorld(xpix, ypix, xw, yw);
	      // And process the pixel coordinates through Instrument map:
	      devpm->toWorld(xpix,ypix,xpix,ypix);
	      Detection* d = new Detection;
	      d->xpix = xpix;
	      d->ypix = ypix;
	      d->xw = xw;
	      d->yw = yw;
	      testPoints.push_back(d);
	    }
	  } // finish collecting test points for this extension
	} // finish extension loop for this exposure
	
	{
	  // Should be able to put this into a subroutine: ????
	  MapMarq mm(testPoints, warp);
	  // Since polynomial fit is linear, should just take one iteration:
	  double var=0.;
	  DVector beta(warp->nParams(), 0.);
	  DVector params(warp->nParams(), 0.);
	  tmv::SymMatrix<double> alpha(warp->nParams(), 0.);
	  mm(params, var, beta, alpha);
	  beta /= alpha;
	  params = beta;
	  beta.setZero();
	  alpha.setZero();
	  var = 0.;

	  mm(params, var, beta, alpha);
	  beta /= alpha;
	  params += beta;
	  warp->setParams(params);
	}
	
	
	// Insert the Exposure map into the collection
	mapCollection.learnMap(*warp);

	// Check for fixed maps
	if (regexMatchAny(fixMapList, warp->getName())) {
	  cerr << "WARNING: Requested fixed parameters for Exposure map "
	       << warp->getName()
	       << " that does not have initial values.  Ignoring request to fix params."
	       << endl;
	}

	delete warp;

	// Clear out the testPoints:
	for (list<Detection*>::iterator i = testPoints.begin();
	     i != testPoints.end();
	     ++i)
	  delete *i;
      } // end solution for this exposure map

    } // End instrument loop

    // Freeze the parameters of all the fixed maps
    mapCollection.setFixed(fixedMapNames, true);

    // Now construct Wcs and a reprojected SubMap for every extension

    // and save SubMap for it
    for (int iext=0; iext<extensions.size(); iext++) {
      Extension& extn = *extensions[iext];
      Exposure& expo = *exposures[extn.exposure];
      int ifield = expo.field;
      if ( expo.instrument < 0) {
	// Tag & reference exposures have no instruments and no fitting
	// being done.  Coordinates are fixed to xpix = xw.
	// Build a simple Wcs that takes its name from the field
	SphericalICRS icrs;
	mapCollection.defineWcs(fieldNames.nameOf(ifield), icrs, IdentityMap().getName(),
				DEGREE);
	extn.wcs = mapCollection.issueWcs(fieldNames.nameOf(ifield));
	// And have this Wcs reproject into field coordinates, then save as a SubMap
	// and re-issue for use in this extension.
	extn.wcs->reprojectTo(*fieldProjections[ifield]);
	mapCollection.learnMap(*extn.wcs);
	extn.map = mapCollection.issueMap(fieldNames.nameOf(ifield));
      } else {
	// Real instrument, make a map combining its exposure with its Device map:
	string wcsName = expo.name + "/" 
	  + instruments[expo.instrument]->deviceNames.nameOf(extn.device);
	list<string> mapElements;
	// The map itself is the device map:
	mapElements.push_back(instruments[expo.instrument]->mapNames[extn.device]);
	// Followed by the exposure map:
	mapElements.push_back(expo.mapName);
	string chainName = wcsName + "_chain";
	mapCollection.defineChain( chainName, mapElements);

	// Create a Wcs that projects into this exposure's projection:
	mapCollection.defineWcs(wcsName, *expo.projection, chainName, DEGREE);
	extn.wcs = mapCollection.issueWcs(wcsName);

	// Reproject this Wcs into the field system and get a SubMap including reprojection:
	extn.wcs->reprojectTo(*fieldProjections[ifield]);
	mapCollection.learnMap(*extn.wcs);
	extn.map = mapCollection.issueMap(wcsName);
      }
    } // Extension loop

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
      for (int i=0; i<seq.size(); i++) {
	if (seq[i]==0 && !extns.empty()) {
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
	extns.push_back(extn[i]);
	objs.push_back(obj[i]);
      } // End loop of catalog entries
      
      // Make a last match out of final entries
      if (!extns.empty()) {
	Detection* d = new Detection;
	d->catalogNumber = extns[0];
	d->objectNumber = objs[0];
	matches.push_back(new Match(d));
	extensions[extns[0]]->keepers.insert(std::pair<long,Detection*>(objs[0], d));
	for (int j=1; j<extns.size(); j++) {
	  Detection* d = new Detection;
	  d->catalogNumber = extns[j];
	  d->objectNumber = objs[j];
	  matches.back()->add(d);
	  extensions[extns[j]]->keepers.insert(std::pair<long, Detection*>(objs[j], d));
	}
      }
    } // End loop over input matched catalogs

    // Now loop over all original catalog bintables, reading the desired rows
    // and collecting needed information into the Detection structures
    for (int iext = 0; iext < extensions.size(); iext++) {
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
      neededColumns.push_back(errKey);

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
	if (errorColumnIsDouble) {
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
	  cout << "--Starting strict tolerance passes" << endl;
	  oldthresh = thresh;
	  nclip = ca.sigmaClip(thresh, false, clipEntireMatch);
	  cout << "Clipped " << nclip
	       << " matches " << endl;
	  continue;
	} else {
	  // Done!
	  break;
	}
      }
      oldthresh = thresh;
      nclip = ca.sigmaClip(thresh, false, clipEntireMatch);
      if (nclip==0 && coarsePasses) {
	// Nothing being clipped; tighten tolerances and re-fit
	coarsePasses = false;
	ca.setRelTolerance(chisqTolerance);
	cout << "No clipping--Starting strict tolerance passes" << endl;
	continue;
      }
      cout << "Clipped " << nclip
	   << " matches " << endl;
      long int dcount=0;
      long int mcount=0;
      
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

    // Save the new header files
    
    // Keep track of those we've already written to:
    set<string> alreadyOpened;

    // accuracy of SCAMP-format fit to real solution:
    const double SCAMPTolerance=0.0001*ARCSEC/DEGREE;
  
    for (int iext=0; iext<extensions.size(); iext++) {
      Extension& extn = *extensions[iext];
      Exposure& expo = *exposures[extn.exposure];
      if (expo.instrument < 0) continue;
      Instrument& inst = *instruments[expo.instrument];
      // Open file for new ASCII headers  ??? Save serialized per exposure too??
      string newHeadName = extn.tpvOutFile;
      bool overwrite = false;
      if ( alreadyOpened.find(newHeadName) == alreadyOpened.end()) {
	overwrite = true;
	alreadyOpened.insert(newHeadName);
      }
      ofstream newHeads(newHeadName.c_str(), overwrite ? 
			ios::out : (ios::out | ios::in | ios::ate));
      if (!newHeads) {
	cerr << "WARNING: could not open new header file " << newHeadName << endl;
	cerr << "...Continuing with program" << endl;
	continue;
      }
      // Fit the TPV system to a projection about the field's pointing center:
      expo.projection->setLonLat(0.,0.);
      Wcs* tpv = fitTPV(inst.domains[extn.device],
			*extn.wcs,
			*expo.projection,
			"NoName",
			SCAMPTolerance);
      newHeads << writeTPV(*tpv);
      delete tpv;
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
	sigw.push_back(sigmaWorld);

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
    cout << "#                   N    DOF    dx    +-    dy    +-   RMS chi_red  \n"
	 << "#                             |.....milliarcsec...........|         "
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
      delete extensions[i];
    // Get rid of exposures
    for (int i=0; i<exposures.size(); i++)
      delete exposures[i];
    // Get rid of instruments
    for (int i=0; i<instruments.size(); i++)
      delete instruments[i];

  } catch (std::runtime_error& m) {
    quit(m,1);
  }
}
    
void
MapMarq::operator()(const DVector& a, 
		    double& chisq, 
		    DVector& beta, 
		    tmv::SymMatrix<double>& alpha) {
  chisq = 0;
  beta.setZero();
  alpha.setZero();
  int np = a.size();
  Assert(pm->nParams()==np);
  Assert(beta.size()==np);
  Assert(alpha.nrows()==np);
  DMatrix derivs(2,np);
  pm->setParams(a);

  for (list<Detection*>::const_iterator i = dets.begin();
       i != dets.end();
       ++i) {
    const Detection* d = *i;
    if (d->isClipped) continue;
    double xmod, ymod;
    pm->toWorldDerivs(d->xpix, d->ypix, xmod, ymod, derivs);
    xmod = d->xw - xmod;
    ymod = d->yw - ymod;

    chisq += xmod*xmod + ymod*ymod;
    for (int i=0; i<np; i++) {
      beta[i] += xmod*derivs(0,i) + ymod*derivs(1,i);
      for (int j=0; j<=i; j++) 
	alpha(i,j)+=derivs(0,i)*derivs(0,j) + derivs(1,i)*derivs(1,j);
    }
  }
}
