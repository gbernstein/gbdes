// Program to use photometric and astrometric solutions to 
// produce catalogs of selected magnitudes and colors for objects matched by WCSFoF.
// An output FOF catalog is produced that augments the input FOF catalog with appropriate
// tags to the magnitude and color catalogs that have been created.
#include <fstream>
#include <sstream>
#include <map>

#include "Std.h"
#include "Astrometry.h"
#include "FitsTable.h"
#include "StringStuff.h"
#include "Pset.h"
#include "PixelMapCollection.h"
#include "PhotoMatch.h"
#include "PhotoMapCollection.h"
#include "PhotoInstrument.h"
#include "TemplateMap.h"

using namespace std;
using namespace stringstuff;
using namespace photometry;
using img::FTable;
using FITS::FitsTable;

string usage=
  "MakeColors: Merge magnitudes for matched detections to create new\n"
  "            catalogs giving magnitudes and colors for each match.\n"
  "usage: MakeColors <match file> <mag file out> [parameter file] [parameter file...]\n"
  "      <match file>:  FITS file with binary tables produced by WCSFoF, will be amended\n"
  "      <mag file out>: Name of output FITS catalog file to hold magnitudes\n"
  "     stdin is read as parameters.";

// Temporary documentation:
// parameter renameInstruments is a string that has format like
// <instrument name>=<map name>, <instrument name> = <map name>, ...
// where the Instrument's PixelMap will be given the map name.  Regex allowed.
// Note that this is assuming that regexes do not include = or , characters.
// whitespace will be stripped from edges of each name.

// parameter magnitudes is a string with format
// <band name>@<instrument name>[, <band name>@<instrument name>] ...
// where <band name> is the name of a magnitude band, and <instrument name> is an
// instrument assumed to create valid mags for this band.  One band can come from
// multiple instruments.  Output magnitude catalog will have one column for each distinct band.
// Regex is ok for instrument specification.

// parameter magnitudeFile gives file name for the output FITS magnitude table.

// parameter colors is a string saying what colors will be calculated.  Format is
// <band name> - <band name> [+ offset] @ <file name> [, <band name> - ....]
// Each comma-separated specifier gives the two bands that will be differenced to produce the
// color that is written into a FITS table with given file name.

// Routine to return clipped mean & uncertainty given data and weights w:
void clipMeanAndError(vector<double>& x, vector<double>& w, double& mean, double& err,
		      double threshold);

// A helper function that strips white space from front/back of a string and replaces
// internal white space with underscores:
void
spaceReplace(string& s) {
  stripWhite(s);
  s = regexReplace("[[:space:]]+","_",s);
}

// Here is the default character at which to split lists given in parameters strings
const char DefaultListSeperator=',';

// Collect the regular expressions that will match instruments that are
// assign to each magnitude band
struct Band {
  Band(): number(-1), hasData(false) {}
  string name;
  int number;
  list<string> regexes;
  bool hasData;
};

class Color {
public:
  Color(): fits(0), data(0), rowCount(0), extensionNumber(-1) {}
  ~Color() {
    if (data) delete data;
    if (fits) delete fits;
  }
  string band1;
  string band2;
  string colorSpec;
  string filename;
  double offset;
  int bandNumber1;
  int bandNumber2;
  FitsTable* fits;	// Pointer to the FitsTable that will hold this color's catalog
  FTable* data;
  long rowCount;	// Number of objects in this color's table so far
  int extensionNumber;	// extension number assigned to this color's catalog;
};

// Collect this information for each desired detection in the input catalogs
class MagPoint {
public:
  double ra;
  double dec;
  double mag;
  double magErr;
};

// Another helper function to split up a string into a list of whitespace-trimmed strings.
// Get rid of any null ones.

list<string> splitArgument(string input, const char listSeperator = DefaultListSeperator) {
  list<string> out = split(input, listSeperator);
  for (list<string>::iterator i = out.begin();
       i != out.end(); )  {
    stripWhite(*i);
    if (i->empty()) {
      i = out.erase(i);
    } else {
      ++i;
    }
  }
  return out;
}


int
main(int argc, char *argv[])
{
  double clipThresh;
  double maxMagError;
  double magSysError;
  string renameInstruments;
  string magnitudes;
  string colors;
  string wcsFiles;
  string photoFiles;
  Pset parameters;
  {
    const int def=PsetMember::hasDefault;
    const int low=PsetMember::hasLowerBound;
    const int up=PsetMember::hasUpperBound;
    const int lowopen = low | PsetMember::openLowerBound;
    const int upopen = up | PsetMember::openUpperBound;

    parameters.addMember("maxMagError",&maxMagError, def | lowopen,
			 "Only output mags with errors below this level", 0.05, 0.);
    parameters.addMember("magSysError",&magSysError, def | lowopen,
			 "Systematic added in quadrature to mag error", 0.001, 0.);
    parameters.addMember("clipThresh",&clipThresh, def | low,
			 "Magnitude clipping threshold (sigma)", 5., 2.);
    parameters.addMember("magnitudes",&magnitudes, def,
			 "bands to output and instruments they come from","");
    parameters.addMember("colors",&colors, def,
			 "colors to tabulate and files they go into","");
    parameters.addMember("renameInstruments",&renameInstruments, def,
			 "list of new names to give to instruments for maps","");
    parameters.addMember("wcsFiles",&wcsFiles, def,
			 "files holding WCS maps to override starting WCSs","");
    parameters.addMember("photoFiles",&photoFiles, def,
			 "files holding single-band photometric solutions for input catalogs","");
}

  // Fractional reduction in RMS required to continue sigma-clipping:
  const double minimumImprovement=0.02;
  const int REF_INSTRUMENT=-1;	// Instrument for reference objects (no fitting) ??? color terms???
  const int TAG_INSTRUMENT=-2;	// Exposure number for tag objects (no fitting nor contrib to stats)

  const string stellarAffinity="STELLAR";

  const double NO_MAG_DATA = -100.;	// Value entered when there is no valid mag or color
  const string colorColumnName = "COLOR";
  const string colorErrorColumnName = "COLOR_ERR";

  try {
    
    // Read parameters
    if (argc<3) {
      cerr << usage << endl;
      cerr << "--------- Default Parameters: ---------" << endl;
      parameters.setDefault();
      parameters.dump(cerr);
      exit(1);
    }
    string inputTables = argv[1];
    string magOutFile = argv[2];

    parameters.setDefault();
    for (int i=3; i<argc; i++) {
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

    // Teach PixelMapCollection about new kinds of PixelMaps:
    astrometry::PixelMapCollection::registerMapType<astrometry::TemplateMap1d>();

    // Teach PhotoMapCollection about new kinds of PhotoMaps:
    // Currently all map types are automatically loaded in PhotoMapCollection constructor,
    // but if new ones are invented, this is the format:
    // PhotoMapCollection::registerMapType<PolyMap>();

    // First is a regex map from instrument names to the names of their PhotoMaps
    RegexReplacements instrumentTranslator;
    {
      list<string> ls = split(renameInstruments, DefaultListSeperator);
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


    // Get the list of files holding astrometric maps and read them.
    list<string> wcsFileList = splitArgument(wcsFiles);
    astrometry::PixelMapCollection astromaps;
    for (list<string>::const_iterator i = wcsFileList.begin();
	 i != wcsFileList.end();
	 ++i) {
      ifstream ifs(i->c_str());
      if (!astromaps.read(ifs))
	cerr << "WARNING: could not read astrometric map file " << *i << endl;
    }

    // Read the photometric solutions now:
    list<string> photoFileList = splitArgument(photoFiles);
    PhotoMapCollection photomaps;
    for (list<string>::const_iterator i = photoFileList.begin();
	 i != photoFileList.end();
	 ++i) {
      ifstream ifs(i->c_str());
      if (!photomaps.read(ifs))
	cerr << "WARNING: could not read astrometric map file " << *i << endl;
    }


    typedef map<string, Band> BandMap;
    BandMap bands;
    {
      list<string> bandspecs = splitArgument(magnitudes);
      for (list<string>::iterator i = bandspecs.begin();
	   i != bandspecs.end();
	   ++i) {
	list<string> bandInst = split(*i, '@');
	if (bandInst.size() != 2) {
	  cerr << "Bad band/instrument pair: " << *i << endl;
	  exit(1);
	}
	string bandName = bandInst.front();
	stripWhite(bandName);
	string instrumentRegex = bandInst.back();
	stripWhite(instrumentRegex);
	if (bandName.empty() || instrumentRegex.empty()) {
	  cerr << "Missing band or instrument: " << *i << endl;
	  exit(1);
	}
	bands[bandName].regexes.push_back(instrumentRegex);
      }
    }

    // Once all the bands have been read, assign them numbers
    {
      int bandCounter = 0;
      for (BandMap::iterator i = bands.begin();
	   i != bands.end();
	   ++i) {
	Band& b = i->second;
	b.name = i->first;
	b.number = bandCounter++;
      }
    }

    // Parse specifications for color.  Make a little struct to describe each color to be
    // extracted.


    list<Color> colorList;
    {
      list<string> colorspecs = splitArgument(colors);
      for (list<string>::iterator i = colorspecs.begin();
	   i != colorspecs.end();
	   ++i) {
	bool success = true;
	Color c;
	list<string> colorfile = split(*i, '@');
	success = (colorfile.size() ==2);
	if (success) {
	  // 2nd entry is filename
	  c.filename = colorfile.back();
	  stripWhite(c.filename);
	  // 1st entry is the specification of the desired color
	  c.colorSpec = colorfile.front();
	  stripWhite(c.colorSpec);

	  // Get possible offset constant from color description at '+'
	  list<string> coloroffset = split(colorfile.front(), '+');
	  success =  (coloroffset.size()>0 && coloroffset.size() <=2);
	  c.offset = 0.;
	  if (success && coloroffset.size() == 2) {
	    // get offset
	    istringstream iss(coloroffset.back());
	    success = (iss >> c.offset);
	  }
	  // split up the color term at -
	  if (success) {
	    list<string> twobands = split(coloroffset.front(), '-');
	    if (twobands.size() == 2) {
	      c.band1 = twobands.front();
	      c.band2 = twobands.back();
	      stripWhite(c.band1);
	      stripWhite(c.band2);
	    } else {
	      success = false;
	    }
	  }
	}
	if (!success) {
	  cerr << "Error in color specification: " << *i << endl;
	  exit(1);
	}

	// Check that both mags in the color are known to us:
	if (!bands.count(c.band1)) {
	  cerr << "Undefined magnitude <" << c.band1 << " requested for color" << endl;
	  exit(1);
	}
	if (!bands.count(c.band2)) {
	  cerr << "Undefined magnitude <" << c.band2 << " requested for color" << endl;
	  exit(1);
	}
	c.bandNumber1 = bands[c.band1].number;
	c.bandNumber2 = bands[c.band2].number;

	colorList.push_back(c);
      } // end color spec loop
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
    
    /**/cerr << "Found " << instrumentHDUs.size() << " instrument HDUs" 
	     << " and " << catalogHDUs.size() << " catalog HDUs" << endl;

    // This flag will be set if we have already opened (and overwritten) the
    // output FITS catalog.
    bool outputCatalogAlreadyOpen = false;

    // Read in all the instrument extensions and their device info.
    vector<Instrument*> instruments(instrumentHDUs.size(),0);
    // This array will say what band number any instrument gives data for.  -1 means none of interest
    vector<int> instrumentBandNumbers(instrumentHDUs.size(), -1);

    for (int i=0; i<instrumentHDUs.size(); i++) {
      FITS::FitsTable ft(inputTables, FITS::ReadOnly, instrumentHDUs[i]);

      string instrumentName;
      int instrumentNumber;
      if (!ft.header()->getValue("Name", instrumentName)
	  || !ft.header()->getValue("Number", instrumentNumber)) {
	cerr << "Could not read name and/or number of instrument at extension "
	     << instrumentHDUs[i] << endl;
      }
      spaceReplace(instrumentName);
      // Is this an instrument we are going to care about?
      string name = instrumentName;
      instrumentTranslator(name);  // Apply any name remapping specified.
      
      // Now see if this PhotoMap instrument name appears
      for (BandMap::const_iterator iBand = bands.begin();
	   iBand != bands.end();
	   ++iBand) {

	if (regexMatchAny(iBand->second.regexes, name))  {
	  // This is an instrument we will use. 
	  // Assign a band number
	  instrumentBandNumbers[i] = iBand->second.number;
	  // Get instrument's devices:
	  FTable ff=ft.extract();
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
	  break;	// Break out of band loop if we find one match.
	} // Done with section when finding a band
      } // Done with band loop
    } // Done with instrument loop

    // Read in the table of exposures
    vector<Exposure*> exposures;
    {
      FITS::FitsTable ft(inputTables, FITS::ReadOnly, "Exposures");
      FTable ff = ft.use();
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
	// Continue in loop if this exposure contains no useful mag info.
	if (instrumentBandNumbers[instrumentNumber[i]] < 0)
	  continue;

	spaceReplace(names[i]);
	// The projection we will use for this exposure:
	astrometry::Gnomonic gn(astrometry::Orientation(astrometry::SphericalICRS(ra[i]*DEGREE,
										  dec[i]*DEGREE)));
	Exposure* expo = new Exposure(names[i],gn);
	expo->field = fieldNumber[i];
	expo->instrument = instrumentNumber[i];
	exposures.push_back(expo);
      }

      // ???? Add the magnitude and color catalogs to the Exposure table ???
    }


    // Read info about all Extensions - we will keep the Table around.
    FTable extensionTable;
    vector<Extension*> extensions;
    FITS::FitsTable ft(inputTables, FITS::ReadOnly, "Extensions");
    extensionTable = ft.use();
    // This array will give band associated with each extension.  -1 for no desired band.
    vector<int> extensionBandNumbers(extensionTable.nrows(), -1);

    // ??? Add the mag & color catalogs to the Extension table ???
    // ??? and assign them extension numbers
    int magCatalogExtensionNumber = -1; // ????

    for (int i=0; i<extensionTable.nrows(); i++) {
      Extension* extn = new Extension;
      int iExposure;
      extensionTable.readCell(iExposure, "ExposureNumber", i);

      if (!exposures[iExposure]) {
	// This extension is not in an exposure of interest.  Skip it.
	delete extn;
	extn = 0;
	extensions.push_back(extn);
	continue;
      }
      extensions.push_back(extn);
      extn->exposure = iExposure;
      int iDevice;
      extensionTable.readCell(iDevice, "DeviceNumber", i);
      extn->device = iDevice;

      Exposure& expo = *exposures[iExposure];
      // Save away the band number for this extension
      extensionBandNumbers[i] = instrumentBandNumbers[expo.instrument];
      // and get the name its maps should be found under.
      string mapName = expo.name + "/" 
	+ instruments[expo.instrument]->deviceNames.nameOf(extn->device);

      string s;
      extensionTable.readCell(s, "WCS", i);
      if (stringstuff::nocaseEqual(s, "ICRS")) {
	// Create a Wcs that just takes input as RA and Dec in degrees;
	astrometry::IdentityMap identity;
	astrometry::SphericalICRS icrs;
	extn->startWcs = new astrometry::Wcs(&identity, icrs, "ICRS_degrees", DEGREE);
      } else {
	try {
	  // See if there is a map to use from the astrometric solution files
	  extn->startWcs = astromaps.issueWcs(mapName);
	} catch (PhotometryError) {
	  // else we read from the input file's starting WCS
	  istringstream iss(s);
	  astrometry::PixelMapCollection pmcTemp;
	  if (!pmcTemp.read(iss)) {
	    cerr << "Did not find WCS format for " << mapName << endl;
	    exit(1);
	  }
	  string wcsName = pmcTemp.allWcsNames().front();
	  extn->startWcs = pmcTemp.cloneWcs(wcsName);
	}
      }
      // destination projection is the Exposure projection, whose coords used for any
      // exposure-level magnitude corrections
      extn->startWcs->reprojectTo(*exposures[extn->exposure]->projection);

      // Get the photometric solution too
      extn->map = photomaps.issueMap(mapName);
      if (extn->map->needsColor()) {
	cerr << "MakeColors is not able to use the PhotoMap <" << mapName
	     << "> because it requires color information." << endl;
	exit(1);
      }
    }

    //////////////////////////////////////////////////////////
    // Read in all the data
    //////////////////////////////////////////////////////////

    // Start by reading all matched catalogs, noting which 
    // telling each Extension which objects it should retrieve from its catalog

    // Keep a set of desired object IDs from each extension
    vector<set<long> > desiredObjects(extensions.size());

    for (int icat = 0; icat < catalogHDUs.size(); icat++) {
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
	stripWhite(affinity);
	if (!stringstuff::nocaseEqual(affinity, stellarAffinity))
	  continue;
      }
      vector<LONGLONG> extn;
      vector<LONGLONG> obj;
      ff.readCells(extn, "Extension");
      ff.readCells(obj, "Object");

      for (int j = 0; j<extn.size(); j++) {
	if (!extensions[extn[j]]) continue;
	desiredObjects[extn[j]].insert(obj[j]);
      }
    } // End loop over input matched catalogs

    // One catalog for each extension, with object ID's as keys
    vector<map<long, MagPoint> > pointMaps(extensions.size());

    // Now loop over all original catalog bintables, reading the desired rows
    // and collecting needed information into the Detection structures
    for (int iext = 0; iext < extensions.size(); iext++) {
      if (!extensions[iext]) continue; // Skip unused 
      // Relevant structures for this extension
      Extension& extn = *extensions[iext];
      Exposure& expo = *exposures[extn.exposure];

      string filename;
      extensionTable.readCell(filename, "Filename", iext);
      /**/cerr << "# Reading object catalog " << iext
	       << "/" << extensions.size()
	       << " from " << filename 
	       << endl;
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

      astrometry::Wcs* startWcs = extn.startWcs;
      photometry::SubMap* photomap = extn.map;

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
      bool magColumnIsDouble = true;
      try {
	double d;
	ff.readCell(d, magKey, 0);
      } catch (img::FTableError& e) {
	magColumnIsDouble = false;
      }

      for (long irow = 0; irow < ff.nrows(); irow++) {
	set<long>::iterator pr = desiredObjects[iext].find(id[irow]);
	if (pr == desiredObjects[iext].end()) continue; // Not a desired object
	// Remove this object's ID from the list of those being sought.
	desiredObjects[iext].erase(pr);

	// Read the PhotometryArguments for this Detection:
	PhotoArguments args;
	double magIn;
	double magErr;

	MagPoint mp;

	ff.readCell(args.xDevice, xKey, irow);
	ff.readCell(args.yDevice, yKey, irow);
	astrometry::SphericalICRS sky = startWcs->toSky(args.xDevice, args.yDevice);
	sky.getLonLat(mp.ra, mp.dec);
	startWcs->toWorld(args.xDevice, args.yDevice,
			  args.xExposure, args.yExposure);
	args.color = 0;	// Note no color information is available!

	// Get the mag input and its error
	if (magColumnIsDouble) {
	  ff.readCell(magIn, magKey, irow);
	} else {
	  float f;
	  ff.readCell(f, magKey, irow);
	  magIn = f;
	}
	if (errorColumnIsDouble) {
	  ff.readCell(magErr, magErrKey, irow);
	} else {
	  float f;
	  ff.readCell(f, magErrKey, irow);
	  magErr = f;
	}

	// Map to output and estimate output error
	mp.mag = photomap->forward(magIn, args);
	mp.magErr = magErr * pow(photomap->derivative(magIn, args), 2.);

	// Add this point's info to our MagPoint catalog
	pointMaps[iext][id[irow]] = mp;
      } // End loop over input rows of this extension

      // Check that we have found all the objects we wanted:
      if (!desiredObjects[iext].empty()) {
	cerr << "Did not find all desired objects in catalog " << filename
	     << " extension " << hduNumber
	     << endl;
	exit(1);
      }
    } // end loop over extensions to read


    // Make output catalogs
    FitsTable magFitsTable(magOutFile, FITS::OverwriteFile);
    FTable magTable = magFitsTable.use();
    magTable.clear();	// should already be empty though
    {
      vector<double> dummy;
      magTable.addColumn(dummy, "RA");
      magTable.addColumn(dummy, "Dec");
      for (BandMap::const_iterator iBand = bands.begin();
	   iBand != bands.end();
	   ++iBand)
	magTable.addColumn(dummy, iBand->second.name);
	magTable.addColumn(dummy, iBand->second.name+"_ERR");
    }
    long magRowCounter = 0;	// Count of objects with assigned magnitudes

    // Make a catalog for each color being calculated
    for (list<Color>::iterator i = colorList.begin();
	 i != colorList.end();
	 ++i) {
      i->fits = new FitsTable(i->filename, FITS::OverwriteFile);
      i->data = &(i->fits->use());
      i->data->clear();	// Should be redundant
      i->data->header()->append("COLORSPEC",i->colorSpec);
      vector<double> dummy;
      i->data->addColumn(dummy, "RA");
      i->data->addColumn(dummy, "Dec");
      i->data->addColumn(dummy, colorColumnName);
      i->data->addColumn(dummy, colorErrorColumnName);
    }


    // Iterate over the match catalogs again, this time calculating mags & colors
    for (int icat = 0; icat < catalogHDUs.size(); icat++) {
      FITS::FitsTable ft(inputTables, FITS::ReadOnly, catalogHDUs[icat]);
      FTable ff = ft.extract();
      {
	// Only calculate mags and colors for stellar objects
	string affinity;
	if (!ff.getHdrValue("Affinity", affinity)) {
	  cerr << "Could not find affinity keyword in header of extension " 
	       << catalogHDUs[icat] 
	       << endl;
	  exit(1);
	}
	stripWhite(affinity);
	if (!stringstuff::nocaseEqual(affinity, stellarAffinity))
	  continue;
      }
      vector<int> seqIn;
      vector<LONGLONG> extnIn;
      vector<LONGLONG> objIn;
      ff.readCells(seqIn, "SequenceNumber");
      ff.readCells(extnIn, "Extension");
      ff.readCells(objIn, "Object");

      // We'll collect output data in these arrays, with mag & color entries 
      vector<int> seqOut;
      vector<LONGLONG> extnOut;
      vector<LONGLONG> objOut;
      seqOut.reserve(seqIn.size());
      extnOut.reserve(extnIn.size());
      objOut.reserve(objIn.size());

      // Arrays in which we'll average the positions and magnitudes:
      double sumRA = 0.;
      double sumDec = 0;
      int nRADec = 0;
      vector<vector<double> > bandWeights(bands.size());
      vector<vector<double> > bandMags(bands.size());
      
      int lastSeq = 0;
      for (int j = 0; j<=extnIn.size(); j++) {
	if (j==extnIn.size() || seqIn[j]==0) {
	  // Calculate and save away all mags and colors for previous match
	  bool hasMagnitudes = false;
	  DVector mags(bands.size());
	  DVector magErrors(bands.size());
	  vector<bool> magValid(bands.size(), false);
	  for (int i=0; i<mags.size(); i++) {
	    clipMeanAndError(bandMags[i], bandWeights[i], mags[i], magErrors[i],
			     clipThresh);
	    magValid[i] = magErrors[i]>0. && magErrors[i] <= maxMagError;
	    if (magValid[i]) 
	      hasMagnitudes = true;
	  }
	  if (hasMagnitudes) {
	    // Get RA & Dec in degrees
	    sumRA /= (nRADec * DEGREE);
	    sumDec /= (nRADec * DEGREE);
	    ff.writeCell(sumRA, "RA", magRowCounter);
	    ff.writeCell(sumDec, "Dec", magRowCounter);
	    // write mags & colors if there are any of use here
	    for (BandMap::const_iterator i= bands.begin();
		 i != bands.end();
		 ++i) {
	      int index = i->second.number;
	      if (magValid[index]) {
		ff.writeCell( mags[index],
			      i->second.name, magRowCounter);
		ff.writeCell( magErrors[index],
			      i->second.name+"_ERR", magRowCounter);
	      } else {
		// Mags with too-large errors get no-data marker
		ff.writeCell( NO_MAG_DATA,
			      i->second.name, magRowCounter);
		ff.writeCell( NO_MAG_DATA,
			      i->second.name+"_ERR", magRowCounter);
	      }
	    } // end mag band loop

	    //  add a seq for the mag catalog entry
	    seqOut.push_back(++lastSeq);
	    extnOut.push_back(magCatalogExtensionNumber);
	    objOut.push_back(magRowCounter++);
	    
	    //       for each possible color
	    for (list<Color>::iterator icolor = colorList.begin();
		 icolor != colorList.end();
		 ++icolor) {
	      //         do we have it?
	      if ( magValid[icolor->bandNumber1] && magValid[icolor->bandNumber2]) {
		//         write to color catalog
		icolor->data->writeCell(mags[icolor->bandNumber1] - mags[icolor->bandNumber2]
					+ icolor->offset,
					colorColumnName, 
					icolor->rowCount);
		icolor->data->writeCell(hypot(magErrors[icolor->bandNumber1], 
					      magErrors[icolor->bandNumber2]),
					colorErrorColumnName, 
					icolor->rowCount);
		icolor->data->writeCell(sumRA, "RA", icolor->rowCount);
		icolor->data->writeCell(sumDec, "Dec", icolor->rowCount);

		//   add seq entry for color catalog
		seqOut.push_back(++lastSeq);
		extnOut.push_back(icolor->extensionNumber);
		objOut.push_back(icolor->rowCount++);
	      } // End processing for a valid color
	    } // end color loop
	  } // end hasMagnitudes

	  //  clear mag & position accumulators
	  sumRA = sumDec = 0.;
	  nRADec = 0;
	  for (int i = 0; i<bandWeights.size(); i++) bandWeights[i].clear();
	  for (int i = 0; i<bandMags.size(); i++) bandMags[i].clear();
	} // End block for completing a matched set

	if (j==extnIn.size()) break;  // Quit if we don't have new data

	// Propagate every detection to output table
	seqOut.push_back(seqIn[j]);
	extnOut.push_back(extnIn[j]);
	objOut.push_back(objIn[j]);

	lastSeq = seqIn[j];

	// See if this detection has a mag we're interested in
	int iBand = extensionBandNumbers[extnIn[j]];
	if (iBand>=0 ) {
	  // Yes; accumulate its information into averages
	  MagPoint mp = pointMaps[extnIn[j]][objIn[j]]; // ?? check that it's here ???
	  sumRA += mp.ra;
	  sumDec += mp.dec;
	  nRADec++;
	  double weight = 1. / (mp.magErr*mp.magErr + magSysError*magSysError);
	  if (mp.magErr==0.) weight = 0.;
	  bandWeights[iBand].push_back(weight);
	  bandMags[iBand].push_back(mp.mag);
	}
      } // End object loop for this match catalog.

      // Write back the new sequences
      ff.writeCells(seqOut, "SequenceNumber", 0);
      ff.writeCells(extnOut, "Extension", 0);
      ff.writeCells(objOut, "Object", 0);

    } // End loop over input matched catalogs



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
