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
#include "Instrument.h"
#include "TemplateMap.h"
#include "FitsImage.h"

#include "Units.h"  // Specifies angular units used
#include "FitSubroutines.h"

using namespace std;
using namespace stringstuff;
using namespace photometry;
using img::FTable;
using FITS::FitsTable;
using photometry::DVector;
using astrometry::WCS_UNIT;

string usage=
  "MagColor: Merge magnitudes for matched detections to create new\n"
  "          catalogs giving magnitudes and colors for each match.\n"
  "usage: MagColor <match file in> <match file out> <mag file out> [parameter file...]\n"
  "   [-parameter[=]value...]\n"
  "      <match file in>:  FITS file with binary tables produced by WCSFoF\n"
  "      <match file out>:  FITS file with binary tables, augmented with match and color info\n"
  "      <mag file out>: Name of output FITS catalog file to hold magnitudes\n"
  "      Program parameters specified as command-line options or read from\n"
  "          parameter file(s) specified on cmd line";

// Temporary documentation:
// parameter "useInstrumentList" is a string that has format like
// <instrument name>,<instrument name>, ...
// and gives regexes that match all instrument that can be used to
// construct magnitudes and colors.  Empty list means use anything.
// Note that this is assuming that regexes do not include = or , characters.
// Whitespace will be stripped from edges of each name.

// Parameter "bands" is a comma-separated list of band names for which
// we want to output a mean magnitude.
// One band can come from multiple instruments.
// Output magnitude catalog will have one column for each distinct band.

// Parameter "colors" is a string saying what colors will be calculated.  Format is
// <band name> - <band name> @ <file name> [, <band name> - ....]
// Each comma-separated specifier gives the two bands that will be differenced to produce the
// color that is written into a FITS table with given file name.

// If any of the photometric maps have color terms, they'll use the FIRST color on the list.
// It will be calculated iteratively from the incoming data.  If requisite bands are
// not present, the mags that have color terms are flagged as invalid.

// Parameter "skipFile" is name of file containing (extension, object) number pairs that will
// be ignored in making colors.

// Collect the regular expressions that will match instruments that are
// assign to each magnitude band
struct Band {
  Band(): number(-1), hasData(false) {}
  string name;
  int number;
  bool hasData;
};

class Color {
public:
  Color(): rowCount(0), extensionNumber(-1) {}
  string band1;
  string band2;
  string colorSpec;
  string filename;
  int bandNumber1;
  int bandNumber2;
  FTable data;
  long rowCount;	// Number of objects in this color's table so far
  int exposureNumber;	// extension number assigned to this color's catalog;
  int extensionNumber;	// extension number assigned to this color's catalog;
};

// Collect this information for each desired detection in the input catalogs
class MagPoint {
public:
  double ra;
  double dec;
  double magIn;
  double magErr;
  photometry::PhotoArguments args;
  photometry::SubMap* photomap;
};

// Routine to return clipped mean & uncertainty in mag for single band, with sigma clipping at threshold.
// Also accumulate weighted sums of ra, dec for un-clipped points.
void clipMeanAndError(vector<MagPoint>& vmp, double color,
		      double threshold, 
		      double& mag, double& magErr,
		      double& sumwRA, double& sumwDec, double& sumw);

int
main(int argc, char *argv[])
{
  double clipThresh;
  double maxMagError;
  double globalSysError;
  string useInstruments;
  string bandspec;
  string colors;
  string wcsFiles;
  string photoFiles;
  string skipFile;
  Pset parameters;
  {
    const int def=PsetMember::hasDefault;
    const int low=PsetMember::hasLowerBound;
    const int up=PsetMember::hasUpperBound;
    const int lowopen = low | PsetMember::openLowerBound;
    const int upopen = up | PsetMember::openUpperBound;

    parameters.addMember("maxMagError",&maxMagError, def | lowopen,
			 "Only output mags with errors below this level (mmag)", 50., 0.);
    parameters.addMember("globalSysError",&globalSysError, def | lowopen,
			 "Systematic added in quadrature to mag error (mmag)", 1., 0.);
    parameters.addMember("clipThresh",&clipThresh, def | low,
			 "Magnitude clipping threshold (sigma)", 5., 2.);
    parameters.addMember("skipFile",&skipFile, def,
			 "Optional file containing extension,object pairs to ignore", "");
    parameters.addMember("bands",&bandspec, def,
			 "bands to output and instruments they come from","");
    parameters.addMember("colors",&colors, def,
			 "colors to tabulate and files they go into","");
    parameters.addMember("useInstruments",&useInstruments, def,
			 "list of instruments to use for mags",".*");
    parameters.addMember("wcsFiles",&wcsFiles, def,
			 "files holding WCS maps to override starting WCSs","");
    parameters.addMember("photoFiles",&photoFiles, def,
			 "files holding single-band photometric solutions for input catalogs","");
  }

  try {
    // Read all the command-line and parameter-file program parameters
    processParameters(parameters, usage, 3, argc, argv);
    string inputTables = argv[1];
    string outputTables = argv[2];
    string magOutFile = argv[3];

    maxMagError *= MMAG;
    globalSysError *= MMAG;
    
    /////////////////////////////////////////////////////
    // Parse all the parameters describing maps etc. 
    /////////////////////////////////////////////////////

    // Teach *MapCollection about new kinds of PhotoMaps and PixelMaps
    loadPhotoMapParser();
    loadPixelMapParser();

    // The list of instruments that we will be using to make magnitudes
    list<string> useInstrumentList = splitArgument(useInstruments);

    // Get the list of files holding astrometric maps and read them.
    list<string> wcsFileList = splitArgument(wcsFiles);
    astrometry::PixelMapCollection astromaps;
    for (auto fname : wcsFileList) {
      ifstream ifs(fname.c_str());
      if (!ifs)
	cerr << "WARNING: could not open astrometric map file " << fname << endl;
      if (!astromaps.read(ifs))
	cerr << "WARNING: could not read astrometric map file " << fname << endl;
    }

    // Read the photometric solutions now:
    list<string> photoFileList = splitArgument(photoFiles);
    PhotoMapCollection photomaps;
    for (auto fname : photoFileList) {
      ifstream ifs(fname.c_str());
      if (!ifs)
	cerr << "WARNING: could not open photometric map file " << fname << endl;
      if (!photomaps.read(ifs))
	cerr << "WARNING: could not read photometric map file " << fname << endl;
    }

    ExtensionObjectSet clipSet(skipFile);

    map<string, Band> bands;
    {
      for (auto bandName : splitArgument(bandspec)) {
	stripWhite(bandName);
	if (bandName.empty()) {
	  cerr << "Missing band in spec <" << bandspec << ">" << endl;
	  exit(1);
	}
	// Make a new band
	bands[bandName];
      }
    }

    // Once all the bands have been read, assign them numbers
    {
      int bandCounter = 0;
      for (auto& pr : bands) {
	Band& b = pr.second;
	b.name = pr.first;
	b.number = bandCounter++;
      }
    }

    // Parse specifications for color.  Make a little struct to describe each color to be
    // extracted.

    list<Color> colorList;
    {
      for (auto colorspec : splitArgument(colors)) {
	Color c;
	list<string> colorfile = split(colorspec, '@');
	bool success = (colorfile.size() ==2);
	if (success) {
	  // 2nd entry is filename
	  c.filename = colorfile.back();
	  stripWhite(c.filename);
	  // 1st entry is the specification of the desired color
	  c.colorSpec = colorfile.front();
	  stripWhite(c.colorSpec);

	  // split up the color term at -
	  list<string> twobands = split(colorfile.front(), '-');
	  if (twobands.size() == 2) {
	    c.band1 = twobands.front();
	    c.band2 = twobands.back();
	    stripWhite(c.band1);
	    stripWhite(c.band2);
	  } else {
	    success = false;
	  }
	}
	if (!success) {
	  cerr << "Error in color specification: " << colorspec << endl;
	  exit(1);
	}

	// Check that both mags in the color are known to us:
	if (!bands.count(c.band1)) {
	  cerr << "Unknown band <" << c.band1 << "> requested for color" << endl;
	  exit(1);
	}
	if (!bands.count(c.band2)) {
	  cerr << "Unknown band <" << c.band2 << "> requested for color" << endl;
	  exit(1);
	}
	c.bandNumber1 = bands[c.band1].number;
	c.bandNumber2 = bands[c.band2].number;

	colorList.push_back(c);
      } // end color spec loop
    } 

    // Pointer to the Color that will be used in photometric mappings.
    const Color* useColor = nullptr;

    /////////////////////////////////////////////////////
    //  Read in properties of all Instruments, Devices, Exposures
    //  and transfer info to the output file
    /////////////////////////////////////////////////////

    // All the names will be stripped of leading/trailing white space, and internal
    // white space replaced with single underscore - this keeps PhotoMap parsing functional.

    // Let's figure out which of our FITS extensions are Instrument or MatchCatalog
    vector<int> instrumentHDUs;
    vector<int> catalogHDUs;
    inventoryFitsTables(inputTables, instrumentHDUs, catalogHDUs);

    // Read the fields table & propagate the info,
    // and discard info as it is not used here.
    {
      NameIndex fieldNames;
      vector<astrometry::SphericalCoords*> fieldProjections;
      vector<double> fieldEpochs;
      readFields(inputTables, outputTables, fieldNames, fieldProjections,
		 fieldEpochs);
      for (auto ptr : fieldProjections) delete ptr;
    }

    // We've already opened (and overwritten) the output file:
    bool outputCatalogAlreadyOpen = true;
    
    // Read in all the instrument extensions and their device info from input
    // FITS file, save useful ones and write to output FITS file.
    vector<unique_ptr<Instrument>> instruments =
      readInstruments(instrumentHDUs, useInstrumentList, inputTables, outputTables,
		      outputCatalogAlreadyOpen);
    
    // Now we're going to get rid of instruments that are not in a desired band.
    for (auto& iptr : instruments) {
      if (iptr && bands.count(iptr->band)==0) {
	      iptr = nullptr;
      }
    }

    // Read in the table of exposures.  Will not use FitSubroutines readExposures()
    // because we need to alter this table to add entries for the mag catalog
    // and color catalogs.
    vector<Exposure*> exposures;
    int magTableExposureNumber = -1;	// This will be the exposure # of the mag catalog

    {
      FITS::FitsTable ft(inputTables, FITS::ReadOnly, "Exposures");
      FTable ff = ft.extract();
      vector<string> names;
      vector<double> ra;
      vector<double> dec;
      vector<int> fieldNumber;
      vector<int> instrumentNumber;
      vector<double> airmass;
      vector<double> exptime;
      // May have MJD and Epoch columns too
      vector<double> mjd;
      vector<string> epoch;
      vector<double> sysvar;  // Systematic error variance per exposure

      ff.readCells(names, "Name");
      /* Here I need to get rid of this column and then re-add it, because the
       * string columns are being stored as fixed-length by Python code, and
       * I may be appending new strings that are longer.  So need to redeclare the
       * column as variable-length.
       */
      ff.eraseColumn("Name");
      ff.addColumn(names,"Name");
      
      ff.readCells(ra, "RA");
      ff.readCells(dec, "Dec");
      ff.readCells(fieldNumber, "fieldNumber");
      ff.readCells(instrumentNumber, "InstrumentNumber");
      ff.readCells(airmass, "Airmass");
      ff.readCells(exptime, "Exptime");

      bool haveEpoch = true;
      try {
	ff.readCells(epoch, "Epoch");
      } catch (img::FTableNonExistentColumn& m) {
	haveEpoch = false;
      }

      bool haveMJD = true;
      try {
	ff.readCells(mjd, "MJD");
      } catch (img::FTableNonExistentColumn& m) {
	haveMJD = false;
      }

      try {
	ff.readCells(sysvar, "syserrmmag");
	for (int i=0; i<sysvar.size(); i++)
	  sysvar[i] = pow(sysvar[i]*MMAG, 2.) + globalSysError*globalSysError;
      } catch (img::FTableError& e) {
	// No per-exposure systematic errors, just use global one.
	sysvar.resize(ra.size());
	for (int i=0; i<sysvar.size(); i++)
	  sysvar[i] = globalSysError*globalSysError;
      }

      for (int i=0; i<names.size(); i++) {
	// Only create an Exposure if this exposure contains useful mag info.
	if (instrumentNumber[i]<0 || !instruments[instrumentNumber[i]] ) {
	  exposures.push_back(nullptr);
	} else {
	  spaceReplace(names[i]);
	  // The projection we will use for this exposure:
	  astrometry::Gnomonic gn(astrometry::Orientation(astrometry::SphericalICRS(ra[i]*WCS_UNIT,
									        dec[i]*WCS_UNIT)));
	  auto expo = new Exposure(names[i],gn);
	  expo->field = fieldNumber[i];
	  expo->instrument = instrumentNumber[i];
	  expo->airmass = airmass[i];
	  expo->exptime = exptime[i];
	  if (haveMJD) expo->mjd = mjd[i];
	  if (haveEpoch) expo->epoch = epoch[i];
	  expo->photometricVariance = sysvar[i];
	  exposures.push_back(expo);
	}
      }

      // Add an exposure for the magnitude catalog.  Mag and color catalogs will be
      // considered tags, and have no instrument.  We will give field number, RA, and Dec of
      // the exposures as 0 since we are putting all colors into one catalog that crosses fields.
      // ???? could split colors into fields

      magTableExposureNumber = ff.nrows();
      ff.writeCell(0., "RA", magTableExposureNumber);
      ff.writeCell(0., "Dec", magTableExposureNumber);
      ff.writeCell(0, "FieldNumber", magTableExposureNumber);
      ff.writeCell(TAG_INSTRUMENT, "InstrumentNumber", magTableExposureNumber);
      ff.writeCell(magOutFile, "Name", magTableExposureNumber);
      ff.writeCell(1., "Airmass", magTableExposureNumber);
      ff.writeCell(1., "Exptime", magTableExposureNumber);
      if (haveMJD) ff.writeCell(0., "MJD", magTableExposureNumber);
      if (haveEpoch) ff.writeCell(string(), "Epoch", magTableExposureNumber);

      // And now each of the color catalogs, which will be its own exposure
      for (auto& c : colorList) {
	c.exposureNumber = ff.nrows();
	ff.writeCell(0., "RA", c.exposureNumber);
	ff.writeCell(0., "Dec", c.exposureNumber);
	ff.writeCell(0, "FieldNumber", c.exposureNumber);
	ff.writeCell(TAG_INSTRUMENT, "InstrumentNumber", c.exposureNumber);
	ff.writeCell(c.filename, "Name", c.exposureNumber);
	ff.writeCell(1., "Airmass", c.exposureNumber);
	ff.writeCell(1., "Exptime", c.exposureNumber);
	if (haveMJD) ff.writeCell(0., "MJD", c.exposureNumber);
	if (haveEpoch) ff.writeCell(string(), "Epoch", c.exposureNumber);
      }

      // Done with exposure table.  Append augmented table to output FITS file
      FITS::FitsTable out(outputTables, FITS::ReadWrite+FITS::Create,
			  "Exposures");
      out.copy(ff);
    }


    // Read info about all Extensions - we will keep the Table around.
    FTable extensionTable;
    {
      FITS::FitsTable ftExtensions(inputTables, FITS::ReadOnly, "Extensions");
      extensionTable = ftExtensions.extract();
    }
    vector<Photo::Extension*> extensions(extensionTable.nrows(),nullptr);

    // This array will give band associated with each extension.  -1 for no desired band.
    vector<int> extensionBandNumbers(extensionTable.nrows(), -1);

    for (int i=0; i<extensionTable.nrows(); i++) {
      int iExposure;
      extensionTable.readCell(iExposure, "Exposure", i);
      if (!exposures[iExposure]) continue;  // Not an extension of interest

      auto extn = new Photo::Extension;
      Exposure& expo = *exposures[iExposure];
      extensions[i] = extn;
      extn->exposure = iExposure;
      int iDevice;
      extensionTable.readCell(iDevice, "Device", i);
      extn->device = iDevice;
      extn->magshift = +2.5*log10(expo.exptime);

      Assert(expo.instrument>=0);  // This exposure should have a real instrument
      Assert(instruments[expo.instrument]);  // that we are using
      string b = instruments[expo.instrument]->band;
      Assert(bands.count(b));  // and a band in use too
      // Save away the band number for this extension
      extensionBandNumbers[i] = bands[b].number;
      // and get the name its maps should be found under.
      string mapName = expo.name + "/" 
	+ instruments[expo.instrument]->deviceNames.nameOf(extn->device);

      string s;
      extensionTable.readCell(s, "WCSIN", i);
      if (stringstuff::nocaseEqual(s, "_ICRS")) {
	// Create a Wcs that just takes input as RA and Dec in std units;
	astrometry::IdentityMap identity;
	astrometry::SphericalICRS icrs;
	extn->startWcs = std::unique_ptr<astrometry::Wcs>(
    new astrometry::Wcs(&identity, icrs, "ICRS_degrees", WCS_UNIT)
  );
      } else {
	try {
	  // See if there is a map to use from the astrometric solution files
	  // Note we are cloning instead of issuing since Extension expects to own the startWcs.
	  extn->startWcs = std::unique_ptr<astrometry::Wcs>(astromaps.cloneWcs(mapName));
	} catch (astrometry::AstrometryError& m) {
	  // else we read from the input file's starting WCS
	  istringstream iss(s);
	  astrometry::PixelMapCollection pmcTemp;
	  if (!pmcTemp.read(iss)) {
	    cerr << "Did not find WCS format for " << mapName << endl;
	    exit(1);
	  }
	  string wcsName = pmcTemp.allWcsNames().front();
	  extn->startWcs = std::unique_ptr<astrometry::Wcs>(pmcTemp.cloneWcs(wcsName));
	}
      }
      // destination projection is the Exposure projection, whose coords used for any
      // exposure-level magnitude corrections
      extn->startWcs->reprojectTo(*expo.projection);

      // Get the photometric solution too
      try {
	extn->map = photomaps.issueMap(mapName);
      } catch (photometry::PhotometryError& m) {
	// No map for available for this exposure, ignore its data henceforth
	cerr << "WARNING: no photometry map available for " << mapName << endl;
	delete extn;
	extensions[i] = nullptr;
	continue;
      }
      if (extn->map->needsColor() && !useColor) {
	// At least one photomap needs color.  Make sure we have a color!
	if (colorList.empty()) {
	  cerr << "PhotoMap <" << mapName << "> has color term but no colors"
	       << " are specified." << endl;
	  exit(1);
	} else {
	  // Set up the first color to be the one used in color terms.
	  useColor = &(colorList.front());
	}
      }
    }

    // Add the mag & color catalogs to the Extension table and assign them extension numbers
    int magTableExtensionNumber = extensionTable.nrows();
    extensionTable.writeCell(magOutFile, "FILENAME", magTableExtensionNumber);
    extensionTable.writeCell(magTableExposureNumber, "EXPOSURE", magTableExtensionNumber);
    // Given a nonsense device number:
    extensionTable.writeCell(-1, "DEVICE", magTableExtensionNumber);
    // Assuming the table is in the first non-primary extension of the file:
    extensionTable.writeCell(1, "EXTENSION", magTableExtensionNumber);
    extensionTable.writeCell(string("_ICRS"), "WCSIN", magTableExtensionNumber);
    // Names of required keys - some are null as there is no relevant file
    extensionTable.writeCell(string("RA"), "XKEY", magTableExtensionNumber);
    extensionTable.writeCell(string("Dec"), "YKEY", magTableExtensionNumber);
    extensionTable.writeCell(string("_ROW"), "IDKEY", magTableExtensionNumber);
    extensionTable.writeCell(string(""), "ERRKEY", magTableExtensionNumber);

    // Get a vector of all column names in this table
    list<string> extantCols;
    {
      vector<string> vCols = extensionTable.listColumns();
      for (auto colname : vCols) extantCols.push_back(colname);
    }

    /** ?? problem writing these because sometimes they are float, sometimes double???
	but clearly some other issue with the string-valued ones...
    // Assign null weights, magweights, keys if they exist in the table
    if (regexMatchAny(extantCols, "Weight"))
      extensionTable.writeCell(0., "Weight", magTableExtensionNumber);
    if (regexMatchAny(extantCols, "MagWeight"))
      extensionTable.writeCell(0., "MagWeight", magTableExtensionNumber);
    if (regexMatchAny(extantCols, "MagKey"))
      extensionTable.writeCell(string(""), "MagKey", magTableExtensionNumber);
    if (regexMatchAny(extantCols, "MagErrKey"))
      extensionTable.writeCell(string(""), "MagErrKey", magTableExtensionNumber);
    ***/
    for (auto& c : colorList) {
      c.extensionNumber = extensionTable.nrows();
      extensionTable.writeCell(c.filename, "FILENAME", c.extensionNumber);
      extensionTable.writeCell(c.exposureNumber, "EXPOSURE", c.extensionNumber);
      // Given a nonsense Device number:
      extensionTable.writeCell(-1, "DEVICE", c.extensionNumber);
      // Assuming the table is in the first non-primary extension of the file:
      extensionTable.writeCell(1, "EXTENSION", c.extensionNumber);
      extensionTable.writeCell(string("_ICRS"), "WCSIN", c.extensionNumber);
      // Names of required keys - some are null as there is no relevant file
      extensionTable.writeCell(string("RA"), "XKEY", c.extensionNumber);
      extensionTable.writeCell(string("Dec"), "YKEY", c.extensionNumber);
      extensionTable.writeCell(string("_ROW"), "IDKEY", c.extensionNumber);
      extensionTable.writeCell(string(""), "ERRKEY", c.extensionNumber);

      /***
      // Assign null weights, magweights, keys if they exist in the table
      if (regexMatchAny(extantCols, "Weight"))
	extensionTable.writeCell(0., "Weight", c.extensionNumber);
      if (regexMatchAny(extantCols, "MagWeight"))
	extensionTable.writeCell(0., "MagWeight", c.extensionNumber);
      if (regexMatchAny(extantCols, "MagKey"))
	extensionTable.writeCell(string(""), "MagKey", c.extensionNumber);
      if (regexMatchAny(extantCols, "MagErrKey"))
	extensionTable.writeCell(string(""), "MagErrKey", c.extensionNumber);
      ***/
      if (regexMatchAny(extantCols, "ColorExpression"))
	extensionTable.writeCell(string("COLOR"), "ColorExpression", c.extensionNumber);
    }      

    // The extension table is now augmented.  Write it to the output FITS file:
    {
      FitsTable out(outputTables, 
		    outputCatalogAlreadyOpen ? FITS::ReadWrite+FITS::Create : FITS::OverwriteFile,
		    "Extensions");
      outputCatalogAlreadyOpen = true;
      out.copy(extensionTable);
    }

    //////////////////////////////////////////////////////////
    // Read in all the data
    //////////////////////////////////////////////////////////

    // Start by reading all matched catalogs,
    // telling each Extension which objects it should retrieve from its catalog

    // Keep a set of desired object IDs from each extension
    vector<set<long> > desiredObjects(extensions.size());

    for (int icat = 0; icat < catalogHDUs.size(); icat++) {
      FITS::FitsTable ft(inputTables, FITS::ReadOnly, catalogHDUs[icat]);
      FTable ff = ft.use();
      {
	// Only stellar affinity has multiband matches:
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
	if (clipSet(extn[j],obj[j])) continue;	// Skip unwanted detections.
	desiredObjects[extn[j]].insert(obj[j]);
      }
    } // End loop over input matched catalogs

    // One catalog for each extension, with object ID's as keys
    vector<map<long, MagPoint> > pointMaps(extensions.size());

    // Now loop over all original catalog bintables, reading the desired rows
    // and collecting needed information into the Detection structures
    // Should be safe to multithread this loop as different threads write
    // only to distinct parts of memory.  Protect the FITS table read though.
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic,60) num_threads(CATALOG_THREADS)
#endif
    for (int iext = 0; iext < extensions.size(); iext++) {
      if (!extensions[iext]) continue; // Skip unused 
      // Relevant structures for this extension
      auto& extn = *extensions[iext];
      auto& expo = *exposures[extn.exposure];

      string filename;
      extensionTable.readCell(filename, "FILENAME", iext);
      /**/if (iext%50==0) cerr << "# Reading object catalog " << iext
			       << "/" << extensions.size()
			       << " from " << filename 
			       << " seeking " << desiredObjects[iext].size()
			       << " objects" << endl;
      int hduNumber;
      extensionTable.readCell(hduNumber, "EXTENSION", iext);
      string xKey;
      extensionTable.readCell(xKey, "XKEY", iext);
      string yKey;
      extensionTable.readCell(yKey, "YKEY", iext);
      string idKey;
      extensionTable.readCell(idKey, "IDKEY", iext);
      string magKey;
      extensionTable.readCell(magKey, "MAGKEY", iext);
      string magErrKey;
      extensionTable.readCell(magErrKey, "MAGERRKEY", iext);
      double weight;
      extensionTable.readCell(weight, "MAGWEIGHT", iext);

      astrometry::Wcs* startWcs = extn.startWcs.get();
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

      // Be willing to get an element of array-valued bintable cell
      // for mag or magerr.  Syntax would be
      // MAGAPER[4]  to get 4th (0-indexed) element of MAGAPER column
      int magKeyElement = elementNumber(magKey);
      int magErrKeyElement = elementNumber(magErrKey);
      neededColumns.push_back(magKey);
      neededColumns.push_back(magErrKey);

      img::FTable ff;
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

      bool magColumnIsDouble = isDouble(ff, magKey, magKeyElement);
      bool magErrColumnIsDouble = isDouble(ff, magErrKey, magErrKeyElement);

      for (long irow = 0; irow < ff.nrows(); irow++) {
	set<long>::iterator pr = desiredObjects[iext].find(id[irow]);
	if (pr == desiredObjects[iext].end()) continue; // Not a desired object
	// Remove this object's ID from the list of those being sought.
	desiredObjects[iext].erase(pr);

	// Read the PhotometryArguments for this Detection:
	MagPoint mp;

	ff.readCell(mp.args.xDevice, xKey, irow);
	ff.readCell(mp.args.yDevice, yKey, irow);
	astrometry::SphericalICRS sky = startWcs->toSky(mp.args.xDevice, mp.args.yDevice);
	sky.getLonLat(mp.ra, mp.dec);
	startWcs->toWorld(mp.args.xDevice, mp.args.yDevice,
			  mp.args.xExposure, mp.args.yExposure);
	mp.args.color = 0.;	// Initialize with zero color.


	// Get the mag input and its error, adjust for exposure time
	mp.magIn = getTableDouble(ff, magKey, magKeyElement, magColumnIsDouble,irow)
	  + extn.magshift;
	mp.magErr = getTableDouble(ff, magErrKey, magErrKeyElement, magErrColumnIsDouble,irow);
	mp.magErr = sqrt(mp.magErr*mp.magErr + expo.photometricVariance);
	// Add exposure systematic error in quadrature
	mp.photomap = photomap;

	// If we have any inf or nan magnitudes, do not propagate this point
	if (isnan(mp.magIn) || isinf(mp.magIn) || isnan(mp.magErr) || isinf(mp.magErr))
	  continue;
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

    // Make output tables
    FitsTable magFitsTable(magOutFile, FITS::OverwriteFile);
    FTable magTable = magFitsTable.use();
    magTable.clear();	// should already be empty though
    {
      vector<double> dummy;
      magTable.addColumn(dummy, "RA");
      magTable.addColumn(dummy, "Dec");
      for (auto iBand : bands) {
	magTable.addColumn(dummy, iBand.second.name);
	magTable.addColumn(dummy, iBand.second.name+"_ERR");
      }
    }
    long magRowCounter = 0;	// Count of objects with assigned magnitudes

    // Make a catalog for each color being calculated
    for (auto c : colorList) {
      c.data.header()->append("COLOR_ID",c.colorSpec);
      vector<double> dummy;
      c.data.addColumn(dummy, "RA");
      c.data.addColumn(dummy, "Dec");
      c.data.addColumn(dummy, colorColumnName);
      c.data.addColumn(dummy, colorErrorColumnName);
    }

    // Iterate over the match catalogs again, this time calculating mags & colors,
    // and adding the magnitude and color table entries to the matches
    // Also will append each one to the output file
    for (int icat = 0; icat < catalogHDUs.size(); icat++) {
      /**/cerr << "Ready to start calculating catalog extension " << inputTables << " " 
	       << catalogHDUs[icat] << endl;
      FITS::FitsTable ft(inputTables, FITS::ReadOnly, catalogHDUs[icat]);
      // Append extension to output and adopt the input data
      FITS::FitsTable out(outputTables,
			  FITS::ReadWrite+FITS::Create,
			  -1);
      out.setName(ft.getName());
      FTable ff = ft.extract();
      out.adopt(ff);
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

      // Place to collect info on all observations of this match, sorted into bands:
      vector<vector<MagPoint> > bandPoints(bands.size());
      
      int lastSeq = 0;
      for (int j = 0; j<=extnIn.size(); j++) {
	if (j==extnIn.size() || seqIn[j]==0) {
	  // We have just completed processing of a match.
	  // Calculate its magnitudes, iterating if there is a color term
	  DVector mags(bands.size());
	  DVector magErrors(bands.size());
	  vector<bool> magValid(bands.size(), false);
	  bool hasMagnitudes = false;
	  // Accumulators for weighted mean of position:
	  double sumwRA;
	  double sumwDec;
	  double sumw;
	  const int colorIterations = useColor ? 3 : 1;
	  for (int colorIter=0; colorIter<colorIterations; colorIter++) {
	    double matchColor = 0.;
	    if (colorIter > 0) {
	      if ( magValid[useColor->bandNumber1] && magValid[useColor->bandNumber2]) {
		// We have the two mags needed for a color:
		matchColor = 
		  mags[useColor->bandNumber1] - mags[useColor->bandNumber2];
	      } else {
		// Cannot use this object if we do not have requisite mags for color terms
		hasMagnitudes = false;
		break;
	      }
	    }

	    hasMagnitudes = false;
	    sumwRA = sumwDec = sumw = 0.;  // clear these accumulators
	    // Calculate mean and error on mag in each band, and accumulate RA/Dec
	    for (int i=0; i<mags.size(); i++) {
	      clipMeanAndError(bandPoints[i], matchColor, clipThresh,
			       mags[i], magErrors[i],
			       sumwRA, sumwDec, sumw);
	      magValid[i] = magErrors[i]>0. && magErrors[i] <= maxMagError;
	      if (magValid[i]) hasMagnitudes = true;
	    }
	  } // End color iteration loop
	  
	  if (hasMagnitudes) {
	    // Get RA & Dec in std units
	    double ra = sumwRA / (sumw * WCS_UNIT);
	    double dec= sumwDec / (sumw * WCS_UNIT);
	    magTable.writeCell(ra, "RA", magRowCounter);
	    magTable.writeCell(dec, "Dec", magRowCounter);
	    // write mags & colors if there are any of use here
	    for (auto iBand : bands) {
	      int index = iBand.second.number;
	      if (magValid[index]) {
		magTable.writeCell( mags[index],
				    iBand.second.name, magRowCounter);
		magTable.writeCell( magErrors[index],
				    iBand.second.name+"_ERR", magRowCounter);
	      } else {
		// Mags with too-large errors get no-data marker
		magTable.writeCell( NO_MAG_DATA,
				    iBand.second.name, magRowCounter);
		magTable.writeCell( NO_MAG_DATA,
				    iBand.second.name+"_ERR", magRowCounter);
	      }
	    } // end mag band loop

	    //  add a seq for the mag catalog entry
	    seqOut.push_back(++lastSeq);
	    extnOut.push_back(magTableExtensionNumber);
	    objOut.push_back(magRowCounter++);
	    
	    //       for each possible color
	    for (auto& c : colorList) {
	      // do we have mags needed for this color?
	      if ( magValid[c.bandNumber1] && magValid[c.bandNumber2]) {
		// yes: write to color catalog
		c.data.writeCell(mags[c.bandNumber1] - mags[c.bandNumber2],
				 colorColumnName, 
				 c.rowCount);
		c.data.writeCell(hypot(magErrors[c.bandNumber1], magErrors[c.bandNumber2]),
				 colorErrorColumnName, 
				 c.rowCount);
		c.data.writeCell(ra, "RA", c.rowCount);
		c.data.writeCell(dec, "Dec", c.rowCount);

		//   add seq entry for color catalog
		seqOut.push_back(++lastSeq);
		extnOut.push_back(c.extensionNumber);
		objOut.push_back(c.rowCount++);
	      } // End processing for a valid color
	    } // end color loop
	  } // end hasMagnitudes

	  //  clear data collection for a new object
	  for (int i = 0; i<bandPoints.size(); i++) bandPoints[i].clear();
	} // End block for completing a matched set

	if (j==extnIn.size()) break;  // Quit if we don't have new data

	// Propagate every detection to output table
	seqOut.push_back(seqIn[j]);
	extnOut.push_back(extnIn[j]);
	objOut.push_back(objIn[j]);

	lastSeq = seqIn[j];

	// See if this detection has a mag we're interested in - if so,
	// add its data for that band.
	int iBand = extensionBandNumbers[extnIn[j]];
	if (iBand>=0 &&
	    pointMaps[extnIn[j]].count(objIn[j]) >0 ) {
	  bandPoints[iBand].push_back(pointMaps[extnIn[j]][objIn[j]]);
	}
      } // End object loop for this match catalog.
	
      // Write back the new sequences
      ff.writeCells(seqOut, "SequenceNumber", 0);
      ff.writeCells(extnOut, "Extension", 0);
      ff.writeCells(objOut, "Object", 0);

    } // End loop over input matched catalogs

    // Write each of the color catalogs out to a file
    for (auto& c : colorList) {
      /**/cerr << "Ready to write color catalog " << c.filename <<endl;
      FitsTable ft(c.filename, FITS::OverwriteFile);
      ft.copy(c.data);
    }

    // Cleanup:

    // Get rid of extensions
    for (auto i : extensions) delete i;
    // Get rid of exposures
    for (auto i : exposures) delete i;

  } catch (std::runtime_error& m) {
    quit(m,1);
  }
}

void 
clipMeanAndError(vector<MagPoint>& vmp, double color,
		 double threshold,
		 double& mag, double& magErr,
		 double& sumwRA, double& sumwDec, double& sumwradec) {

  vector<bool> use(vmp.size(), true);
  vector<double> x(vmp.size());
  vector<double> w(vmp.size());
  for (int i=0; i<vmp.size(); i++) {
    MagPoint& mp = vmp[i];
    mp.args.color = color;
    x[i] = mp.photomap->forward(mp.magIn, mp.args);
    double e = mp.magErr * abs(mp.photomap->derivative(mp.magIn, mp.args));
    w[i] = 1./(e*e);
  }

  bool anyclip = true;
  const int MaxIterations = MAX(int(x.size()/3),3);
  int iteration = 0;
  double sumw, sumxw;
  double threshsq = threshold*threshold;
  do {
    iteration++;
    sumw=0.;
    sumxw = 0.;
    for (int i=0; i<x.size(); i++) {
      if (use[i]) {
	sumw += w[i];
	sumxw += x[i]*w[i];
      }
    }
    mag = magErr = 0.;
    if (sumxw <=0. || iteration > MaxIterations) {
      // Return with signal for no valid data or too many clipping iterations
      return;
    } else {
      mag = sumxw / sumw;
      magErr = sqrt(1./sumxw);
    }
    anyclip = false;
    double worstsq = 0.;
    int worsti=-1;
    for (int i=0; i<x.size(); i++) {
      if (use[i] && 
	  (x[i]-mag)*(x[i]-mag)*w[i] > worstsq) {
	worstsq = (x[i]-mag)*(x[i]-mag)*w[i];
	worsti = i;
      }
    }
    if (worstsq > threshsq) {
      anyclip = true;
      use[worsti] = false;
    }
  } while (anyclip);

  // Accumulate unclipped detections' positions
  for (int i=0; i<x.size(); i++) {
    if (use[i]) {
      sumwRA += w[i] * vmp[i].ra;
      sumwDec += w[i] * vmp[i].dec;
      sumwradec += w[i];
    }
  }

  return;
}
