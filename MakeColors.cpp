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
#include "FitsImage.h"

#include "FitSubroutines.h"

using namespace std;
using namespace stringstuff;
using namespace photometry;
using img::FTable;
using FITS::FitsTable;

string usage=
  "MakeColors: Merge magnitudes for matched detections to create new\n"
  "            catalogs giving magnitudes and colors for each match.\n"
  "usage: MakeColors <match file in> <match file out> <mag file out> [parameter file...]\n"
  "      <match file in>:  FITS file with binary tables produced by WCSFoF\n"
  "      <match file out>:  FITS file with binary tables, augmented with match and color info\n"
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

// parameter clipFile is name of file containing (extension, object) number pairs that will
// be ignored (clipped) in making colors.


// Routine to return clipped mean & uncertainty given data and weights w:
void clipMeanAndError(vector<double>& x, vector<double>& w, double& mean, double& err,
		      double threshold);


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
  Color(): rowCount(0), extensionNumber(-1) {}
  string band1;
  string band2;
  string colorSpec;
  string filename;
  double offset;
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
  double mag;
  double magErr;
};


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
  string skipFile;
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
    parameters.addMember("skipFile",&skipFile, def,
			 "Optional file containing extension,object pairs to ignore", "");
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

  try {
    // Read all the command-line and parameter-file program parameters
    processParameters(parameters, usage, 3, argc, argv);
    string inputTables = argv[1];
    string outputTables = argv[2];
    string magOutFile = argv[3];

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
    RegexReplacements instrumentTranslator = parseTranslator(renameInstruments,
							     "renameInstruments");

    // Get the list of files holding astrometric maps and read them.
    list<string> wcsFileList = splitArgument(wcsFiles);
    astrometry::PixelMapCollection astromaps;
    for (list<string>::const_iterator i = wcsFileList.begin();
	 i != wcsFileList.end();
	 ++i) {
      ifstream ifs(i->c_str());
      if (!ifs)
	cerr << "WARNING: could not open astrometric map file " << *i << endl;
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
      if (!ifs)
	cerr << "WARNING: could not open photometric map file " << *i << endl;
      if (!photomaps.read(ifs))
	cerr << "WARNING: could not read photometric map file " << *i << endl;
    }

    ExtensionObjectSet clipSet(skipFile);

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
    //  and transfer info to the output file
    /////////////////////////////////////////////////////

    // All the names will be stripped of leading/trailing white space, and internal
    // white space replaced with single underscore - this keeps PhotoMap parsing functional.

    // First time we write to output file, overwrite the whole file:
    bool outputIsOpen = false;
    
    // Let's figure out which of our FITS extensions are Instrument or MatchCatalog
    vector<int> instrumentHDUs;
    vector<int> catalogHDUs;
    {
      // Find which extensions are instrument tables.  Copy everything except
      // Exposures, Extensions, and MatchCatalog tables to the output file.
      FITS::FitsFile ff(inputTables);
      for (int i=1; i<ff.HDUCount(); i++) {
	string hduName;
	{
	  FITS::Hdu h(inputTables, FITS::HDUAny, i);
	  hduName = h.getName();
	}
	if (stringstuff::nocaseEqual(hduName, "Instrument")) {
	  instrumentHDUs.push_back(i);
	  // Write all instrument info to output:
	  FITS::FitsTable ft(inputTables, FITS::ReadOnly, i);
	  // Append new table to output:
	  FITS::FitsTable out(outputTables,
			      outputIsOpen ? FITS::ReadWrite+FITS::Create : FITS::OverwriteFile,
			      -1);
	  outputIsOpen = true;
	  out.setName("Instrument");
	  Assert( stringstuff::nocaseEqual(ft.getName(), "Instrument"));
	  FTable ff=ft.extract();
	  out.adopt(ff);
	} else if (stringstuff::nocaseEqual(hduName, "MatchCatalog")) {
	  catalogHDUs.push_back(i);
	} else if (stringstuff::nocaseEqual(hduName, "Extensions")
		   || stringstuff::nocaseEqual(hduName, "Exposures") ) {
	  // Do nothing - we will alter these and write them back later.
	} else {
	  // Append any other tables to the output
	  FITS::FitsTable ft(inputTables, FITS::ReadOnly, i);
	  FITS::FitsTable out(outputTables, 
			      outputIsOpen ? FITS::ReadWrite+FITS::Create : FITS::OverwriteFile,
			      -1);
	  outputIsOpen = true;
	  out.setName(ft.getName());
	  FTable ff=ft.extract();
	  out.adopt(ff);
	}
      }
    }
    
    /**/cerr << "Found " << instrumentHDUs.size() << " instrument HDUs" 
	     << " and " << catalogHDUs.size() << " catalog HDUs" << endl;

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
	  instrumentBandNumbers[instrumentNumber] = iBand->second.number;
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
      /**/cerr << "Done with instrument " << instrumentName 
	       << " band " << instrumentBandNumbers[instrumentNumber] <<  endl;
    } // Done with instrument loop

    // Read in the table of exposures
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
      ff.readCells(names, "Name");
      ff.readCells(ra, "RA");
      ff.readCells(dec, "Dec");
      ff.readCells(fieldNumber, "fieldNumber");
      ff.readCells(instrumentNumber, "InstrumentNumber");

      for (int i=0; i<names.size(); i++) {
	// Only create an Exposure if this exposure contains no useful mag info.
	if (instrumentNumber[i]<0 || instrumentBandNumbers[instrumentNumber[i]] < 0) {
	  exposures.push_back(0);
	} else {
	  spaceReplace(names[i]);
	  // The projection we will use for this exposure:
	  astrometry::Gnomonic gn(astrometry::Orientation(astrometry::SphericalICRS(ra[i]*DEGREE,
										    dec[i]*DEGREE)));
	  Exposure* expo = new Exposure(names[i],gn);
	  expo->field = fieldNumber[i];
	  expo->instrument = instrumentNumber[i];
	  exposures.push_back(expo);
	  /**/cerr << "Exposure " << names[i] 
		   << " using instrument " << expo->instrument 
		   << " band " << instrumentBandNumbers[expo->instrument]
		   << endl;
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

      // And now each of the color catalogs, which will be its own exposure
      for (list<Color>::iterator iColor = colorList.begin();
	   iColor != colorList.end();
	   ++iColor) {
	iColor->exposureNumber = ff.nrows();
	ff.writeCell(0., "RA", iColor->exposureNumber);
	ff.writeCell(0., "Dec", iColor->exposureNumber);
	ff.writeCell(0, "FieldNumber", iColor->exposureNumber);
	ff.writeCell(TAG_INSTRUMENT, "InstrumentNumber", iColor->exposureNumber);
	ff.writeCell(iColor->filename, "Name", iColor->exposureNumber);
      }

      // Done with exposure table.  Append augmented table to output FITS file
      FITS::FitsTable out(outputTables, 
			  outputIsOpen ? FITS::ReadWrite+FITS::Create : FITS::OverwriteFile,
			  "Exposures");
      outputIsOpen = true;
      out.copy(ff);
    }


    // Read info about all Extensions - we will keep the Table around.
    vector<Extension*> extensions;
    FTable extensionTable;
    {
      FITS::FitsTable ftExtensions(inputTables, FITS::ReadOnly, "Extensions");
      extensionTable = ftExtensions.extract();
    }

    // This array will give band associated with each extension.  -1 for no desired band.
    vector<int> extensionBandNumbers(extensionTable.nrows(), -1);


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
	  // Note we are cloning instead of issuing since Extension expects to own the startWcs.
	  extn->startWcs = astromaps.cloneWcs(mapName);
	} catch (astrometry::AstrometryError& m) {
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

    // Add the mag & color catalogs to the Extension table and assign them extension numbers
    int magTableExtensionNumber = extensionTable.nrows();
    extensionTable.writeCell(magOutFile, "Filename", magTableExtensionNumber);
    extensionTable.writeCell(magTableExposureNumber, "ExposureNumber", magTableExtensionNumber);
    // Given a nonsense filenumber and Device number:
    extensionTable.writeCell(-1, "FileNumber", magTableExtensionNumber);
    extensionTable.writeCell(-1, "DeviceNumber", magTableExtensionNumber);
    // Assuming the table is in the first non-primary extension of the file:
    extensionTable.writeCell(1, "HDUNumber", magTableExtensionNumber);
    extensionTable.writeCell(string("ICRS"), "WCS", magTableExtensionNumber);
    // Names of required keys - some are null as there is no relevant file
    extensionTable.writeCell(string("RA"), "XKEY", magTableExtensionNumber);
    extensionTable.writeCell(string("Dec"), "YKEY", magTableExtensionNumber);
    extensionTable.writeCell(string("_ROW"), "IDKEY", magTableExtensionNumber);
    extensionTable.writeCell(string(""), "ERRKEY", magTableExtensionNumber);

    // Assign null weights, magweights, keys if they exist in the table
    // Start by getting a vector of all column names in this table
    list<string> extantCols;
    {
      vector<string> vCols = extensionTable.listColumns();
      for (int j=0; j<vCols.size(); j++) extantCols.push_back(vCols[j]);
    }

    if (regexMatchAny(extantCols, "Weight"))
      extensionTable.writeCell(0., "Weight", magTableExtensionNumber);
    if (regexMatchAny(extantCols, "MagWeight"))
      extensionTable.writeCell(0., "MagWeight", magTableExtensionNumber);
    if (regexMatchAny(extantCols, "Airmass"))
      extensionTable.writeCell(0., "Airmass", magTableExtensionNumber);
    if (regexMatchAny(extantCols, "MagKey"))
      extensionTable.writeCell(string(""), "MagKey", magTableExtensionNumber);
    if (regexMatchAny(extantCols, "MagErrKey"))
      extensionTable.writeCell(string(""), "MagErrKey", magTableExtensionNumber);

    for (list<Color>::iterator iColor = colorList.begin();
	 iColor != colorList.end();
	 ++iColor) {
      iColor->extensionNumber = extensionTable.nrows();
      extensionTable.writeCell(iColor->filename, "Filename", iColor->extensionNumber);
      extensionTable.writeCell(iColor->exposureNumber, "ExposureNumber", iColor->extensionNumber);
      // Given a nonsense filenumber and Device number:
      extensionTable.writeCell(-1, "FileNumber", iColor->extensionNumber);
      extensionTable.writeCell(-1, "DeviceNumber", iColor->extensionNumber);
      // Assuming the table is in the first non-primary extension of the file:
      extensionTable.writeCell(1, "HDUNumber", iColor->extensionNumber);
      extensionTable.writeCell(string("ICRS"), "WCS", iColor->extensionNumber);
      // Names of required keys - some are null as there is no relevant file
      extensionTable.writeCell(string("RA"), "XKEY", iColor->extensionNumber);
      extensionTable.writeCell(string("Dec"), "YKEY", iColor->extensionNumber);
      extensionTable.writeCell(string("_ROW"), "IDKEY", iColor->extensionNumber);
      extensionTable.writeCell(string(""), "ERRKEY", iColor->extensionNumber);

      // Assign null weights, magweights, keys if they exist in the table
      if (regexMatchAny(extantCols, "Weight"))
	extensionTable.writeCell(0., "Weight", iColor->extensionNumber);
      if (regexMatchAny(extantCols, "MagWeight"))
	extensionTable.writeCell(0., "MagWeight", iColor->extensionNumber);
      if (regexMatchAny(extantCols, "Airmass"))
	extensionTable.writeCell(0., "Airmass", iColor->extensionNumber);
      if (regexMatchAny(extantCols, "MagKey"))
	extensionTable.writeCell(string(""), "MagKey", iColor->extensionNumber);
      if (regexMatchAny(extantCols, "MagErrKey"))
	extensionTable.writeCell(string(""), "MagErrKey", iColor->extensionNumber);
      if (regexMatchAny(extantCols, "ColorExpression"))
	extensionTable.writeCell(string("COLOR"), "ColorExpression", iColor->extensionNumber);
    }      

    // The extension table is now augmented.  Write it to the output FITS file:
    {
      FitsTable out(outputTables, 
		    outputIsOpen ? FITS::ReadWrite+FITS::Create : FITS::OverwriteFile,
		    "Extensions");
      outputIsOpen = true;
      out.copy(extensionTable);
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
	if (clipSet(extn[j],obj[j])) continue;	// Skip unwanted detections.
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
      /**/if (iext%10==0) cerr << "# Reading object catalog " << iext
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

    // Make output tables
    FitsTable magFitsTable(magOutFile, FITS::OverwriteFile);
    FTable magTable = magFitsTable.use();
    magTable.clear();	// should already be empty though
    {
      vector<double> dummy;
      magTable.addColumn(dummy, "RA");
      magTable.addColumn(dummy, "Dec");
      for (BandMap::const_iterator iBand = bands.begin();
	   iBand != bands.end();
	   ++iBand) {
	magTable.addColumn(dummy, iBand->second.name);
	magTable.addColumn(dummy, iBand->second.name+"_ERR");
      }
    }
    long magRowCounter = 0;	// Count of objects with assigned magnitudes

    // Make a catalog for each color being calculated
    for (list<Color>::iterator i = colorList.begin();
	 i != colorList.end();
	 ++i) {
      i->data.header()->append("COLOR_ID",i->colorSpec);
      vector<double> dummy;
      i->data.addColumn(dummy, "RA");
      i->data.addColumn(dummy, "Dec");
      i->data.addColumn(dummy, colorColumnName);
      i->data.addColumn(dummy, colorErrorColumnName);
    }

    // Iterate over the match catalogs again, this time calculating mags & colors,
    // and adding the magnitude and color table entries to the matches
    // Also will append each one to the output file
    for (int icat = 0; icat < catalogHDUs.size(); icat++) {
      /**/cerr << "Ready to start calculating catalog extension" << inputTables << " " 
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
	    magTable.writeCell(sumRA, "RA", magRowCounter);
	    magTable.writeCell(sumDec, "Dec", magRowCounter);
	    // write mags & colors if there are any of use here
	    for (BandMap::const_iterator i= bands.begin();
		 i != bands.end();
		 ++i) {
	      int index = i->second.number;
	      if (magValid[index]) {
		magTable.writeCell( mags[index],
				    i->second.name, magRowCounter);
		magTable.writeCell( magErrors[index],
				    i->second.name+"_ERR", magRowCounter);
	      } else {
		// Mags with too-large errors get no-data marker
		magTable.writeCell( NO_MAG_DATA,
				    i->second.name, magRowCounter);
		magTable.writeCell( NO_MAG_DATA,
				    i->second.name+"_ERR", magRowCounter);
	      }
	    } // end mag band loop

	    //  add a seq for the mag catalog entry
	    seqOut.push_back(++lastSeq);
	    extnOut.push_back(magTableExtensionNumber);
	    objOut.push_back(magRowCounter++);
	    
	    //       for each possible color
	    for (list<Color>::iterator icolor = colorList.begin();
		 icolor != colorList.end();
		 ++icolor) {
	      //         do we have it?
	      if ( magValid[icolor->bandNumber1] && magValid[icolor->bandNumber2]) {
		//         write to color catalog
		icolor->data.writeCell(mags[icolor->bandNumber1] - mags[icolor->bandNumber2]
					+ icolor->offset,
					colorColumnName, 
					icolor->rowCount);
		icolor->data.writeCell(hypot(magErrors[icolor->bandNumber1], 
					      magErrors[icolor->bandNumber2]),
					colorErrorColumnName, 
					icolor->rowCount);
		icolor->data.writeCell(sumRA, "RA", icolor->rowCount);
		icolor->data.writeCell(sumDec, "Dec", icolor->rowCount);

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

    // Write each of the color catalogs out to a file
    for (list<Color>::iterator iColor = colorList.begin();
	 iColor!=colorList.end();
	 ++iColor) {
      /**/cerr << "Ready to write color catalog " << iColor->filename <<endl;
      FitsTable ft(iColor->filename, FITS::OverwriteFile);
      ft.copy(iColor->data);
    }

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

void 
clipMeanAndError(vector<double>& x, vector<double>& w, 
		 double& mean, double& err,
		 double threshold) {
  Assert(x.size() == w.size());

  vector<bool> use(x.size(), true);
  bool anyclip = true;
  const int MaxIterations = 10;
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
    if (sumxw <=0.) {
      mean = err = 0.;
      return;
    } else {
      mean = sumxw / sumw;
      err = sqrt(1./sumxw);
      if (iteration > MaxIterations) return;
    }
    anyclip = false;
    double worstsq = 0.;
    int worsti=-1;
    for (int i=0; i<x.size(); i++) {
      if (use[i] && 
	  (x[i]-mean)*(x[i]-mean)*w[i] > worstsq) {
	worstsq = (x[i]-mean)*(x[i]-mean)*w[i];
	worsti = i;
      }
    }
    if (worstsq > threshsq) {
      anyclip = true;
      use[worsti] = false;
    }
  } while (anyclip);
  return;
}
