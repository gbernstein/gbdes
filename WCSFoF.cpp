/* ???? TODO: 
   Allow WCS to come from name of a Wcs in a serialized file.
   Better yet, separate file for additional FITS keywords.
   Change so that FITS extensions are not overwritten
*/
// Program to match catalogs based on world coordinates.

#include "Std.h"
#include "FoF.h"
#include "FitsTable.h"
#include "Astrometry.h"
#include "ReadLdacHeader.h"
#include "NameIndex.h"
#include "Pset.h"
#include <map>
#include <iostream>
#include "PixelMapCollection.h"
#include "TPVMap.h"
#include "StringStuff.h"
#include "ExtensionAttribute.h"

using namespace std;
using namespace img;
using namespace FITS;
using stringstuff::stripWhite;

string usage="WCSFoF: Match objects using FoF algorithm on world coordinate system\n"
             "WCSFoF <field specs> <exposure specs> [parameter file] [parameter file...]\n"
	     "      <field specs>:  file holding <field name> <central coords> <extent>\n"
	     "                      one field per line, where extent is max distance   \n"
             "                      (in degrees) N, S, E, or W of any point from center.\n"
             "      <exposure specs>: FITS file with binary table of input file info...\n"
             "Program will look for additional parameters at stdin";
  
/// Struct that will hold the info about each point that matcher (and subsequent programs)
/// will need
struct Point {
  Point(double x_, double y_, long ext, long obj, long expo):
    x(2), extensionNumber(ext), objectNumber(obj), exposureNumber(expo) {
    x[0]=x_; x[1]=y_;
  }
  vector<double> x;
  long extensionNumber;  // "extension" is an individual input FITS bintable.
  long objectNumber;     // gives the object number of this Point in its input catalog
  long exposureNumber;   // Which exposure the catalog is from.
  const vector<double> getX() const {return x;}
};


// *** Here is the typedef that says we'll be collecting & matching Points
// *** using the Catalog template class in FoF.h.  Change this typedef if
// *** a better 2d matching structure is desired!
typedef fof::Catalog<Point,2> PointCat;

struct Field {
  string name;
  // Coordinate system that has lat=lon=0 at field center and does projection to use for matching
  astrometry::SphericalCoords* projection;
  double extent;
  double matchRadius;
  // Map from affinity name to its catalog of matches:
  typedef map<string, PointCat*> CatMap;
  CatMap catalogs;
  Field(): projection(0) {}
  ~Field() {
    for (CatMap::iterator i=catalogs.begin();
	 i != catalogs.end();
	 ++i) 
      delete i->second;
    if (projection) delete projection;
  }
  PointCat* catalogFor(const string& affinity) {
    if (catalogs.find(affinity)==catalogs.end()) {
      vector<double> lower(2, -extent);
      vector<double> upper(2, extent);
      // *** This line creates new point catalogs.  Might need to be changed if
      // *** we decide to use a different catalog structure.
      catalogs.insert( std::pair<string, PointCat*>(affinity,
						    new PointCat(lower,
								 upper,
								 matchRadius)));
    }
    return catalogs[affinity];
  }
private:
  // Hide:
  Field(const Field& rhs);
  void operator=(const Field& rhs);
};

// Right now a device is just a region of pixel coordinates, plus a name
struct Device: public Bounds<double> {
  string name;
};

// Instrument is a collection of Devices, with a name
struct Instrument: public NameIndex {
  Instrument(const string& name_): name(name_) {}
  string name;
  vector<Device> devices;
  Device& operator[](const string& extName) {
    int index = indexOf(extName);
    if (index<0) throw std::runtime_error("Nonexistent Device " + extName
					  + " in Instrument " + name);
    return devices[index];
  }
  Device& operator[](int index) {return devices[index];}
  int newDevice(const string& newName) {
    devices.push_back(Device());
    devices.back().name = newName;
    return NameIndex::append(newName);
  }
  // Hide this from base class
  void append(string s);
};

int
main(int argc,
     char *argv[])
{
  // Read in parameters
  double matchRadius;
  bool useAffinities;
  string outCatalogName;
  int minMatches;
  bool allowSelfMatches;
  string renameInstruments;
  string stringAttributes;
  string intAttributes;
  string doubleAttributes;

  Pset parameters;
  {
    const int def=PsetMember::hasDefault;
    const int low=PsetMember::hasLowerBound;
    const int up=PsetMember::hasUpperBound;
    const int lowopen = low | PsetMember::openLowerBound;
    const int upopen = up | PsetMember::openUpperBound;

    parameters.addMember("matchRadius",&matchRadius, def | lowopen,
			 "Matching tolerance (arcsec)", 1., 0.);
    parameters.addMember("useAffinities",&useAffinities, def,
			 "Disallow star-galaxy matches & inter-filter galaxy matches?", true);
    parameters.addMember("outName",&outCatalogName, def,
			 "filename for FITS output catalog", "match.cat");
    parameters.addMember("minMatch",&minMatches, def | low,
			 "Minimum number of detections for usable match", 2, 2);
    parameters.addMember("selfMatch",&allowSelfMatches, def,
			 "Retain matches that have 2 elements from same exposure?", false);
    parameters.addMember("renameInstruments",&renameInstruments, def,
			 "list of new names to give to instruments for maps","");
    parameters.addMember("stringAttributes",&stringAttributes, def,
			 "list of string-valued extension attributes","errKey,tpvOut");
    parameters.addMember("intAttributes",&intAttributes, def,
			 "list of string-valued extension attributes","");
    parameters.addMember("doubleAttributes",&doubleAttributes, def,
			 "list of string-valued extension attributes","weight");
  }

  ////////////////////////////////////////////////
  // Read parameters
  ////////////////////////////////////////////////

  if (argc<3) {
    cerr << usage << endl;
    cerr << "--------- Default Parameters: ---------" << endl;
    parameters.setDefault();
    parameters.dump(cerr);
    exit(1);
  }
  string fieldSpecs = argv[1];
  string exposureSpecs = argv[2];

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

  // Convert matching radius to our system units for world coords (degrees)
  matchRadius *= ARCSEC/DEGREE;

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

  try {

    ////////////////////////////////////////////////
    // Read in fields and save away orientations of tangent plane systems for each
    ////////////////////////////////////////////////
    NameIndex fieldNames;
    vector<Field*> fields;

    {
      ifstream ifs(fieldSpecs.c_str());
      if (!ifs) {
	cerr << "Can't open field specification file " << fieldSpecs << endl;
	exit(1);
      }

      string buffer;
      while (stringstuff::getlineNoComment(ifs, buffer)) {
	string name;
	astrometry::SphericalICRS pole;
	double extent;
	{
	  istringstream iss(buffer);
	  if (!(iss >> name >> pole >> extent)) {
	    cerr << "Bad field spec <" << buffer << ">" << endl;
	    exit(1);
	  }
	  // Check for duplicate field name:
	  if (fieldNames.has(name)) {
	    cerr << "Duplicate field name: " << name << endl;
	    exit(1);
	  }
	  Field* f = new Field;
	  f->name = name;
	  astrometry::Orientation orient(pole);
	  f->projection = new astrometry::Gnomonic(orient);
	  f->extent = extent;
	  f->matchRadius = matchRadius;
	  fields.push_back(f);
	  fieldNames.append(name);
	  Assert(fields.size()==fieldNames.size());
	}
      }
    } // Done with input fields file
    if (fields.empty()) {
      cerr << "No fields were given" << endl;
      exit(1);
    }

    /**/cerr << "Read fields" << endl;

    // Open table containing information on input files
    FTable fileTable;
    {
      FitsTable ft(exposureSpecs,FITS::ReadOnly);
      fileTable = ft.extract();
    }

    // Start a list of the instruments
    NameIndex instrumentNames;
    vector<Instrument> instruments;
    const int REFERENCE_INSTRUMENT=-1;
    const int TAG_INSTRUMENT=-2;
    const int NO_INSTRUMENT=-3;

    // And the exposures
    NameIndex exposureNames;
    vector<astrometry::SphericalICRS> exposurePointings;
    vector<int> exposureFields; // Record one field per exposure
    vector<int> exposureInstruments; // Record one instrument per exposure

    // And a list holding all Points being matched
    list<Point> allPoints;

    // Make a table into which we will stuff info about every extension 
    // of every catalog file we read:
    FTable extensionTable;

    const string stellarAffinity="STELLAR";

    // Now assemble all of the catalog attributes we will read from input and/or write to output:
    list<ExtensionAttributeBase*> attributes;

    ExtensionAttribute<string>* filenameAttr = 
      new ExtensionAttribute<string>("Filename", ExtensionAttributeBase::ReadWrite);
    attributes.push_back(filenameAttr);

    ExtensionAttribute<int>* extensionAttr = 
      new ExtensionAttribute<int>("Extension", ExtensionAttributeBase::ReadOnly, -1);
    attributes.push_back(extensionAttr);

    ExtensionAttribute<string>* fieldAttr = 
      new ExtensionAttribute<string>("Field", ExtensionAttributeBase::ReadOnly, "_NEAREST");
    attributes.push_back(fieldAttr);

    ExtensionAttribute<string>* exposureAttr = 
      new ExtensionAttribute<string>("Exposure", ExtensionAttributeBase::ReadOnly);
    attributes.push_back(exposureAttr);

    ExtensionAttribute<string>* instrumentAttr = 
      new ExtensionAttribute<string>("Instrument", ExtensionAttributeBase::ReadOnly);
    attributes.push_back(instrumentAttr);

    ExtensionAttribute<string>* deviceAttr = 
      new ExtensionAttribute<string>("Device", ExtensionAttributeBase::ReadOnly);
    attributes.push_back(deviceAttr);

    ExtensionAttribute<string>* raAttr = 
      new ExtensionAttribute<string>("RA", ExtensionAttributeBase::ReadOnly);
    attributes.push_back(raAttr);

    ExtensionAttribute<string>* decAttr = 
      new ExtensionAttribute<string>("Dec", ExtensionAttributeBase::ReadOnly);
    attributes.push_back(decAttr);

    ExtensionAttribute<string>* affinityAttr = 
      new ExtensionAttribute<string>("Affinity", ExtensionAttributeBase::ReadOnly,
				     stellarAffinity);
    attributes.push_back(affinityAttr);

    ExtensionAttribute<string>* selectAttr = 
      new ExtensionAttribute<string>("Select", ExtensionAttributeBase::ReadOnly);
    attributes.push_back(selectAttr);

    ExtensionAttribute<string>* star_selectAttr = 
      new ExtensionAttribute<string>("Star_select", ExtensionAttributeBase::ReadOnly);
    attributes.push_back(star_selectAttr);

    ExtensionAttribute<string>* xKeyAttr = 
      new ExtensionAttribute<string>("xKey", ExtensionAttributeBase::ReadWrite);
    attributes.push_back(xKeyAttr);

    ExtensionAttribute<string>* yKeyAttr = 
      new ExtensionAttribute<string>("yKey", ExtensionAttributeBase::ReadWrite);
    attributes.push_back(yKeyAttr);

    ExtensionAttribute<string>* idKeyAttr = 
      new ExtensionAttribute<string>("idKey", ExtensionAttributeBase::ReadWrite,
				     "_ROW");
    attributes.push_back(idKeyAttr);

    ExtensionAttribute<string>* wcsfileAttr = 
      new ExtensionAttribute<string>("WCSFile", ExtensionAttributeBase::ReadOnly);
    attributes.push_back(wcsfileAttr);

    // Extension attributes to be determined in this code and then output:
    ExtensionAttribute<int>* fileNumberAttr = 
      new ExtensionAttribute<int>("FileNumber", ExtensionAttributeBase::WriteOnly);
    attributes.push_back(fileNumberAttr);

    ExtensionAttribute<int>* hduNumberAttr = 
      new ExtensionAttribute<int>("HDUNumber", ExtensionAttributeBase::WriteOnly);
    attributes.push_back(hduNumberAttr);

    ExtensionAttribute<int>* exposureNumberAttr = 
      new ExtensionAttribute<int>("ExposureNumber", ExtensionAttributeBase::WriteOnly);
    attributes.push_back(exposureNumberAttr);

    ExtensionAttribute<int>* deviceNumberAttr = 
      new ExtensionAttribute<int>("DeviceNumber", ExtensionAttributeBase::WriteOnly);
    attributes.push_back(deviceNumberAttr);

    ExtensionAttribute<string>* wcsAttr = 
      new ExtensionAttribute<string>("WCS", ExtensionAttributeBase::WriteOnly);
    attributes.push_back(wcsAttr);

    // Now create ExtensionAttributes for any requested optional columns to be passed along
    {
      list<string> names = stringstuff::split(stringAttributes, listSeperator);
      for (list<string>::iterator i = names.begin();
	   i != names.end();
	   ++i) {
	string colName = *i;
	stringstuff::stripWhite(colName);
	if (colName.empty()) continue;
	attributes.push_back(new ExtensionAttribute<string>(colName, ExtensionAttributeBase::ReadWrite));
      }
    }
    {
      list<string> names = stringstuff::split(intAttributes, listSeperator);
      for (list<string>::iterator i = names.begin();
	   i != names.end();
	   ++i) {
	string colName = *i;
	stringstuff::stripWhite(colName);
	if (colName.empty()) continue;
	attributes.push_back(new ExtensionAttribute<int>(colName, ExtensionAttributeBase::ReadWrite));
      }
    }
    {
      list<string> names = stringstuff::split(doubleAttributes, listSeperator);
      for (list<string>::iterator i = names.begin();
	   i != names.end();
	   ++i) {
	string colName = *i;
	stringstuff::stripWhite(colName);
	if (colName.empty()) continue;
	attributes.push_back(new ExtensionAttribute<double>(colName, ExtensionAttributeBase::ReadWrite));
      }
    }

    // Create necessary columns in the Extension table:
    for (list<ExtensionAttributeBase*>::const_iterator i = attributes.begin();
	 i != attributes.end();
	 ++i) 
      (*i)->makeOutputColumn(extensionTable);

    long extensionNumber = 0; // cumulative counter for all FITS tables read


    // If WCS is coming from a serialized PixelMapCollection, I'll keep last-used one around
    // to avoid re-parsing it all the time:
    string currentPixelMapCollectionFileName = "";
    astrometry::PixelMapCollection* inputPmc=0;

    // Loop over input files
    for (int iFile = 0; iFile < fileTable.nrows(); iFile++) {

      // First read in everything of interest from the fileTable
      for (list<ExtensionAttributeBase*>::iterator i = attributes.begin();
	   i != attributes.end();
	   ++i) 
	(*i)->readInputTable(fileTable, iFile);

      string filename = filenameAttr->getValue();
      stripWhite(filename);
      if (filename.empty()) {
	cerr << "Missing filename for row number " << iFile << endl;
	exit(1);
      }

      /**/cerr << "Reading file " << iFile << "/" << fileTable.nrows()
	       << " " << filename << endl;

      // Open primary extension of the file to get its header
      img::Header primaryHeader;
      try {
	Hdu primary(filename, FITS::HDUAny, 0, FITS::ReadOnly);
	primaryHeader.copyFrom(*primary.header());
      } catch (FITSError& m) {
	// Could try non-FITS file format here...
	quit(m,1);
      }
      

      // Open the WCS file that holds additional keywords, if any
      string wcsFile = wcsfileAttr->getValue();
      stripWhite(wcsFile);
	
      bool usePixelMapCollection = false;   // True if we have a serialized PMC as input
      bool useTPVInput = false;	// True if might have serialized FITS Header with TPV as input
      bool xyAreRaDec = stringstuff::nocaseEqual(wcsFile, "ICRS")
	|| stringstuff::nocaseEqual(wcsFile, "@_ICRS");
      ifstream wcsStream;
      if (! (wcsFile.empty() || xyAreRaDec)) {
	// See if input file is a serialized PixelMapCollection
	if (wcsFile==currentPixelMapCollectionFileName) {
	  // Already have the correct opened file
	  usePixelMapCollection = true;
	} else {
	  // Try opening file as new serialized PMC
	  wcsStream.open(wcsFile.c_str());
	  if (!wcsStream.is_open()) {
	    cerr << "Could not open WCS information file <" << wcsFile 
		 << ">" << endl;
	    exit(1);
	  }
	  astrometry::PixelMapCollection* tryPMC = new astrometry::PixelMapCollection;
	  if (tryPMC->read(wcsStream)) {
	    // Successfully found serialized file:
	    if (inputPmc) delete inputPmc;
	    inputPmc = tryPMC;
	    currentPixelMapCollectionFileName = wcsFile;
	    usePixelMapCollection = true;
	  } else {
	    // Rewind input and try later as a TPV file:
	    wcsStream.seekg(0);
	    wcsStream.clear();
	    useTPVInput = true;
	  }
	}
      }

      // Look at just the specified HDU, or if extension number is negative (=default), use all possible
      // FITSTable extensions.
      int firstHdu = extensionAttr->getValue();
      int lastHdu = firstHdu;
      if (lastHdu < 0) {
	firstHdu = 1;
	FitsFile ff(filename);
	lastHdu = ff.HDUCount()-1;
      }
      bool isLDAC = false;
      
      for (int hduNumber = firstHdu; 
	   hduNumber<=lastHdu;
	   hduNumber++) {

	bool haveLDACHeader = false;
	Header localHeader = primaryHeader;
	FTable ft;
	{
	  FitsTable fitstab(filename, FITS::ReadOnly, hduNumber);
	  // If this is an LDAC image header, move on to next extension, which should be data
	  if (stringstuff::nocaseEqual(fitstab.getName(), "LDAC_IMHEAD")) {
	    isLDAC = true;
	    if (hduNumber==lastHdu) {
	      cerr << "Found only LDAC_HEADER at last HDU " << hduNumber
		   << " of file " << filename << endl;
	      exit(1);
	    }
	    continue;
	  }

	  // If this LDAC file, need to get header from previous extension
	  if (stringstuff::nocaseEqual(fitstab.getName(), "LDAC_OBJECTS")) 
	    isLDAC = true;

	  if (isLDAC) {
	    img::Header ldacHead = ReadLdacHeader(filename, hduNumber-1);
	    localHeader += ldacHead;
	  }

	  // Read the entire catalog
	  ft = fitstab.extract();
	} // done with the FITS extension

	localHeader += *ft.header();

	// Now update from the header any attributes that referred to header entries:
	for (list<ExtensionAttributeBase*>::iterator i = attributes.begin();
	     i != attributes.end();
	     ++i) 
	  (*i)->checkHeader(localHeader);
	
	if (stringstuff::nocaseEqual(idKeyAttr->getValue(), "_ROW")) {
	  // Want the input table row number as ID column, so make such a column
	  vector<long> vi(ft.nrows());
	  for (int i=0; i<vi.size(); i++) vi[i] = i;
	  ft.addColumn(vi, "_ROW");
	}

	// Assign an exposure number:
	string thisExposure = exposureAttr->getValue();
	if (!exposureNames.has(thisExposure)) {
	  // Create a new Exposure upon first encounter of it:
	  exposureNames.append(thisExposure);

	  // Assign an instrument, new one if not yet existent:
	  string thisInstrument = instrumentAttr->getValue();
	  if (thisInstrument.empty()) {
	    cerr << "Missing instrument description for file " << filename << endl;
	    exit(1);
	  }
	  int instrumentNumber=NO_INSTRUMENT;
	  if (stringstuff::nocaseEqual(thisInstrument, "REFERENCE"))
	    instrumentNumber = REFERENCE_INSTRUMENT;
	  else if (stringstuff::nocaseEqual(thisInstrument, "TAG"))
	    instrumentNumber = TAG_INSTRUMENT;
	  else {
	    instrumentTranslator(thisInstrument);
	    if (!instrumentNames.has(thisInstrument)) {
	      instrumentNames.append(thisInstrument);
	      instruments.push_back(Instrument(thisInstrument));
	    }
	    instrumentNumber = instrumentNames.indexOf(thisInstrument);
	    Assert(instrumentNumber >= 0);
	  }
	  exposureInstruments.push_back(instrumentNumber);

	  // Get center coordinates of this extension as center of exposure:
	  string thisRA = raAttr->getValue();
	  string thisDec = decAttr->getValue();
	  stripWhite(thisRA);
	  stripWhite(thisDec);
	  if (thisRA.empty() || thisDec.empty() ) {
	    if (xyAreRaDec) {
	      // OK to not have "exposure" coordinates, put in dummy:
	      thisRA = "00:00:00";
	      thisDec= "+00:00:00";
	    } else {
	      cerr << "Missing RA/Dec for exposure " << thisExposure
		   << " from file " << filename
		   << endl;
	      exit(1);
	    }
	  }
	  string coords = thisRA + " " + thisDec;
	  astrometry::SphericalICRS thisRADec;
	  istringstream iss(coords);
	  if (!(iss >> thisRADec)) {
	    cerr << "Error reading RA & Dec for file " << filename
		 << " from string <" << coords << ">" << endl;
	    exit(1);
	  }
	  exposurePointings.push_back(thisRADec);

	  // Assign a field to this exposure:
	  string thisField = fieldAttr->getValue();
	  if (thisField.empty()) {
	    cerr << "Missing field name for file " << filename << endl;
	    exit(1);
	  }
	  if (stringstuff::nocaseEqual(thisField, "_NEAREST")) {
	    // Assign exposure to nearest field:
	    Assert(!fields.empty());
	    vector<Field*>::const_iterator i = fields.begin();
	    (*i)->projection->setLonLat(0.,0.);
	    double minDistance = thisRADec.distance(*(*i)->projection);
	    thisField = (*i)->name;
	    for ( ; i != fields.end(); ++i) {
	      (*i)->projection->setLonLat(0.,0.);
	      if (thisRADec.distance(*(*i)->projection) < minDistance) {
		minDistance = thisRADec.distance(*(*i)->projection);
		thisField = (*i)->name;
	      }
	    }
	  }
	  int fieldNumber = fieldNames.indexOf(thisField);
	  if (fieldNumber < 0) {
	    cerr << "Did not find information for field <" << thisField << ">" << endl;
	    exit(1);
	  }
	  exposureFields.push_back(fieldNumber);

	  Assert(exposurePointings.size()==exposureNames.size());
	  Assert(exposurePointings.size()==exposureFields.size());
	  Assert(exposurePointings.size()==exposureInstruments.size());
	} // end creation of new Exposure information
	int exposureNumber = exposureNames.indexOf(thisExposure);

	// Get field and instrument from the exposure:
	int fieldNumber = exposureFields[exposureNumber];
	int instrumentNumber = exposureInstruments[exposureNumber];

	// Error if different extensions put same exposure in different fields or instruments:
	{
	  string field = fieldAttr->getValue();
	  if ( !field.empty() && !stringstuff::nocaseEqual(field, "_NEAREST")
	       && !stringstuff::nocaseEqual(field, fieldNames.nameOf(exposureFields[exposureNumber]))) {
	    cerr << "Conflicting field assignments in file " << filename 
		 << ":\n exposure has <" << fieldNames.nameOf(exposureFields[exposureNumber])
		 << ">\n extension has <" << field << ">"
		 << endl;
	    exit(1);
	  }
	} // done checking field agreement
	{
	  string instrument = instrumentAttr->getValue();
	  if ( !instrument.empty()) {
	    if ( stringstuff::nocaseEqual(instrument, "REFERENCE")) {
	      if (instrumentNumber != REFERENCE_INSTRUMENT ) {
		cerr << "Conflicting instrument assignments in file " << filename 
		     << ":\n extension is REFERENCE but exposure is not"
		     << endl;
		exit(1);
	      }
	    } else if (stringstuff::nocaseEqual(instrument, "TAG")) {
	      if (instrumentNumber != TAG_INSTRUMENT) {
		cerr << "Conflicting instrument assignments in file " << filename 
		     << ":\n extension is TAG but exposure is not"
		     << endl;
		exit(1);
	      }
	    } else {
	      instrumentTranslator(instrument);
	      if (!stringstuff::nocaseEqual(instrument, 
					    instrumentNames.nameOf(instrumentNumber))) {
		cerr << "Conflicting instrument assignments in file " << filename 
		     << ":\n exposure has <" << instrumentNames.nameOf(instrumentNumber) 
		     << ">\n extension has <" << instrument << ">"
		     << endl;
		exit(1);
	      }
	    }
	  }
	}  // Done checking instrument agreement

	// Assign a device
	int deviceNumber = -1;
	if (instrumentNumber >= 0) {
	  string thisDevice = deviceAttr->getValue();
	  if (thisDevice.empty()) {
	    cerr << "Empty device field for " << filename << endl;
	    exit(1);
	  }
	  if (!instruments[instrumentNumber].has(thisDevice)) {
	    instruments[instrumentNumber].newDevice(thisDevice);
	  }
	  deviceNumber = instruments[instrumentNumber].indexOf(thisDevice);
	  Assert(deviceNumber >= 0);
	}

	// Set for output the critical info about this extension:
	fileNumberAttr->setValue(iFile);
	exposureNumberAttr->setValue(exposureNumber);
	hduNumberAttr->setValue(hduNumber);
	deviceNumberAttr->setValue(deviceNumber);

	// Now get the WCS for this extension
	astrometry::Wcs* wcs = 0;
	if (xyAreRaDec) {
	  // Do nothing here
	} else if (usePixelMapCollection) {
	  Assert(exposureInstruments[exposureNumber] >= 0);
	  Assert(inputPmc);
	  Instrument& inst = instruments[exposureInstruments[exposureNumber]];
	  string wcsName = exposureNames.nameOf(exposureNumber) + "/"
	    + inst[deviceNumber].name;
	  // Replace any white space with underscore:
	  wcsName = stringstuff::regexReplace("[[:space:]]+","_",wcsName);
	  wcs = inputPmc->issueWcs(wcsName)->duplicate();
	} else if (useTPVInput) {
	  // Get the WCS from an external header file
	  Header wcsHead;
	  if (!(wcsStream >> wcsHead)) {
	    cerr << "Error reading WCS/header from file " << wcsFile << endl;
	    exit(1);
	  }
	  wcs = astrometry::readTPV(wcsHead);
	  if (!wcs) {
	    cerr << "Failed reading TPV WCS from file " << wcsFile << endl;
	    exit(1);
	  }
	} else {
	  // WCS should be in the headers we already have
	  wcs = astrometry::readTPV(localHeader);
	  if (!wcs) {
	    cerr << "Failed reading TPV WCS for HDU " << hduNumber
		 << " in file " << filename
		 << endl;
	    exit(1);
	  }
	}

	// Now save away serialized WCS, and set projection to be the field coordinates
	string wcsDump;
	if (xyAreRaDec) {
	  wcsDump = "ICRS";
	} else {
	  Assert(wcs);
	  astrometry::PixelMapCollection pmc;
	  pmc.learnWcs(*wcs);
	  ostringstream oss;
	  pmc.writeWcs(oss,wcs->getName());
	  wcsDump = oss.str();
	  wcs->reprojectTo(*fields[fieldNumber]->projection);
	}
	wcsAttr->setValue(wcsDump);

	// Record attributes of this extension to an output table
	for (list<ExtensionAttributeBase*>::iterator i = attributes.begin();
	     i != attributes.end();
	     ++i) 
	  (*i)->writeOutputTable(extensionTable, extensionNumber);

	// Go ahead and filter the input rows and get the columns we want
	{
	  string selectionExpression = selectAttr->getValue();
	  stripWhite(selectionExpression);
	  if (!selectionExpression.empty())
	    ft.filterRows(selectionExpression);
	}
	vector<double> vx;
	string key = xKeyAttr->getValue();
	try {
	  ft.readCells(vx, key);
	} catch (FTableError& m) {
	  // Trap for using float column in source file instead of long:
	  vector<float> vf;
	  ft.readCells(vf, key);
	  vx.resize(vf.size());
	  for (int i=0; i<vf.size(); i++) vx[i]=vf[i];
	}
	vector<double> vy;
	key = yKeyAttr->getValue();
	try {
	  ft.readCells(vy, key);
	} catch (FTableError& m) {
	  // Trap for using float column in source file instead of long:
	  vector<float> vf;
	  ft.readCells(vf, key);
	  vy.resize(vf.size());
	  for (int i=0; i<vf.size(); i++) vy[i]=vf[i];
	}
	vector<long> vid;
	key = idKeyAttr->getValue();
	try {
	  ft.readCells(vid, key);
	} catch (FTableError& m) {
	  // Trap for using int column in source file instead of long:
	  vector<int> vidint;
	  ft.readCells(vidint, key);
	  vid.reserve(vidint.size());
	  vid.insert(vid.begin(), vidint.begin(), vidint.end());
	}
	vector<bool> isStar(ft.nrows(), true);
	{
	  string starExpression = star_selectAttr->getValue();
	  stripWhite(starExpression);
	  if (!starExpression.empty())
	    ft.evaluate(isStar, starExpression);
	}

	
	// Select the Catalog(s) to which points from this extension will be added
	string thisAffinity = affinityAttr->getValue();
	if (thisAffinity.empty()) {
	  cerr << "Missing affinity for " << filename << endl;
	  exit(1);
	}
	PointCat* starCatalog = fields[fieldNumber]->catalogFor(stellarAffinity);
	PointCat* galaxyCatalog = (stringstuff::nocaseEqual(stellarAffinity,
							    thisAffinity) ?
				   starCatalog : 
				   fields[fieldNumber]->catalogFor(thisAffinity));
	
	// Now loops over objects in the catalog
	for (int iObj = 0; iObj < vx.size(); iObj++) {
	  double xpix = vx[iObj];
	  double ypix = vy[iObj];
	  // Expand bounds of device:
	  if (instrumentNumber >=0) 
	    instruments[instrumentNumber][deviceNumber]+=Position<double>(xpix,ypix);
	  
	  //  maps coords to field's tangent plane
	  double xw, yw;
	  if (xyAreRaDec) {
	    fields[fieldNumber]->projection->convertFrom(astrometry::SphericalICRS(vx[iObj]*DEGREE, 
										  vy[iObj]*DEGREE ));
	    fields[fieldNumber]->projection->getLonLat(xw, yw);
	    // Matching program will speak degrees, not radians:
	    xw /= DEGREE;
	    yw /= DEGREE;
	  } else {
	    Assert(wcs);
	    wcs->toWorld(xpix, ypix, xw, yw);
	  }
	  // ??? Note that frequent small memory allocations are making memory mgmt and overhead
	  // exceed the actual Point memory usage.  Also the allPoints list is really superfluous
	  // (just extra link pointers), 16B per Point.  Probably a good deal of overhead
	  // in the Catalog structures too.  In fact up to ~300 B per Point now!
	  allPoints.push_back(Point(xw, yw, extensionNumber, vid[iObj], exposureNumber));
	  if (isStar[iObj])
	    starCatalog->add(allPoints.back());
	  else
	    galaxyCatalog->add(allPoints.back());
	} // end object loop

	if (wcs) delete wcs;
	
	extensionNumber++;
      } // end input extension loop
    } // end input file loop

    if (inputPmc) delete inputPmc;
    // Delete all the ExtensionAttributes:
    for (list<ExtensionAttributeBase*>::iterator i = attributes.begin();
	 i != attributes.end();
	 ++i) 
      delete *i;

    cerr << "*** Read " << allPoints.size() << " objects" << endl;

    //  Open output catalog & put input file table into output table.
    {
      FitsTable ft(outCatalogName, FITS::ReadWrite + FITS::OverwriteFile, "Files");
      ft.copy(fileTable);
    }
    //  Write fields - field centers are given in degrees.
    {
      FitsTable ft(outCatalogName, FITS::ReadWrite + FITS::Create, "Fields");
      FTable ff=ft.use();
      vector<string> names;
      vector<double> ra;
      vector<double> dec;
      for (int i=0; i<fields.size(); i++) {
	names.push_back(fields[i]->name);
	fields[i]->projection->setLonLat(0.,0.);
	astrometry::SphericalICRS pole(*fields[i]->projection);
	double r,d;
	pole.getLonLat(r,d);
	ra.push_back( r/DEGREE);
	dec.push_back( d/DEGREE);
      }
      ff.addColumn(names, "Name");
      ff.addColumn(ra, "RA");
      ff.addColumn(dec, "Dec");
    }

    //  Write exposure information to output table
    {
      FitsTable ft(outCatalogName, FITS::ReadWrite + FITS::Create, "Exposures");
      FTable ff=ft.use();
      vector<string> names;
      vector<double> ra;
      vector<double> dec;
      for (int i=0; i<exposureNames.size(); i++) {
	names.push_back(exposureNames.nameOf(i));
	astrometry::SphericalICRS pole = exposurePointings[i];
	double r,d;
	pole.getLonLat(r,d);
	ra.push_back( r/DEGREE);
	dec.push_back( d/DEGREE);
      }
      ff.addColumn(names, "Name");
      ff.addColumn(ra, "RA");
      ff.addColumn(dec, "Dec");
      ff.addColumn(exposureFields, "FieldNumber");
      ff.addColumn(exposureInstruments, "InstrumentNumber");
    }

    //  Write instrument information to output table
    for (int i=0; i<instruments.size(); i++) {
      // Append new table for each instrument
      FitsTable ft(outCatalogName, FITS::ReadWrite + FITS::Create, -1);
      ft.setName("Instrument");
      ft.setVersion(i+1); // might use this later, or put into extension name???
      FTable ff = ft.use();
      Instrument& inst = instruments[i];
      ff.header()->replace("Name", inst.name, "Instrument name");
      ff.header()->replace("Number", i, "Instrument number");
      const int nDevices = inst.size();
      vector<string> devnames(nDevices);
      vector<double> vxmin(nDevices);
      vector<double> vxmax(nDevices);
      vector<double> vymin(nDevices);
      vector<double> vymax(nDevices);
      for (int j=0; j<nDevices; j++) {
	devnames[j] = inst[j].name;
	vxmin[j] = floor(inst[j].getXMin());
	vxmax[j] = ceil(inst[j].getXMax());
	vymin[j] = floor(inst[j].getYMin());
	vymax[j] = ceil(inst[j].getYMax());
      }
      ff.addColumn(devnames, "Name");
      ff.addColumn(vxmin, "XMin");
      ff.addColumn(vxmax, "XMax");
      ff.addColumn(vymin, "YMin");
      ff.addColumn(vymax, "YMax");
    }

    //  Write the extension catalog to output file
    {
      FitsTable ft(outCatalogName, FITS::ReadWrite + FITS::Create, "Extensions");
      ft.copy(extensionTable);
    }

    long matchCount = 0;
    long pointCount = 0;
    //  Loop over fields
    for (int iField=0; iField<fields.size(); iField++) {
      //  loop over affinities per field
      for (Field::CatMap::iterator iAffinity=fields[iField]->catalogs.begin();
	   iAffinity != fields[iField]->catalogs.end();
	   ++iAffinity) {
	const PointCat* pcat= iAffinity->second;
	cerr << "Catalog for field " << fields[iField]->name
	     << " affinity " << iAffinity->first
	     << " with " << pcat->size() << " matches "
	     << endl;
	if (pcat->empty()) continue;
	FitsTable ft(outCatalogName, FITS::ReadWrite + FITS::Create, -1);
	ft.setName("MatchCatalog");
	FTable ff=ft.use();
	ff.header()->replace("Field", fields[iField]->name, "Field name");
	ff.header()->replace("FieldNum", iField, "Field number");
	ff.header()->replace("Affinity", iAffinity->first, "Affinity name");

	vector<int> sequence;
	vector<long> extn;
	vector<long> obj;
	long matches=0;
	// Now loop through matches in this catalog
	for (PointCat::const_iterator j=pcat->begin();
	     j != pcat->end();
	     ++j) {
	  // Skip any Match that is below minimum match size
	  if ((*j)->size() < minMatches) continue;
	  bool selfMatch = false;
	  if (!allowSelfMatches) {
	    set<long> itsExposures;
	    for (list<const Point*>::const_iterator k=(*j)->begin();
		 k != (*j)->end();
		 ++k) 
	      if ( !itsExposures.insert((*k)->exposureNumber).second) {
		selfMatch = true;
		break;
	      }
	  }
	  if (selfMatch) continue;

	  // Add elements of this match to the output vectors
	  int seq=0;
	  ++matches;
	  ++matchCount;
	  for (list<const Point*>::const_iterator k=(*j)->begin();
	       k != (*j)->end();
	       ++k, ++seq) {
	    sequence.push_back(seq);
	    extn.push_back((*k)->extensionNumber);
	    obj.push_back((*k)->objectNumber);
	  }
	} // end match loop
	cerr << "...Kept " << matches << " matches with " << sequence.size()
	     << " points." << endl;
	pointCount += sequence.size();
	ff.addColumn(sequence, "SequenceNumber");
	ff.addColumn(extn, "Extension");
	ff.addColumn(obj, "Object");
      } // end Affinity loop
    } // end Field loop
    
    cerr << "Total of " << matchCount
	 << " matches with " << pointCount
	 << " points." << endl;

    // Clean up
    for (int i=0; i<fields.size(); i++)
      delete fields[i];

  } catch (std::runtime_error &m) {
    quit(m,1);
  }
}
