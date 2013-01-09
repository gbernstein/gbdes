// Program to match catalogs based on world coordinates.

#include "Std.h"
#include "FoF.h"
#include "FitsTable.h"
#include "Astrometry.h"
#include "SCAMPMap.h"
#include "ReadLdacHeader.h"
#include "NameIndex.h"
#include "Pset.h"
#include <map>
#include <iostream>

using namespace std;
using namespace img;
using namespace FITS;

string usage="WCSFoF: Match objects using FoF algorithm on world coordinate system\n"
             "WCSFoF <field specs> <exposure specs> [parameter file] [parameter file...]\n"
	     "      <field specs>:  file holding <field name> <central coords> <extent>\n"
	     "                      one field per line, where extent is max distance   \n"
             "                      (in degrees) N, S, E, or W of any point from center.\n"
             "      <exposure specs>: FITS file with binary table of input file info...";

// Strip leading and trailing white space from a string:
void stripWhite(string& s) {
  while (!s.empty() && std::isspace(s[0]))
    s.erase(0,1);
  while (!s.empty() && std::isspace(s[s.size()-1]))
    s.erase(s.size()-1);
}
// Parse the inValue of string: if it starts with @ sign, specifies to read a value
// from keyword of a header:
string maybeFromHeader(const string& inValue, const Header& h) {
  if (inValue.empty() || inValue[0]!='@') return inValue;
  string result;
  if (!h.getValue(inValue.substr(1), result)) {
    // Try reading an integer and convering to a string
    int i;
    if (!h.getValue(inValue.substr(1), i)) {
      throw std::runtime_error("Could not find string-valued header keyword " + 
			       inValue.substr(1));
    }
    ostringstream oss;
    oss << i;
    result = oss.str();
  }
  return result;
}

// Struct that will hold the info about each point that matcher (and subsequent progs)
// will need
struct Point {
  Point(double x_, double y_, long ext, long obj, long expo):
    x(2), extensionNumber(ext), objectNumber(obj), exposureNumber(expo) {
    x[0]=x_; x[1]=y_;
  }
  vector<double> x;
  long extensionNumber;
  long objectNumber;
  long exposureNumber;
  const vector<double> getX() const {return x;}
};


typedef fof::Catalog<Point,2> PointCat;

struct Field {
  string name;
  astrometry::Orientation orient;
  double extent;
  double matchRadius;
  // Map from affinity name to its catalog of matches:
  typedef map<string, PointCat*> CatMap;
  CatMap catalogs;
  ~Field() {
    for (CatMap::iterator i=catalogs.begin();
	 i != catalogs.end();
	 ++i) 
      delete i->second;
  }
  PointCat* catalogFor(const string& affinity) {
    if (catalogs.find(affinity)==catalogs.end()) {
      vector<double> lower(2, -extent);
      vector<double> upper(2, extent);
      catalogs.insert( std::pair<string, PointCat*>(affinity,
						    new PointCat(lower,
								 upper,
								 matchRadius)));
    }
    return catalogs[affinity];
  }
};

// Right now an device is just a region of pixel coordinates, plus a name
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
  }

  ////////////////////////////////////////////////
  // Read parameters
  ////////////////////////////////////////////////

  if (argc<3) {
    cerr << usage << endl;
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

  try {

    ////////////////////////////////////////////////
    // Read in fields and save away orientations of tangent plane systems for each
    ////////////////////////////////////////////////
    NameIndex fieldNames;
    vector<Field> fields;

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
	  Field f;
	  f.name = name;
	  f.orient = astrometry::Orientation(pole);
	  f.extent = extent;
	  f.matchRadius = matchRadius;
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

    /**/cerr << "Read fileTable " << fileTable.ncols() << endl;

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

    const string filenameCol="FILENAME";
    const string extensionCol = "EXTENSION";
    const string fieldCol="FIELD";
    const string exposureCol="EXPOSURE";
    const string raCol="RA";
    const string decCol="DEC";
    const string instrumentCol="INSTRUMENT";
    const string deviceCol = "DEVICE";
    const string affinityCol="AFFINITY";
    const string inSelectCol="SELECT";
    const string starSelectCol="STAR_SELECT";
    const string wcsFileCol="WCSFILE";
    const string wcsOutFileCol="WCSOUT";
    const string xkeyCol="XKEY";
    const string ykeyCol="YKEY";
    const string idkeyCol="IDKEY";
    const string errorkeyCol="ERRORKEY";
    const string weightCol="WEIGHT";

    const string globalAffinity="GLOBAL";

    /**/cerr << "Start input file loop" << endl;

    // Make a table into which we will stuff info about every extension 
    // of every catalog file we read:
    FTable extensionTable;
    {
      vector<string> vs;
      vector<int> vi;
      vector<double> vd;
      extensionTable.addColumn(vs, "Filename");
      extensionTable.addColumn(vi, "FileNumber");
      extensionTable.addColumn(vi, "HDUNumber");
      extensionTable.addColumn(vi, "ExposureNumber");
      extensionTable.addColumn(vi, "DeviceNumber");
      extensionTable.addColumn(vs, "WCS");
      extensionTable.addColumn(vs, "xKey");
      extensionTable.addColumn(vs, "yKey");
      extensionTable.addColumn(vs, "idKey");
      extensionTable.addColumn(vs, "errKey");
      extensionTable.addColumn(vd, "Weight");
      extensionTable.addColumn(vs, "WCSOut");
    }

    long extensionNumber = 0; // cumulative counter for all FITS tables read

    // Loop over input files
    for (int iFile = 0; iFile < fileTable.nrows(); iFile++) {
      string filename;
      fileTable.readCell(filename, filenameCol, iFile);
      int extension;
      fileTable.readCell(extension, extensionCol, iFile);
      string fieldName;
      fileTable.readCell(fieldName, fieldCol, iFile);
      string exposureName;
      fileTable.readCell(exposureName, exposureCol, iFile);
      string ra;
      fileTable.readCell(ra, raCol, iFile);
      string dec;
      fileTable.readCell(dec, decCol, iFile);
      string instrumentName;
      fileTable.readCell(instrumentName, instrumentCol, iFile);
      string deviceName;
      fileTable.readCell(deviceName, deviceCol, iFile);
      string wcsFile;
      fileTable.readCell(wcsFile, wcsFileCol, iFile);
      string xkey;
      fileTable.readCell(xkey, xkeyCol, iFile);
      string ykey;
      fileTable.readCell(ykey, ykeyCol, iFile);
      string idkey;
      try {
	fileTable.readCell(idkey, idkeyCol, iFile);
      } catch (FTableNonExistentColumn& m) {
	// Default to row number:
	idkey = "@_ROW";
      }

      string errorkey;
      try {
	fileTable.readCell(errorkey, errorkeyCol, iFile);
      } catch (FTableNonExistentColumn& m) {
	errorkey.clear();
      }

      string wcsOut;
      try {
	fileTable.readCell(wcsOut, wcsOutFileCol, iFile);
      } catch (FTableNonExistentColumn& m) {
	wcsOut.clear();
      }

      string weightString;
      try {
	fileTable.readCell(weightString, weightCol, iFile);
      } catch (FTableNonExistentColumn& m) {
	weightString.clear();
      }

      // The affinity and selections we can take defaults for:
      string selectionExpression;
      try {
	fileTable.readCell(selectionExpression,inSelectCol, iFile);
	stripWhite(selectionExpression);
      } catch (FTableNonExistentColumn& m) {
	selectionExpression.clear();
      }
      string affinityName;
      try {
	fileTable.readCell(affinityName, affinityCol, iFile);
	stripWhite(affinityName);
      } catch (FTableNonExistentColumn& m) {
	affinityName = globalAffinity;
      }
      string starExpression;
      try {
	fileTable.readCell(starExpression,starSelectCol, iFile);
	stripWhite(starExpression);
      } catch (FTableNonExistentColumn& m) {
	starExpression.clear();	// everything becomes a star
      }

      /**/cerr << "Read file " << iFile << "/" << fileTable.nrows()
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
      // ???? Could check to see if WCS file name is stored in a header keyword ???
      // Strip leading & trailing whitespace from wcsFile:
      stripWhite(wcsFile);
	
      bool xyAreRaDec = stringstuff::nocaseEqual(wcsFile, "ICRS")
	|| stringstuff::nocaseEqual(wcsFile, "@_ICRS");
      ifstream wcsStream;
      if (! (wcsFile.empty() || xyAreRaDec)) {
	wcsStream.open(wcsFile.c_str());
	if (!wcsStream.is_open()) {
	  cerr << "Could not open WCS information file <" << wcsFile 
	       << ">" << endl;
	  exit(1);
	}
      }

      int firstHdu = extension >= 0 ? extension : 1;
      int lastHdu = extension;
      if (extension < 0) {
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

	Header wcsHead;
	if (wcsStream.is_open()) {
	  if (!(wcsStream >> wcsHead)) {
	    cerr << "Error reading WCS/header from file " << wcsFile << endl;
	    exit(1);
	  }
	  // Header values from "WCS" file supersede those in the catalog's file:
	  localHeader += wcsHead;
	}

	string thisXKey = maybeFromHeader(xkey, localHeader);
	string thisYKey = maybeFromHeader(ykey, localHeader);
	string thisAffinity = maybeFromHeader(affinityName, localHeader);
	if (thisAffinity.empty()) thisAffinity = globalAffinity;
	// we will not look up the selection filters in the header

	string thisIdKey;
	if (stringstuff::nocaseEqual(idkey, "@_ROW")) {
	  // Want the input table row number as ID column, so make such a column
	  vector<long> vi(ft.nrows());
	  for (int i=0; i<vi.size(); i++) vi[i] = i;
	  ft.addColumn(vi, "_ROW");
	  thisIdKey = "_ROW";
	} else {
	  thisIdKey = maybeFromHeader(idkey, localHeader);
	}

	// Couple of things that WCS fitting program will want if they're here:
	string thisErrorKey = maybeFromHeader(errorkey, localHeader);
	string thisWcsOutFile = maybeFromHeader(wcsOut, localHeader);

	// Read catalog's weight: it either defaults to 1. if string is empty,
	// or it's a header keyword to look up, or it should be a number
	double weight;
	if (weightString.empty()) {
	  weight = 1.;
	} else if (weightString[0]=='@') {
	  string keyword = weightString.substr(1);
	  if (!localHeader.getValue(keyword, weight)) {
	    cerr << "Could not find weight in double-valued header keyword " 
		 << keyword << endl;
	    exit(1);
	  }
	} else {
	  istringstream iss(weightString);
	  if (!(iss >> weight)) {
	    cerr << "Bad weight value string <" << weightString << ">" << endl;
	    exit(1);
	  }
	}

	// Go ahead and filter the input rows and get the columns we want
	if (!selectionExpression.empty())
	  ft.filterRows(selectionExpression);
	vector<double> vx;
	try {
	  ft.readCells(vx, thisXKey);
	} catch (FTableError& m) {
	  // Trap for using float column in source file instead of long:
	  vector<float> vf;
	  ft.readCells(vf, thisXKey);
	  vx.resize(vf.size());
	  for (int i=0; i<vf.size(); i++) vx[i]=vf[i];
	}
	vector<double> vy;
	try {
	  ft.readCells(vy, thisYKey);
	} catch (FTableError& m) {
	  // Trap for using float column in source file instead of long:
	  vector<float> vf;
	  ft.readCells(vf, thisYKey);
	  vy.resize(vf.size());
	  for (int i=0; i<vf.size(); i++) vy[i]=vf[i];
	}
	vector<long> vid;
	try {
	  ft.readCells(vid, thisIdKey);
	} catch (FTableError& m) {
	  // Trap for using int column in source file instead of long:
	  vector<int> vidint;
	  ft.readCells(vidint, thisIdKey);
	  vid.reserve(vidint.size());
	  vid.insert(vid.begin(), vidint.begin(), vidint.end());
	}
	vector<bool> isStar(ft.nrows(), true);
	if (!starExpression.empty())
	  ft.evaluate(isStar, starExpression);

	// Determine center coordinates of this extension:
	string thisRA = maybeFromHeader(ra, localHeader);
	string thisDec = maybeFromHeader(dec, localHeader);
	stripWhite(thisRA);
	stripWhite(thisDec);
	astrometry::SphericalICRS thisRADec;
	// Will default to 0 RA, 0 Dec if not given.
	// Coords only needed if first Device for exposure or using @_NEAREST field
	if (!thisRA.empty() || !thisDec.empty()) {
	  string coords = thisRA + " " + thisDec;
	  istringstream iss(coords);
	  if (!(iss >> thisRADec)) {
	    cerr << "Error reading RA & Dec for file " << filename
		 << " from string <" << coords << ">" << endl;
	    exit(1);
	  }
	}

	// Assign a field:
	string thisField;
	if (stringstuff::nocaseEqual(fieldName, "@_NEAREST")) {
	  // Assign exposure to nearest field:
	  Assert(!fields.empty());
	  vector<Field>::const_iterator i = fields.begin();
	  double minDistance = thisRADec.distance(i->orient.getPole());
	  thisField = i->name;
	  for ( ; i != fields.end(); ++i) {
	    if (thisRADec.distance(i->orient.getPole()) < minDistance) {
	      minDistance = thisRADec.distance(i->orient.getPole());
	      thisField = i->name;
	    }
	  }
	} else {
	  thisField = maybeFromHeader(fieldName, localHeader);
	}
	int fieldNumber = fieldNames.indexOf(thisField);
	if (fieldNumber < 0) {
	  cerr << "Did not find information for field <" << thisField << ">" << endl;
	  exit(1);
	}


	// Assign an instrument, new one if not yet existent:
	string thisInstrument = maybeFromHeader(instrumentName, localHeader);
	int instrumentNumber=NO_INSTRUMENT;
	if (stringstuff::nocaseEqual(thisInstrument, "REFERENCE"))
	  instrumentNumber = REFERENCE_INSTRUMENT;
	else if (stringstuff::nocaseEqual(thisInstrument, "TAG"))
	  instrumentNumber = TAG_INSTRUMENT;
	else {
	  if (!instrumentNames.has(thisInstrument)) {
	    instrumentNames.append(thisInstrument);
	    instruments.push_back(Instrument(thisInstrument));
	  }
	  instrumentNumber = instrumentNames.indexOf(thisInstrument);
	  Assert(instrumentNumber >= 0);
	}

	// Assign an exposure number:
	string thisExposure = maybeFromHeader(exposureName, localHeader);
	if (!exposureNames.has(thisExposure)) {
	  exposureNames.append(thisExposure);
	  exposurePointings.push_back(thisRADec);
	  exposureFields.push_back(fieldNumber);
	  exposureInstruments.push_back(instrumentNumber);
	  Assert(exposurePointings.size()==exposureNames.size());
	  Assert(exposurePointings.size()==exposureFields.size());
	  Assert(exposurePointings.size()==exposureInstruments.size());
	} 
	int exposureNumber = exposureNames.indexOf(thisExposure);
	// Error if different extensions put same exposure in different fields or instruments:
	if (fieldNumber != exposureFields[exposureNumber]) {
	  cerr << "Conflicting field assignments for exposure " << thisExposure << endl;
	  exit(1);
	}
	if (instrumentNumber != exposureInstruments[exposureNumber]) {
	  cerr << "Conflicting instrument assignment " << thisInstrument 
	       << " for exposure " << thisExposure << endl;
	  exit(1);
	}

	// Henceforth should use RA/Dec from the exposure.

	// Assign a device
	string thisDevice="";
	int deviceNumber = -1;
	if (instrumentNumber >= 0) {
	  thisDevice = maybeFromHeader(deviceName, localHeader);
	  if (!instruments[instrumentNumber].has(thisDevice)) 
	    instruments[instrumentNumber].newDevice(thisDevice);
	  deviceNumber = instruments[instrumentNumber].indexOf(thisDevice);
	  Assert(deviceNumber >= 0);
	}

	// Read in a WCS, from pixel coords to the tangent plane for this Field:
	astrometry::PixelMap* map = 0;
	if (!xyAreRaDec) {
	  if (wcsStream.is_open()) {
	    // Try the WCS file first
	    try {
	      map = new astrometry::SCAMPMap(wcsHead, &(fields[fieldNumber].orient));
	    } catch (astrometry::AstrometryError& m) {
	      map = 0;
	    }
	  } 
	  if (!map) {
	    // Try the local header if not
	    try {
	      map = new astrometry::SCAMPMap(localHeader, &(fields[fieldNumber].orient));
	    } catch (astrometry::AstrometryError& m) {
	      map = 0;
	    }
	  }
	  if (!map) {
	    // If still nothing, quit.   
	    cerr << "Could not generate WCS for HDU " << hduNumber
		 << " in file " << filename
		 << endl;
	    exit(1);
	  }
	}

	// Record information about this extension to an output table

	extensionTable.writeCell(filename, "Filename", extensionNumber);
	extensionTable.writeCell(iFile, "FileNumber", extensionNumber);
	extensionTable.writeCell(hduNumber, "HDUNumber", extensionNumber);
	extensionTable.writeCell(exposureNumber, "ExposureNumber", extensionNumber);
	extensionTable.writeCell(deviceNumber, "DeviceNumber", extensionNumber);
	extensionTable.writeCell(thisXKey, "xKey", extensionNumber);
	extensionTable.writeCell(thisYKey, "yKey", extensionNumber);
	extensionTable.writeCell(thisIdKey, "idKey", extensionNumber);
	extensionTable.writeCell(thisErrorKey, "errKey", extensionNumber);
	extensionTable.writeCell(weight, "Weight", extensionNumber);
	extensionTable.writeCell(thisWcsOutFile, "WCSOut", extensionNumber);

	{
	  string wcsDump;
	  if (!map) {
	    wcsDump = "ICRS";
	  } else {
	    const astrometry::SCAMPMap* sm = dynamic_cast<const astrometry::SCAMPMap*> (map);
	    if (!sm) {
	      cerr << "Only know how to serialize SCAMPMap now, file/extn "
		   << filename << " / " << iFile << endl;
	      exit(1);
	    }
	    img::Header hh = sm->writeHeader();
	    ostringstream oss;
	    oss << hh;
	    wcsDump = oss.str();
	  }
	  extensionTable.writeCell(wcsDump, "WCS", extensionNumber);
	}
	
	// Select the Catalog(s) to which points from this extension will be added
	PointCat* starCatalog = fields[fieldNumber].catalogFor(globalAffinity);
	PointCat* galaxyCatalog = (stringstuff::nocaseEqual(globalAffinity,
							    thisAffinity) ?
				   starCatalog : 
				   fields[fieldNumber].catalogFor(thisAffinity));
	
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
	    astrometry::SphericalICRS radec( vx[iObj]*DEGREE, vy[iObj]*DEGREE );
	    astrometry::TangentPlane tp(radec,  fields[fieldNumber].orient);
	    tp.getLonLat(xw, yw);
	    // Matching program will speak degrees, not radians:
	    xw /= DEGREE;
	    yw /= DEGREE;
	  } else {
	    Assert(map);
	    map->toWorld(xpix, ypix, xw, yw);
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

	if (map) delete map;
	
	extensionNumber++;
      } // end input extension loop
    } // end input file loop

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
	names.push_back(fields[i].name);
	astrometry::SphericalICRS pole = fields[i].orient.getPole();
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
      for (Field::CatMap::iterator iAffinity=fields[iField].catalogs.begin();
	   iAffinity != fields[iField].catalogs.end();
	   ++iAffinity) {
	const PointCat* pcat= iAffinity->second;
	cerr << "Catalog for field " << fields[iField].name
	     << " affinity " << iAffinity->first
	     << " with " << pcat->size() << " matches "
	     << endl;
	if (pcat->empty()) continue;
	FitsTable ft(outCatalogName, FITS::ReadWrite + FITS::Create, -1);
	ft.setName("MatchCatalog");
	FTable ff=ft.use();
	ff.header()->replace("Field", fields[iField].name, "Field name");
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

  } catch (std::runtime_error &m) {
    quit(m,1);
  }
}
