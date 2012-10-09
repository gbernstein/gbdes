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

// Parse the inValue of string: if it starts with @ sign, specifies to read a value
// from keyword of a header:
string maybeFromHeader(const string& inValue, const Header& h) {
  if (inValue.empty() || inValue[0]!='@') return inValue;
  string result;
  if (!h.getValue(inValue.substr(1), result)) 
    throw std::runtime_error("Could not find string-valued header keyword " + 
			     inValue.substr(1));
  return result;
}

// Struct that will hold the info about each point that matcher (and subsequent progs)
// will need
struct Point {
  Point(double x_, double y_, long ext, long obj):
    x(2), extensionNumber(ext), objectNumber(obj) {x[0]=x_; x[1]=y_;}
  vector<double> x;
  long extensionNumber;
  long objectNumber;
  const vector<double> getX() const {return x;}
};


typedef fof::Catalog<Point,2> PointCat;

struct Field {
  string name;
  astrometry::Orientation orient;
  double extent;
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
								 extent)));
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

  // Convert matching radius to our system units (radians)
  matchRadius *= ARCSEC;

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
	  f.extent = extent * DEGREE;
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
    const string xkeyCol="XKEY";
    const string ykeyCol="YKEY";
    const string idkeyCol="IDKEY";

    const string globalAffinity="GLOBAL";

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
      fileTable.readCell(idkey, idkeyCol, iFile);
      // The affinity and selections we can take defaults for:
      string selectionExpression;
      try {
	fileTable.readCell(selectionExpression,inSelectCol, iFile);
      } catch (FTableNonExistentColumn& m) {
	selectionExpression.clear();
      }
      string affinityName;
      try {
	fileTable.readCell(affinityName, affinityCol, iFile);
      } catch (FTableNonExistentColumn& m) {
	affinityName = globalAffinity;
      }
      string starExpression;
      try {
	fileTable.readCell(starExpression,starSelectCol, iFile);
      } catch (FTableNonExistentColumn& m) {
	starExpression.clear();	// everything becomes a star
      }

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
      bool xyAreRaDec = stringstuff::nocaseEqual(wcsFile, "ICRS")
	|| stringstuff::nocaseEqual(wcsFile, "@_ICRS");
      ifstream wcsStream;
      if (! (wcsFile.empty() || xyAreRaDec)) {
	wcsStream.open(wcsFile.c_str());
	if (!wcsStream.is_open()) {
	  cerr << "Could not open WCS information file " << wcsFile << endl;
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

	// Go ahead and filter the input rows and get the columns we want
	if (!selectionExpression.empty())
	  ft.filterRows(selectionExpression);
	vector<double> vx;
	ft.readCells(vx, thisXKey);
	vector<double> vy;
	ft.readCells(vy, thisYKey);
	vector<long> vid;
	ft.readCells(vid, thisIdKey);
	vector<bool> isStar(ft.nrows(), true);
	if (!starExpression.empty())
	  ft.evaluate(isStar, starExpression);

	// Determine center coordinates of this extension:
	// ??? only need these if finding nearest field?
	string thisRA = maybeFromHeader(ra, localHeader);
	string thisDec = maybeFromHeader(dec, localHeader);
	astrometry::SphericalICRS thisRADec;
	{
	  string coords = thisRA + " " + thisDec;
	  istringstream iss(coords);
	  if (!(iss >> thisRADec)) {
	    cerr << "Error reading RA & Dec for file " << filename
		 << " from string <" << coords << ">" << endl;
	    exit(1);
	  }
	}

	// Assign an exposure number:
	string thisExposure = maybeFromHeader(exposureName, localHeader);
	if (!exposureNames.has(thisExposure)) {
	  exposureNames.append(thisExposure);
	  exposurePointings.push_back(thisRADec);
	  Assert(exposurePointings.size()==exposureNames.size());
	}
	int exposureNumber = exposureNames.indexOf(thisExposure);
	// Reposition all devices of an exposure to same RA/Dec ???
	thisRADec = exposurePointings[exposureNumber];

	// And assign a field:
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

	// Assign a device
	string thisDevice = maybeFromHeader(deviceName, localHeader);
	if (!instruments[instrumentNumber].has(thisDevice)) 
	  instruments[instrumentNumber].newDevice(thisDevice);
	int deviceNumber = instruments[instrumentNumber].indexOf(thisDevice);
	Assert(deviceNumber >= 0);

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

	// ??? Record information about this extension to an output catalog???
	long extensionNumber = 0; // ???
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
	  allPoints.push_back(Point(xw, yw, extensionNumber, vid[iObj]));
	  if (isStar[iObj])
	    starCatalog->add(allPoints.back());
	  else
	    galaxyCatalog->add(allPoints.back());
	} // end object loop

	if (map) delete map;
	
      } // end input extension loop
    } // end input file loop

    //  Open output catalog
    //  Put input file table into output table?
    //  Write instrument information to output table?
    //  Loop over fields
    //  .. loop over affinities per field
    //  ....purge <minMatches and self-matches
    //  ....write output extension, with field center and affinity in header
    //  ....write match catalog: sequence, input file, input extension, input ID


  } catch (std::runtime_error &m) {
    quit(m,1);
  }
}
