// Program to match catalogs based on world coordinates.

#include "Std.h"
#include "FoF.h"
#include "FitsTable.h"
#include "Astrometry.h"
#include "NameIndex.h"
#include "Pset.h"
#include <map>
#include <iostream>
#include "PixelMapCollection.h"
#include "TPVMap.h"
#include "TemplateMap.h"
#include "StringStuff.h"
#include "Units.h"

#include "FitSubroutines.h"

using namespace std;
using namespace img;
using namespace FITS;
using stringstuff::stripWhite;

using astrometry::WCS_UNIT;  // we will want all WCS's to operate in same units

string usage=
  "Match objects using FoF algorithm on world coordinate system\n"
  "WCSFoF <configfile> [parameter file] [parameter file...]\n"
  "   [-parameter[=]value...]\n"
  "      <configfile>: FITS file with binary tables describing input catalogs\n"
  "      Program parameters specified as command-line options or read from\n"
  "          parameter file(s) specified on cmd line";
  
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

// Instr is a collection of Devices, with a name
struct Instr: public vector<Device> {
  Instr(const string& name_=""): name(name_) {}
  string name;
  Instr(const FTable& ft) {
    ft.header()->getValue("Name",name);
    for (int i=0; i<ft.nrows(); i++) {
      Device d;
      ft.readCell(d.name,"Name",i);
      double xmin, xmax, ymin, ymax;
      ft.readCell(xmin,"XMin",i);
      ft.readCell(xmax,"XMax",i);
      ft.readCell(ymin,"YMin",i);
      ft.readCell(ymax,"YMax",i);
      if (xmin!=0. || xmax!=0. || ymin!=0. || ymax!=0.) {
	// initialize bounds if there are some, otherwise leave undefined
	d.setXMin(xmin);
	d.setXMax(xmax);
	d.setYMin(ymin);
	d.setYMax(ymax);
      }
      push_back(d);
    }
  }
};

struct Expo {
public:
  string name;
  astrometry::SphericalICRS pointing;
  int field;
  int instrument;
  void read(const FTable& ft, int irow) {
    double ra, dec;
    ft.readCell(name, "Name", irow);
    ft.readCell(ra, "RA", irow);
    ft.readCell(dec, "Dec", irow);
    ft.readCell(field, "FieldNumber", irow);
    ft.readCell(instrument, "InstrumentNumber", irow);
    pointing.setRADec(ra*WCS_UNIT, dec*WCS_UNIT);
  }
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
			 "filename for FITS output catalog", "match.fof");
    parameters.addMember("minMatch",&minMatches, def | low,
			 "Minimum number of detections for usable match", 2, 2);
    parameters.addMember("selfMatch",&allowSelfMatches, def,
			 "Retain matches that have 2 elements from same exposure?", false);
  }

  ////////////////////////////////////////////////
  // Read parameters
  ////////////////////////////////////////////////

  // Read all the command-line and parameter-file program parameters
  processParameters(parameters, usage, 1, argc, argv);
  string configFile = argv[1];

  // Convert matching radius to our system units for world coords
  matchRadius *= ARCSEC/WCS_UNIT;

  try {
    // Teach PixelMapCollection about all types of PixelMaps it might need to deserialize
    loadPixelMapParser();

    ////////////////////////////////////////////////
    // Read in field table and save away orientations of tangent plane systems for each
    ////////////////////////////////////////////////
    FTable fieldTable;
    try {
      fieldTable = FitsTable(configFile, FITS::ReadOnly, "Fields").extract();
    } catch (std::runtime_error& e) {
      quit(e,1);
    }

    vector<Field*> fields;

    for (int ifield = 0; ifield < fieldTable.nrows(); ifield++) {
      string name;
      double ra;
      double dec;
      double extent;
      fieldTable.readCell(name, "NAME", ifield);
      fieldTable.readCell(ra, "RA", ifield);
      fieldTable.readCell(dec, "DEC", ifield);
      fieldTable.readCell(extent, "RADIUS", ifield);

      Field* f = new Field;
      f->name = name;
      astrometry::Orientation orient(astrometry::SphericalICRS(ra*WCS_UNIT, dec*WCS_UNIT));
      f->projection = new astrometry::Gnomonic(orient);
      f->extent = extent;
      f->matchRadius = matchRadius;
      fields.push_back(f);
    } // Done reading fields
    Assert(!fields.empty());

    // Now read in all the instrument tables
    list<int> instrumentHDUs;
    {
      // Find which extensions are instrument tables
      FITS::FitsFile ff(configFile);
      for (int i=1; i<ff.HDUCount(); i++) {
	FITS::Hdu h(configFile, FITS::HDUAny, i);
	if (stringstuff::nocaseEqual(h.getName(), "Instrument"))
	  instrumentHDUs.push_back(i);
      }
    }
    vector<FTable> instrumentTables(instrumentHDUs.size());
    vector<Instr*> instruments(instrumentHDUs.size(),0);

    for (list<int>::const_iterator i=instrumentHDUs.begin();
	 i != instrumentHDUs.end();
	 ++i) {
      FTable ft = FitsTable(configFile, FITS::ReadOnly, *i).extract();
      string instrumentName;
      int instrumentNumber;
      if (!ft.header()->getValue("Name", instrumentName)
	  || !ft.header()->getValue("Number", instrumentNumber)) {
	cerr << "Could not read name and/or number of instrument at extension "
	     << *i << endl;
      }
      // Get rid of spaces in the name
      spaceReplace(instrumentName);
      Assert(instrumentNumber < instrumentTables.size() && instrumentNumber>=0);
      ft.header()->setValue("Name",instrumentName);

      // Now save the table and make an Instrument structure
      instrumentTables[instrumentNumber] = ft;
      instruments[instrumentNumber] = new Instr(ft);
    }
    // Check that all Instruments were read
    for (int i=0; i<instruments.size(); i++)
      if (!instruments[i]) {
	cerr << "Failed to read instrument number " << i << endl;
	exit(1);
      }


    // Now get information on exposures
    FTable exposureTable;
    try {
      exposureTable = FitsTable(configFile, FITS::ReadOnly, "Exposures").extract();
    } catch (std::runtime_error& e) {
      quit(e,1);
    }
      
    vector<Expo*> exposures;
    for (int i=0; i<exposureTable.nrows(); i++) {
      Expo* e = new Expo;
      e->read(exposureTable, i);
      exposures.push_back(e);
    }
    Assert(exposures.size() == exposureTable.nrows());
	   
    // And a list holding all Points being matched.  Note that catalogs do not make copies
    // of the Points, they just hold pointers to them.  So this is a warehouse.
    // If we wanted to save the memory of the list links we could just let all the Points
    // be zombies that are destroyed when the program ends.
    list<Point> allPoints;

    // Table of information about every catalog file to read:
    FTable extensionTable;
    try {
      extensionTable = FitsTable(configFile, FITS::ReadOnly, "Extensions").extract();
    } catch (std::runtime_error& e) {
      quit(e,1);
    }

    // If WCS is coming from a serialized PixelMapCollection, I'll keep last-used one around
    // to avoid re-parsing it all the time:
    string currentPixelMapCollectionFileName = "";
    astrometry::PixelMapCollection* inputPmc=0;

    // Loop over input catalogs
    for (long iextn = 0; iextn < extensionTable.nrows(); iextn++) {
      string filename;
      extensionTable.readCell(filename, "FILENAME", iextn);
      stripWhite(filename);
      int hduNumber;
      extensionTable.readCell(hduNumber, "EXTENSION", iextn);

      if (iextn%10==0) cerr << "# Reading object catalog " << iextn
			    << "/" << extensionTable.nrows()
			    << " in " << filename 
			    << " HDU #" << hduNumber
			    << endl;
      FTable ft = FitsTable(filename, FITS::ReadOnly, hduNumber).extract();

      string idKey;
      extensionTable.readCell(idKey, "IDKEY",iextn);
      
      if (stringstuff::nocaseEqual(idKey, "_ROW")) {
	// Want the input table row number as ID column, so make such a column
	vector<long> vi(ft.nrows());
	for (int i=0; i<vi.size(); i++) vi[i] = i;
	ft.addColumn(vi, "_ROW");
      }

      // Get exposure number, instrument, field, and device numbers
      int exposureNumber;
      extensionTable.readCell(exposureNumber, "Exposure", iextn);
      int instrumentNumber = exposures[exposureNumber]->instrument;
      int fieldNumber = exposures[exposureNumber]->field;
      int deviceNumber;
      extensionTable.readCell(deviceNumber, "Device", iextn);

      // Now get the WCS for this extension
      astrometry::Wcs* wcs = 0;
      string wcsin;
      extensionTable.readCell(wcsin, "WCSIN", iextn);
      list<string> wcssplit = stringstuff::split(wcsin,'@');
      if (wcssplit.size()==2) {
	// @ sign signals that we are getting a map from a PMC
	string wcsName = wcssplit.front();
	string pmcFile = wcssplit.back();
	// Replace any white space in wcs name with underscore:
	wcsName = stringstuff::regexReplace("[[:space:]]+","_",wcsName);

	if (pmcFile!=currentPixelMapCollectionFileName) {
	  // Need to open this PMC file
	  // Try opening file as new serialized PMC
	  ifstream wcsStream(pmcFile.c_str());
	  if (!wcsStream.is_open()) {
	    cerr << "Could not open PixelMapCollection file <" << pmcFile 
		 << ">" << endl;
	    exit(1);
	  }
	  astrometry::PixelMapCollection* tryPMC = new astrometry::PixelMapCollection;
	  if (tryPMC->read(wcsStream)) {
	    // Successfully found serialized file:
	    if (inputPmc) delete inputPmc;
	    inputPmc = tryPMC;
	    currentPixelMapCollectionFileName = pmcFile;
	  } else {
	    cerr << "Failed deserializing PixelMapCollection file <" << pmcFile 
		 << ">" << endl;
	    exit(1);
	  }
	} // Done reading the pmcFile

	// Read the wcs from the collection
	wcs = inputPmc->issueWcs(wcsName)->duplicate();
      } else if (stringstuff::nocaseEqual(wcsin, "_ICRS")) {
	// No mapping will be needed.  Do nothing.
      } else {
	// Assume that wcsin is a FITS WCS specification
	Header wcsHead;
	//	wcsin += "\nEND\n";  //*** hack to fix problem in astropy.io.fits
	istringstream wcsStream(wcsin);
	if (!(wcsStream >> wcsHead)) {
	  cerr << "Error reading WCS/header from WCSIN at extension " << iextn << endl;
	  exit(1);
	}
	wcs = astrometry::readTPV(wcsHead);
	if (!wcs) {
	  cerr << "Failed reading TPV WCS from extension " << iextn << endl;
	  exit(1);
	}
      }

      // Now replace WCSIN in table with a new serialized one, 
      // and set projection to be the field coordinates
      if (wcs) {
	string wcsDump;
	astrometry::PixelMapCollection pmc;
	pmc.learnWcs(*wcs);
	ostringstream oss;
	pmc.writeWcs(oss,wcs->getName());
	wcsDump = oss.str();
	wcs->reprojectTo(*fields[fieldNumber]->projection);
	extensionTable.writeCell(wcsDump, "WCSIN", iextn);
      }

      // Now begin reading the input data

      // Filter the input rows and get the columns we want
      {
	string selectionExpression;
	extensionTable.readCell(selectionExpression, "SELECT", iextn);
	stripWhite(selectionExpression);
	if (!selectionExpression.empty())
	  ft.filterRows(selectionExpression);
      }
      vector<double> vx;
      string key;
      extensionTable.readCell(key, "XKEY", iextn);
      try {
	ft.readCells(vx, key);
      } catch (FTableError& m) {
	// Trap for using float column in source file instead of double:
	vector<float> vf;
	ft.readCells(vf, key);
	vx.resize(vf.size());
	for (int i=0; i<vf.size(); i++) vx[i]=vf[i];
      }
      vector<double> vy;
      extensionTable.readCell(key, "YKEY", iextn);
      try {
	ft.readCells(vy, key);
      } catch (FTableError& m) {
	// Trap for using float column in source file instead of double:
	vector<float> vf;
	ft.readCells(vf, key);
	vy.resize(vf.size());
	for (int i=0; i<vf.size(); i++) vy[i]=vf[i];
      }
      vector<long> vid;
      extensionTable.readCell(key, "IDKEY", iextn);
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
	string starExpression;
	extensionTable.readCell(starExpression, "STARSELECT", iextn);
	stripWhite(starExpression);
	if (!starExpression.empty())
	  ft.evaluate(isStar, starExpression);
      }

      // Select the Catalog(s) to which points from this extension will be added
      string thisAffinity;
      extensionTable.readCell(thisAffinity, "AFFINITY", iextn);

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
	  (*instruments[instrumentNumber])[deviceNumber]+=Position<double>(xpix,ypix);
	  
	//  maps coords to field's tangent plane
	double xw, yw;
	if (wcs) {
	  Assert(wcs);
	  wcs->toWorld(xpix, ypix, xw, yw);
	} else {
	  // Already have RA, Dec, just project them
	  // Assume they are in our std units (see Units.h)
	  fields[fieldNumber]->projection->convertFrom(astrometry::SphericalICRS(vx[iObj]*WCS_UNIT, 
										 vy[iObj]*WCS_UNIT ));
	  fields[fieldNumber]->projection->getLonLat(xw, yw);
	  // Matching program will use our std unit (Units.h), not radians:
	  xw /= WCS_UNIT;
	  yw /= WCS_UNIT;
	}

	allPoints.push_back(Point(xw, yw, iextn, vid[iObj], exposureNumber));

	// Now we match!!!!
	if (isStar[iObj])
	  starCatalog->add(allPoints.back());
	else
	  galaxyCatalog->add(allPoints.back());
      } // end object loop

      if (wcs) delete wcs;
    } // end input extension loop

    if (inputPmc) delete inputPmc;

    cerr << "*** Read " << allPoints.size() << " objects" << endl;

    //  Write all of our tables to output file
    {
      // Fields
      FitsTable ft(outCatalogName, FITS::ReadWrite + FITS::OverwriteFile, "Fields");
      ft.copy(fieldTable);
    }
    {
      // Exposures
      FitsTable ft(outCatalogName, FITS::ReadWrite + FITS::Create, "Exposures");
      ft.copy(exposureTable);
    }

    for (int i=0; i<instruments.size(); i++) {
      // Instrument tables
      // Update the device bounds in the table first
      const Instr& inst = *instruments[i];
      int nDevices = inst.size();
      vector<double> vxmin(nDevices);
      vector<double> vxmax(nDevices);
      vector<double> vymin(nDevices);
      vector<double> vymax(nDevices);
      for (int j=0; j<nDevices; j++) {
	vxmin[j] = floor(inst[j].getXMin());
	vxmax[j] = ceil(inst[j].getXMax());
	vymin[j] = floor(inst[j].getYMin());
	vymax[j] = ceil(inst[j].getYMax());
      }
      instrumentTables[i].writeCells(vxmin, "XMin");
      instrumentTables[i].writeCells(vxmax, "XMax");
      instrumentTables[i].writeCells(vymin, "YMin");
      instrumentTables[i].writeCells(vymax, "YMax");

      FitsTable ft(outCatalogName, FITS::ReadWrite + FITS::Create, -1);
      ft.setName("Instrument");
      ft.setVersion(i+1); // might use this later, or put into extension name???
      ft.copy(instrumentTables[i]);
    }

    {
      // Extension table
      FitsTable ft(outCatalogName, FITS::ReadWrite + FITS::Create, "Extensions");
      ft.copy(extensionTable);
    }

    // Now write all the match catalogs
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
    cerr << "Cleaning fields: " << endl;
    for (int i=0; i<fields.size(); i++)
      delete fields[i];
    cerr << "Cleaning instruments: " << endl;
    for (int i=0; i<instruments.size(); i++)
      delete instruments[i];
    cerr << "Cleaning exposures: " << endl;
    for (int i=0; i<exposures.size(); i++)
      delete exposures[i];
    cerr << "Done!: " << endl;

  } catch (std::runtime_error &m) {
    quit(m,1);
  }
}
