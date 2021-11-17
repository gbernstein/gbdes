#ifndef WCSFOF_MATCH_H
#define WCSFOF_MATCH_H

#include "FoF.h"
#include "Astrometry.h"
//#include "Bounds.h"
#include "FitsTable.h"
#include "Units.h"
#include "FitSubroutines.h"


#include <vector>

using namespace std;
using namespace img;

using astrometry::WCS_UNIT; 

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

class FoFClass {
  public:
    FoFClass();
    double matchRadius = 1.0;
    bool useAffinities = true;
    //string outCatalogName;
    int minMatches = 2;
    bool allowSelfMatches = false;
    vector<Field*> fields;
    vector<Instr*> instruments;
    vector<Expo*> exposures;
    FTable extensionTable;
    list<Point> allPoints;
    string currentPixelMapCollectionFileName = "";
    astrometry::PixelMapCollection* inputPmc=0;
    vector<int> sequence;
    vector<long> extn;
    vector<long> obj;
    int a = 5;
    int b = 4;
    int c;
    void addTest();


    void getExtensionInfo(long iextn, string thisAffinity, int exposureNumber, int instrumentNumber,
                          int fieldNumber, int deviceNumber, vector<bool> isStar, vector<double> vx, vector<double> vy, vector<long> vid);
    void getWCS(long iextn, int fieldNumber, astrometry::Wcs* wcs);
    void reprojectWCS(astrometry::Wcs* wcs, int fieldNumber);
    void addCatalog(astrometry::Wcs* wcs, string thisAffinity, int exposureNumber, int fieldNumber, int instrumentNumber,
                    int deviceNumber, long iextn, vector<bool> isStar, vector<double> vx, vector<double> vy, vector<long> vid);
    void writeMatches(string outCatalogName);
    void sortMatches(int fieldNumber);//, vector<int> sequence, vector<long> extn, vector<long> obj);
};

void friendOfFriend(double matchRadius, bool useAffinities, string outCatalogName, int minMatches,
                    bool allowSelfMatches, vector<Field*> fields, vector<Instr*> instruments,
                    vector<Expo*> exposures, FTable extensionTable);

#endif