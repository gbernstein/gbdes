#ifndef WCSFOF_MATCH_H
#define WCSFOF_MATCH_H

#include "FoF.h"
#include "Astrometry.h"
#include "FitsTable.h"
#include "Units.h"
#include "FitSubroutines.h"

#include <vector>

/// Struct that will hold the info about each point that matcher (and subsequent programs)
/// will need
struct Point {
    Point(double x_, double y_, long ext, long obj, long expo)
            : x(2), extensionNumber(ext), objectNumber(obj), exposureNumber(expo) {
        x[0] = x_;
        x[1] = y_;
    }
    std::vector<double> x;
    long extensionNumber;  // "extension" is an individual input FITS bintable.
    long objectNumber;     // gives the object number of this Point in its input catalog
    long exposureNumber;   // Which exposure the catalog is from.
    const std::vector<double> getX() const { return x; }
};

// *** Here is the typedef that says we'll be collecting & matching Points
// *** using the Catalog template class in FoF.h.  Change this typedef if
// *** a better 2d matching structure is desired!
typedef fof::Catalog<Point, 2> PointCat;

struct Field {
    Field(const Field &) = delete;
    Field &operator=(Field const &) = delete;
    Field(Field &&) = delete;
    Field &operator=(Field &&) = delete;

    std::string name;
    // Coordinate system that has lat=lon=0 at field center and does projection to use for matching
    std::unique_ptr<astrometry::SphericalCoords> projection;
    double extent;
    double matchRadius;
    // Map from affinity name to its catalog of matches:
    typedef std::map<std::string, std::unique_ptr<PointCat>> CatMap;
    CatMap catalogs;
    Field() : projection() {}

    PointCat &catalogFor(const std::string &affinity) {
        if (catalogs.find(affinity) == catalogs.end()) {
            std::vector<double> lower(2, -extent);
            std::vector<double> upper(2, extent);
            // *** This line creates new point catalogs.  Might need to be changed if
            // *** we decide to use a different catalog structure.
            catalogs.insert(std::pair<std::string, std::unique_ptr<PointCat>>(
                    affinity, new PointCat(lower, upper, matchRadius)));
        }
        return *catalogs[affinity];
    }
};

// Right now a device is just a region of pixel coordinates, plus a name
struct Device : public Bounds<double> {
    std::string name;
};

// Instr is a collection of Devices, with a name
struct Instr : public std::vector<Device> {
    Instr(const std::string &name_ = "") : name(name_) {}
    std::string name;
    Instr(const img::FTable &ft) {
        ft.header()->getValue("Name", name);
        for (int i = 0; i < ft.nrows(); i++) {
            Device d;
            ft.readCell(d.name, "Name", i);
            double xmin, xmax, ymin, ymax;
            ft.readCell(xmin, "XMin", i);
            ft.readCell(xmax, "XMax", i);
            ft.readCell(ymin, "YMin", i);
            ft.readCell(ymax, "YMax", i);
            if (xmin != 0. || xmax != 0. || ymin != 0. || ymax != 0.) {
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
    std::string name;
    astrometry::SphericalICRS pointing;
    int field;
    int instrument;
    void read(const img::FTable &ft, int irow) {
        double ra, dec;
        ft.readCell(name, "Name", irow);
        ft.readCell(ra, "RA", irow);
        ft.readCell(dec, "Dec", irow);
        ft.readCell(field, "FieldNumber", irow);
        ft.readCell(instrument, "InstrumentNumber", irow);
        pointing.setRADec(ra * astrometry::WCS_UNIT, dec * astrometry::WCS_UNIT);
    }
};

class FoFClass {
public:
    FoFClass();
    FoFClass(Fields &fields_, std::vector<std::shared_ptr<Instrument>> instruments_,
             ExposuresHelper exposures_, std::vector<double> fieldExtents, double matchRadius);

    std::vector<std::unique_ptr<Field>> fields;
    std::vector<std::unique_ptr<Instrument>> instruments;
    std::vector<std::unique_ptr<Exposure>> exposures;
    img::FTable extensionTable;
    std::list<Point> allPoints;
    std::string currentPixelMapCollectionFileName;
    std::vector<int> sequence;
    std::vector<long> extn;
    std::vector<long> obj;

    void getExtensionInfo(long iextn, std::string &thisAffinity, int &exposureNumber, int &instrumentNumber,
                          int &fieldNumber, int &deviceNumber, std::vector<bool> &isStar,
                          std::vector<double> &vx, std::vector<double> &vy, std::vector<long> &vid);
    std::unique_ptr<astrometry::Wcs> getWCS(long iextn, int fieldNumber);
    void reprojectWCS(astrometry::Wcs &wcs, int fieldNumber);
    void addCatalog(const astrometry::Wcs &wcs, std::string thisAffinity, int exposureNumber, int fieldNumber,
                    int instrumentNumber, int deviceNumber, long iextn, std::vector<bool> isStar,
                    std::vector<double> vx, std::vector<double> vy, std::vector<long> vid);
    void writeMatches(std::string outCatalogName, int minMatches = 2, bool allowSelfMatches = false);
    void sortMatches(int fieldNumber, int minMatches = 2, bool allowSelfMatches = false);
};
#endif
