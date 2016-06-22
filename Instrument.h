#ifndef INSTRUMENT_H
#define INSTRUMENT_H
#include "Std.h"
#include "Astrometry.h"
#include "Match.h"
#include <map>
#include "NameIndex.h"
#include "PixelMapCollection.h"

// Definitions are: "pixel coords" -> "instrument coords" via Instrument map
// instrument coords -> "world coords" via Exposure map (world coords are in a projection nominally about
//             the exposure's pointing)
// world coords -> ICRS coords via the Exposure's projection
// "field coords" -> ICRS coords via the Field's projection (field coords are sky coords in a
//              projection that is common to all Exposures in a field).
// world coords -> field coords via the Exposure's reprojection map 

class Instrument {
public:
  Instrument(string name_=""):
    name(name_), nDevices(0) {}
  string name;
  int nDevices;
  NameIndex deviceNames;	// Names of all devices that exist for this instrument
  vector<string> mapNames;	// Names of instrument PixelMaps for each device
  // Keep track of range of pixel coords per device
  vector<Bounds<double> > domains;	// Rectangles bounding pixel coords of objects
  void addDevice(string devName, const Bounds<double>& devBounds=Bounds<double>()) {
    deviceNames.append(devName);
    domains.push_back(devBounds);
    mapNames.push_back("");
    nDevices = deviceNames.size();
    Assert(mapNames.size()==nDevices);
    Assert(domains.size()==nDevices);
  }
  // Everything cleaned up in destructor:
  ~Instrument() {}
  // No copying:
  Instrument(const Instrument& rhs) =delete;
  void operator=(const Instrument& rhs) =delete;
};

// Class that represents a single pointing of the telescope:

class Exposure {
public:
  Exposure(const string& name_, const astrometry::SphericalCoords& coords): 
    name(name_), projection(coords.duplicate()) {}
  ~Exposure() {delete projection;}
  string name;   // The exposure map will have this name too.
  astrometry::SphericalCoords* projection;	// Projection relating world coords to sky for exposure
  int  field;
  int  instrument;
  // No copying:
  Exposure(const Exposure& rhs) =delete;
  void operator=(const Exposure& rhs) =delete;
};

// Class that represents an catalog of objects from a single device on single exposure.
// Will have originated from a single bintable HDU that we can access
class Extension {
public:
  Extension(): map(0), wcs(0), startWcs(0) {}
  int exposure;
  int device;
  string wcsName;      // Name of final WCS (and map into field coordinates)
  string basemapName;  // Name of map from pixel coords to exposure's projected coords
  astrometry::SubMap* map;	  // The total map from pixel coordinates to field coordinates.
  astrometry::Wcs* wcs;       // Wcs from pixel coordinates to sky coordinates = basemap + field projection
  astrometry::Wcs* startWcs;  // Input Wcs for this extension (owned by this class)
  std::map<long, astrometry::Detection*> keepers;  // The objects from this catalog that we will use
  ~Extension() {
    if (startWcs) delete startWcs;
  }
private:
};

#endif
