#ifndef INSTRUMENT_H
#define INSTRUMENT_H
#include <map>
#include "Std.h"
#include "Bounds.h"
#include "Astrometry.h"
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
    name(name_), nDevices(0), band(name_) {}
  string name;
  string band;		// "band" could be more generic than instrument name.
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
  name(name_), projection(coords.duplicate()), mjd(0.), exptime(1.), airmass(1.),
  apcorr(0.) {}
  ~Exposure() {delete projection;}
  string name;   // The exposure map will have this name too.
  astrometry::SphericalCoords* projection;	// Projection relating world coords to sky for exposure
  int  field;
  int  instrument;
  double airmass;
  double exptime;
  double mjd;
  double apcorr;
  string epoch;
  // No copying:
  Exposure(const Exposure& rhs) =delete;
  void operator=(const Exposure& rhs) =delete;
};

// Class that represents an catalog of objects from a single device on single exposure.
// Will have originated from a single bintable HDU that we can access.
// The template argument are SubMap, Detection from either astrometry or photometry.
template <class T1, class T2>
class ExtensionBase {
public:
  ExtensionBase(): map(nullptr), wcs(nullptr), startWcs(nullptr), needsColor(false) {}
  int exposure;
  int device;
  double airmass;      // airmass and apcorr used for nightly priors
  double apcorr;
  double magshift;	// Additive adjustment to all incoming mags (exposure time)

  string wcsName;      // Name of final WCS (and map into field coordinates)
  string mapName;      // Name of photometry map or astrometric map into exposure coords
  T1* map;	       // The map from pixel coordinates to field coordinates.
  astrometry::Wcs* wcs;       // Wcs from pixel coordinates to sky coordinates = basemap + field projection
  astrometry::Wcs* startWcs;  // Input Wcs for this extension (owned by this class)
  bool needsColor;	// Save info on whether map requires color information.

  std::map<long, T2*> keepers; // The objects from this Extension catalog that we will use
  ~ExtensionBase() {
    if (startWcs) delete startWcs;
  }
};

// Class that represents an catalog of objects from a single device on single exposure
// that will be used solely to extract color information.
// The template argument is Match from either astrometry or photometry
template <class T>
class ColorExtensionBase {
public:
  ColorExtensionBase() {}
  int priority;	// Rank of this catalog in heirarchy of colors.  Lower value takes priority.
  // The objects from this catalog that we will use, and the matches they give colors for:
  std::map<long, T*> keepers;  
};

#endif
