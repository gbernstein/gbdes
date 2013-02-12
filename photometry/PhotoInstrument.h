#ifndef INSTRUMENT_H
#define INSTRUMENT_H
#include "Std.h"
#include "Astrometry.h"
#include "Match.h"
#include <map>
#include "NameIndex.h"
#include "PixelMap.h"
#include "PhotoMapCollection.h"

namespace photometry {

  class Instrument {
  public:
    Instrument(string name_=""):
      name(name_), nDevices(0) {}
    string name;
    int nDevices;
    NameIndex deviceNames;	// Names of all devices that exist for this instrument
    vector<string> mapNames;	// Names of instrument PhotoMaps for each device
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
  private:
    // Hide:
    Instrument(const Instrument& rhs);
    void operator=(const Instrument& rhs);
  };

  // Class that represents a single pointing of the telescope:

  class Exposure {
  public:
    Exposure(const string& name_, const SphericalCoords& coords): 
      name(name_), projection(coords.duplicate()) {}
    ~Exposure() {delete projection;}
    string name;
    SphericalCoords* projection;	// Projection relating world coords to sky for this exposure
    int  field;
    int  instrument;
    string mapName;	// name of PhotoMap for this Exposure
  private:
    // Hide:
    Exposure(const Exposure& rhs);
    void operator=(const Exposure& rhs);
  };

  // Class that represents an catalog of objects from a single device on single exposure.
  // Will have originated from a single bintable HDU that we can access
  class Extension {
  public:
    Extension(): map(0), wcs(0), startWcs(0) {}
    int exposure;
    int device;
    SubMap* map;	  // The total magnitude transformation for this exposure
    astrometry::Wcs* startWcs;  // Input Wcs for this extension (owned by this class)
    std::map<long, Detection*> keepers;  // The objects from this catalog that we will use
    ~Extension() {
      if (startWcs) delete startWcs;
    }
  private:
  };

} // end namespace photometry

#endif
