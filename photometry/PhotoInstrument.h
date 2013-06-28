#ifndef INSTRUMENT_H
#define INSTRUMENT_H
#include <map>
#include "Std.h"
#include "NameIndex.h"
#include "Astrometry.h"
#include "Wcs.h"
#include "PhotoMatch.h"
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
    Exposure(const string& name_, const astrometry::SphericalCoords& coords): 
      name(name_), projection(coords.duplicate()) {}
    ~Exposure() {delete projection;}
    string name;
    astrometry::SphericalCoords* projection;	// Projection relating world coords to sky for this exposure
    int  field;
    int  instrument;
    double airmass;	// Airmass at reference position
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
    Extension(): map(0), startWcs(0) {}
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

  // Class that represents an catalog of objects from a single device on single exposure
  // that will be used solely to extract color information.
  class ColorExtension {
  public:
    ColorExtension() {}
    int priority;	// Rank of this catalog in heirarchy of colors.  Lower value takes priority.
    // The objects from this catalog that we will use, and the matches they give colors for:
    std::map<long, Match*> keepers;  
  };

} // end namespace photometry

#endif
