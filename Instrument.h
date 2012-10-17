#ifndef INSTRUMENT_H
#define INSTRUMENT_H
#include "Std.h"
#include "Astrometry.h"
#include "Match.h"
#include <map>
#include "NameIndex.h"
#include "PixelMapCollection.h"

class Instrument {
public:
  Instrument(string name_):
    name(name_), nDevices(0) {}
  string name;
  int nDevices;
  NameIndex deviceNames;	// Names of all devices that exist for this instrument
  vector<int> exposures;	// Which exposures use this Instrument
  vector<PixelMapChain> maps;	// Instrument parts of PixelMaps for each device - key chains
  vector<PixelMap*> pixelMaps;	// Instrument parts of PixelMaps for each device - callable map
  // Keep track of range of pixel coords per device
  vector<Bounds<double> > domains;	// Rectangles bounding pixel coords of objects
  void addDevice(string devName, const Bounds<double>& devBounds=Bounds<double>()) {
    deviceNames.append(devName);
    maps.push_back(PixelMapChain());
    pixelMaps.push_back(0);
    domains.push_back(devBounds);
    nDevices = deviceNames.size();
    Assert(maps.size()==nDevices);
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
  Exposure(const string& name_, const Orientation& orient_): 
    name(name_), orient(orient_), reproject(-1) {}
  string name;
  Orientation orient;	// Telescope pointing for this one
  int  field;
  int  instrument;
  PixelMapKey reproject; // Map from exposure's TangentPlane to field's TangentPlane
  PixelMapChain warp;		// Exposure portion of PixelMap
private:
  // Hide:
  Exposure(const Exposure& rhs);
  void operator=(const Exposure& rhs);
};

// Class that represents an catalog of objects from a single device on single exposure.
// Will have originated from a single bintable HDU that we can access
class Extension {
public:
  Extension(): extensionMap(0),  startpm(0) {}
  int exposure;
  int device;
  SubMap* extensionMap; // Compounded maps for each extension (owned by PixelMapCollection)
  PixelMap* startpm;  // Input PixelMap for this extension (owned by this class)
  map<long, Detection*> keepers;  // The objects from this catalog that we will use
  ~Extension() {
    if (startpm) delete startpm;
  }
private:
};

#endif
