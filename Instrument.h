#ifndef INSTRUMENT_H
#define INSTRUMENT_H
#include "Std.h"
#include "Astrometry.h"
#include "Match.h"
#include <map>
#include "NameIndex.h"
#include "PixelMapCollection.h"

class Exposure;

class Instrument {
public:
  Instrument(string name_):
    name(name_), nDevices(0) {}
  string name;
  int nDevices;
  NameIndex extensionNames;	// Names of all extensions that exist for this instrument
  vector<int> exposures;	// Which exposures use this Instrument
  vector<PixelMapChain> maps;	// Instrument parts of PixelMaps for each extension - key chains
  vector<PixelMap*> pixelMaps;	// Instrument parts of PixelMaps for each extension - callable map
  // Keep track of range of pixel coords per extension
  vector<Bounds<double> > domains;	// Rectangles bounding pixel coords of objects
  void addDevice(string extName) {
    extensionNames.append(extName);
    maps.push_back(PixelMapChain());
    pixelMaps.push_back(0);
    domains.push_back(Bounds<double>());
    nDevices = extensionNames.size();
    Assert(maps.size()==nDevices);
    Assert(domains.size()==nDevices);
  }
  // Everything cleaned up in destructor:
  ~Instrument() {}
private:
  // Hide:
  Instrument(const Instrument& rhs) {}
  void operator=(const Instrument& rhs) {}
};

class Exposure {
public:
  Exposure(): orient(0), reproject(-1) {}
  string name;
  int field;
  int instrument;
  Orientation *orient;	// Telescope pointing for this one
  PixelMapKey reproject; // Map from exposure's TangentPlane to field's TangentPlane
  PixelMapChain warp;		// Exposure portion of PixelMap
  vector<SubMap*> extensionMaps; // Compounded maps for each extension (owned by PixelMapCollection)
  vector<PixelMap*> startpm;  // Input PixelMap for this extension (owned by this class)
  // Everything cleaned up in destructor:
  ~Exposure() {
    if (orient) delete orient;
    for (int i=0; i<startpm.size(); i++)
      if (startpm[i]) delete startpm[i];
  }
private:
  // Hide:
  Exposure(const Exposure& rhs) {}
  void operator=(const Exposure& rhs) {}
};

#endif
