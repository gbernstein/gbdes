// Define all the PixelMapCollection routines.

#include "PixelMapCollection.h"
#include "PolyMap.h"

using namespace astrometry;

//////////////////////////////////////////////////////////////
// Static data
//////////////////////////////////////////////////////////////
vector<string>
PixelMapCollection::mapTypeNames;

vector<PixelMapCollection::Creator>
PixelMapCollection::mapCreators;

bool
PixelMapCollection::creatorsInitialized=false;

template <class MapType>
void
PixelMapCollection::registerMapType() {
  mapTypeNames.push_back(MapType::mapType());
  mapCreators.push_back(MapType::create);
}

void
PixelMapCollection::PixelMapTypeInitialize() {
  if (creatorsInitialized) return;
  creatorsInitialized = true;
  registerMapType<IdentityMap>();
  registerMapType<ReprojectionMap>();
  registerMapType<PolyMap>();
  registerMapType<LinearMap>();
}

//////////////////////////////////////////////////////////////
// Construction / destruction
//////////////////////////////////////////////////////////////

PixelMapCollection::PixelMapCollection(): parameterCount(0) {
  PixelMapTypeInitialize();
}

PixelMapCollection::~PixelMapCollection() {
  // Destroy all of the compounded maps and submaps we've made:
  for (MapIter i = mapElements.begin(); i!=mapElements.end(); ++i) {
    if (i->second.atom) {
      delete i->second.atom;
      i->second.atom = 0;
    }
    if (i->second.realization) {
      delete i->second.realization;
      i->second.realization = 0;
    }
  }
  for (WcsIter i = wcsElements.begin(); i!=wcsElements.end(); ++i) {
    if (i->second.nativeCoords) {
      delete i->second.nativeCoords;
      i->second.nativeCoords = 0;
    }
    if (i->second.realization) {
      delete i->second.realization;
      i->second.realization = 0;
    }
  }
}

//////////////////////////////////////////////////////////////
// Parameter manipluation
//////////////////////////////////////////////////////////////

void
PixelMapCollection::setParams(const DVector& p) {
  Assert(p.size()==parameterCount);
  for (MapIter i = mapElements.begin(); i!=mapElements.end(); ++i) {
    MapElement& map = i->second;
    if (map.isFixed) continue;
    int nSub = map.nParams;
    if (nSub<=0) continue;
    DVector subp(nSub, 0.);
    subp = p.subVector(map.startIndex, map.startIndex+nSub);
    Assert(map.atom)
    map.atom->setParams(subp);
  }
}

DVector
PixelMapCollection::getParams() const {
  DVector p(parameterCount, 888.);
  for (ConstMapIter i = mapElements.begin(); i!=mapElements.end(); ++i) {
    const MapElement& map = i->second;
    if (map.isFixed) continue;
    int nSub = map.nParams;
    if (nSub<=0) continue;
    Assert(map.atom);
    p.subVector(map.startIndex, map.startIndex+nSub) = 
      map.atom->getParams().subVector(0,nSub);
  }
  return p;
}

void 
PixelMapCollection::setFixed(list<string> nameList, bool isFixed) {
  set<string> fixThese;
  for (list<string>::const_iterator iName = nameList.begin();
       iName != nameList.end();
       ++iName) {
    set<string> addThese = dependencies(*iName);
    fixThese.insert(addThese.begin(), addThese.end());
  }
  for (set<string>::iterator iName=fixThese.begin();
       iName != fixThese.end();
       ++iName)
    mapElements[*iName].isFixed = isFixed;

  // Now calculate new parameter indices and propagate new configuration to all SubMaps
  rebuildParameterVector();
}

bool
PixelMapCollection::getFixed(string name) const {
  ConstMapIter i = mapElements.find(name);
  if (i==mapElements.end())
    throw AstrometryError("No member of PixelMapCollection named " + name);
  return i->second.isFixed;
}

void
PixelMapCollection::rebuildParameterVector() {
  // First assign new starting indices to every atomic map element

  // Restart the parameter counting:
  parameterCount = 0;
  for (MapIter i = mapElements.begin(); i!=mapElements.end(); ++i) {
    MapElement& map = i->second;
    // Only atomic map components go into the big parameter vector
    if (!(map.atom)) continue;
    int nMap = map.atom->nParams();
    if (map.isFixed || nMap==0) {
      // No free parameters here.
      map.nParams = 0;
    } else {
      // Have some parameters; append them to master list.
      map.nParams = nMap;
      map.startIndex = parameterCount;
      parameterCount += nMap;
    }
  }

  // Then go through all extant SubMaps and update their versions of vectors & parameter counts.
  for (MapIter i = mapElements.begin(); i!=mapElements.end(); ++i) {
    MapElement& map = i->second;
    // Is there a SubMap here that we need to update?
    SubMap* sub = map.realization;
    if (!sub) continue;

    // Loop through all the map elements that this SubMap uses
    for (int iElement = 0; iElement < sub->vMaps.size(); iElement++) {
      string elementName = sub->vMaps[iElement]->getName();
      ConstMapIter j = mapElements.find(elementName);
      if (j==mapElements.end()) 
	throw AstrometryError("PixelMapCollection::rebuildParameterVector could not find"
			      " map element with name " + elementName);
      if (!j->second.atom)
	throw AstrometryError("PixelMapCollection element " + elementName + " is not atomic"
			      " in rebuildParameterVector.");
      sub->vStartIndices[iElement] = j->second.startIndex;
      sub->vNSubParams[iElement] = j->second.nParams;
    }
    sub->countFreeParameters();
  }
}

set<string>
PixelMapCollection::dependencies(string target,
				 const set<string>& ancestors) const {
  // If target is on the ancestor list, there is a circular dependence
  if (ancestors.count(target) == 0)
    throw AstrometryError("Circular dependence in PixelMapCollection at element " + target);

  // Make a new ancestor list with target on it
  // (I'm sure there is a more clever algorithm for finding circular dependence!)
  set<string> ancestorsPlusTarget(ancestors);
  ancestorsPlusTarget.insert(target);

  // Find target MapElement and add target to output dependency list
  ConstMapIter iTarget = mapElements.find(target);
  if (iTarget == mapElements.end()) 
    throw AstrometryError("PixelMapCollection has no element " + target + " in dependencies()");
  set<string> output;
  output.insert(target);

  // Call this routine recursively for all the map elements that the target depends on
  for (list<string>::const_iterator i = iTarget->second.subordinateMaps.begin();
       i != iTarget->second.subordinateMaps.begin();
       ++i) {
    set<string> subs = dependencies(*i, ancestorsPlusTarget);
    output.insert(subs.begin(), subs.end());
  }

  return output;
}

//////////////////////////////////////////////////////////////
// Member maintenance
//////////////////////////////////////////////////////////////

vector<string>
PixelMapCollection::allMapNames() const {
  vector<string> output;
  for (ConstMapIter i = mapElements.begin();
       i != mapElements.end();
       ++i) 
    output.push_back(i->first);
  return output;
}

vector<string>
PixelMapCollection::allWcsNames() const {
  vector<string> output;
  for (ConstWcsIter i = wcsElements.begin();
       i != wcsElements.end();
       ++i) 
    output.push_back(i->first);
  return output;
}

void 
PixelMapCollection::addMap(PixelMap* pm, bool duplicateNamesAreExceptions) {
  // ???
  // ??? need to distinguish compound, submaps, and wcs's from atomic types.
}

void 
PixelMapCollection::addWcs(Wcs* pm, bool duplicateNamesAreExceptions) {
  // ???
}

// Define a new pixelMap that is compounding of a list of other PixelMaps.  Order
// of the list is from pixel to world coordinates.
void 
PixelMapCollection::defineChain(string chainName, const list<string>& elements) {
  if (mapExists(chainName)) 
    throw AstrometryError("PixelMapCollection::defineChain with duplicate name: " 
			  + chainName);
  // Check that elements exist
  for (list<string>::const_iterator i = elements.begin();
       i != elements.end();
       ++i) 
    if (!mapExists(*i))
      throw AstrometryError("PixelMapCorrection::defineChain with unknown pixel map element: "
			    + *i);
  mapElements.insert(std::pair<string, MapElement>(chainName, MapElement()));
  mapElements[chainName].subordinateMaps = elements;
  // Note that adding a new chain does not change parameter vector assignments
}

// Define a WCS system to be the given PixelMap followed by projection to the
// sky described by the nativeCoords system.
void 
PixelMapCollection::defineWcs(string wcsName, const SphericalCoords& nativeCoords, 
			      string mapName, double wScale) {
  // Check that Wcs name does not exist
  if (mapExists(wcsName)) 
    throw AstrometryError("PixelMapCollection::defineWcs with duplicate name: " 
			  + wcsName);
  // and that the PixelMap does:
  if (!mapExists(mapName)) 
      throw AstrometryError("PixelMapCorrection::defineWcs with undefined pixelMap element: "
			    + mapName);
  wcsElements.insert(std::pair<string, WcsElement>(wcsName, WcsElement()));
  wcsElements[wcsName].mapName = mapName;
  wcsElements[wcsName].nativeCoords = nativeCoords.duplicate();
  wcsElements[wcsName].wScale = wScale;
}

// Return pointer to a SubMap realizing the named coordinate transformation
SubMap* 
PixelMapCollection::issueMap(string mapName) {
  if (!wcsExists(mapName))
    throw AstrometryError("PixelMapCollection::issueMap requested for unknown PixelMap: "
			  +  mapName);
  MapElement& el = mapElements[mapName];

  if (!el.realization) {
    // Create a realization if one does not exist
    list<string> atomList = orderAtoms(mapName);
    vector<int> startIndices(atomList.size(), 0);
    vector<int> nParams(atomList.size(), 0);

    list<PixelMap*> atoms;
    int index=0;
    for (list<string>::const_iterator i = atomList.begin();
	 i != atomList.end();
	 ++i, ++index) {
      Assert(mapElements[*i].atom);	// All elements should be atomic
      atoms.push_back(mapElements[*i].atom);
      // fill in its indices into master vector:
      startIndices[index] = mapElements[*i].startIndex;
      nParams[index] = mapElements[*i].isFixed ? 0 : mapElements[*i].nParams;
    }
    SubMap* sm = new SubMap(atoms, mapName);
    sm->vStartIndices = startIndices;
    sm->vNSubParams = nParams;
    sm->countFreeParameters();
    el.realization = sm;
  }
  return el.realization;
}

// Return a pointer to a Wcs built from a SubMap and realizing the named coord system
Wcs* 
PixelMapCollection::issueWcs(string wcsName) {
  if (!wcsExists(wcsName))
    throw AstrometryError("PixelMapCollection::issueWcs requested for unknown Wcs: "
			  +  wcsName);
  WcsElement& el = wcsElements[wcsName];

  if (!el.realization) {
    // Create a realization if one does not exist
    SubMap* sm = issueMap(el.mapName);
    el.realization = new Wcs(sm, *el.nativeCoords, wcsName, el.wScale);
  }
  return el.realization;
}

list<string> 
PixelMapCollection::orderAtoms(string mapName,
			       const set<string>& ancestors) const {
  // ???
}

//////////////////////////////////////////////////////////////
// (De-) Serialization
//////////////////////////////////////////////////////////////

void 
PixelMapCollection::read(istream& is, string namePrefix) {
  // ???
}
void 
PixelMapCollection::write(ostream& os) const {
  // ???
}
void 
PixelMapCollection::writeMap(ostream& os, string name) const {
  // ???
}
void 
PixelMapCollection::writeWcs(ostream& os, string name) const {
  // ???
}
