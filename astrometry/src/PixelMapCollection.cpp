
// Define all the PixelMapCollection routines.

#include <set>
#include "PixelMapCollection.h"
#include "PolyMap.h"
#include "StringStuff.h"
#include <typeinfo>

using namespace astrometry;

#ifdef USE_YAML  // Don't build any of this without YAML

// The word that we expect all serialized collections to have on their first line
const string
PixelMapCollection::magicKey = "PixelMapCollection";

//////////////////////////////////////////////////////////////
// Static data
//////////////////////////////////////////////////////////////
map<string,PixelMapCollection::Creator>
PixelMapCollection::mapCreators;

bool
PixelMapCollection::creatorsInitialized=false;

void
PixelMapCollection::PixelMapTypeInitialize() {
  if (creatorsInitialized)
    return;
  creatorsInitialized = true;
  registerMapType<IdentityMap>();
  registerMapType<ReprojectionMap>();
  registerMapType<PolyMap>();
  registerMapType<LinearMap>();
  registerMapType<ConstantMap>();
}

//////////////////////////////////////////////////////////////
// Construction / destruction
//////////////////////////////////////////////////////////////

PixelMapCollection::PixelMapCollection() {
  PixelMapTypeInitialize();
  rebuildParameterVector(); // Sets the various counts to zero
}

PixelMapCollection::~PixelMapCollection() {
  // Destroy all of the compounded maps and submaps we've made:
  for (auto& i : mapElements) {
    if (i.second.atom) {
      delete i.second.atom;
      i.second.atom = nullptr;
    }
    if (i.second.realization) {
      delete i.second.realization;
      i.second.realization = nullptr;
    }
  }
  for (auto& i : wcsElements) {
    if (i.second.nativeCoords) {
      delete i.second.nativeCoords;
      i.second.nativeCoords = nullptr;
    }
    if (i.second.realization) {
      delete i.second.realization;
      i.second.realization = nullptr;
    }
  }
}

void
PixelMapCollection::removeMap(string mapName) {
  auto it = mapElements.find(mapName);
  if (it == mapElements.end())
    return; // Do nothing if there is no such map
  if (it->second.atom)
    delete it->second.atom;
  if (it->second.realization)
    delete it->second.realization;
  mapElements.erase(it);
}

void
PixelMapCollection::removeWcs(string wcsName) {
  auto it = wcsElements.find(wcsName);
  if (it == wcsElements.end())
    return; // Do nothing if there is no such wcs
  if (it->second.nativeCoords) {
    delete it->second.nativeCoords;
    it->second.nativeCoords=nullptr;
  }
  if (it->second.realization) {
    delete it->second.realization;
    it->second.realization = nullptr;
  }
  wcsElements.erase(it);
}

//////////////////////////////////////////////////////////////
// Parameter manipluation
//////////////////////////////////////////////////////////////

void
PixelMapCollection::setParams(const DVector& p) {
  Assert(p.size()==parameterCount);
  for (auto& melpair : mapElements) {
    MapElement& map = melpair.second;
    // Clear any defaulted flags
    map.isDefaulted = false;
    if (map.isFixed) continue;
    int nSub = map.nParams;
    if (nSub<=0) continue;
    DVector subp(nSub, 0.);
    subp = p.subVector(map.startIndex, map.startIndex+nSub);
    Assert(map.atom);
    map.atom->setParams(subp);
  }
}

DVector
PixelMapCollection::getParams() const {
  DVector p(parameterCount, 888.);
  for (auto& melpair : mapElements) {
    const MapElement& map = melpair.second;
    if (map.isFixed) continue;
    int nSub = map.nParams;
    if (nSub<=0) continue;
    if (!map.atom) 
      cerr << "mapElement is not atomic: " << melpair.first << " params: " << nSub << endl;
    Assert(map.atom);
    p.subVector(map.startIndex, map.startIndex+nSub) = 
      map.atom->getParams().subVector(0,nSub);
  }
  return p;
}

std::map<std::string, astrometry::DVector>
PixelMapCollection::getParamDict() const {
  map<string, linalg::Vector<double>> outMap;
  for (auto& melpair : mapElements) {
    const MapElement& map = melpair.second;
    if (map.isFixed) continue;
    int nSub = map.nParams;
    if (nSub<=0) continue;
    if (!map.atom) 
      cerr << "mapElement is not atomic: " << melpair.first << " params: " << nSub << endl;
    Assert(map.atom);
    DVector p = map.atom->getParams().subVector(0,nSub);
    string mapName = melpair.first;
    outMap.insert(std::pair<string, astrometry::DVector>(mapName, p));
    }
  return outMap;
}

std::string
PixelMapCollection::getMapType(string mapName) const {
  std::string mapType;
  auto melpair = mapElements.find(mapName);
  if (melpair == mapElements.end())
    throw AstrometryError("Could not find mapName " + mapName);
  
  const MapElement& mel = melpair->second;
  if (mel.atom) {
    mapType = mel.atom->getType();
  }
  else {
    mapType = SubMap::type();
  }
  return mapType;
}

DVector
PixelMapCollection::getWcsNativeCoords(string wcsName, bool degrees) const {
  auto wcspair = wcsElements.find(wcsName);
  if (wcspair == wcsElements.end())
    throw AstrometryError("Could not find wcsName " + wcsName);
  
  const WcsElement& wcs = wcspair->second;
  const SphericalCustomBase *coords = dynamic_cast<const SphericalCustomBase*>(wcs.nativeCoords);
  if (coords == nullptr) {
    throw AstrometryError("NativeCoords subclass does not have 'orient' property");
  }
  SphericalICRS pole = coords->getOrient()->getPole();
  DVector lonLat(2);
  pole.getRADec(lonLat[0], lonLat[1]);
  if (degrees) {
    lonLat[0] /= DEGREE;
    lonLat[1] /= DEGREE;
  }
  return lonLat;
}

void
PixelMapCollection::copyParamsFrom(const PixelMap& pm) {
  auto melpair = mapElements.find(pm.getName());
  if (melpair == mapElements.end())
    throw AstrometryError("Attempt to copyParamsFrom non-existent map element " + pm.getName());
  if (!melpair->second.atom)
    throw AstrometryError("Attempt to copyParamsFrom non-atomic map element " + pm.getName());
  melpair->second.atom->setParams(pm.getParams());
  melpair->second.isDefaulted = false;
}

void 
PixelMapCollection::setFixed(string name, bool isFixed) {
  set<string> s;
  s.insert(name);
  setFixed(s, isFixed);
}

void 
PixelMapCollection::setFixed(set<string> nameList, bool isFixed) {
  set<string> fixThese;
  for (auto iName : nameList) {
    set<string> addThese = dependencies(iName);
    fixThese.insert(addThese.begin(), addThese.end());
  }
  for (auto iName : fixThese) {
    auto& mel = mapElements[iName];
    // Update isFixed if this atom has free parameters
    if (mel.atom && mel.atom->nParams()>0) {
      mel.isFixed = isFixed;
    }
  }
  // Now calculate new parameter indices and propagate new configuration to all SubMaps
  rebuildParameterVector();
}

bool
PixelMapCollection::getFixed(string name) const {
  auto i = mapElements.find(name);
  if (i==mapElements.end()) {
    cerr << "No member of PixelMapCollection named " << name << endl;
    throw AstrometryError("No member of PixelMapCollection named " + name);
  }
  auto& mel = i->second;
  if (mel.atom) {
    return mel.isFixed;
  } else {
    for (auto j : mel.subordinateMaps) {
      if (!getFixed(j)) {
	// Any free member means compound is free
	return false;
      }
    }
    return true;  // Fixed if all subordinates are
  }
}

bool
PixelMapCollection::getDefaulted(string name) const {
  auto i = mapElements.find(name);
  if (i==mapElements.end())
    throw AstrometryError("No member of PixelMapCollection named " + name);
  auto& mel = i->second;
  if (mel.atom) {
    return mel.isDefaulted;
  } else {
    for (auto j : mel.subordinateMaps) {
      if (getDefaulted(j)) {
	// Any defaulted member means compound map is defaulted
	return true;
      }
    }
    return false;
  }
}

bool
PixelMapCollection::isAtomic(string name) const {
  auto i = mapElements.find(name);
  if (i==mapElements.end())
    throw AstrometryError("No member of PixelMapCollection named " + name);
  if (i->second.atom) {
    return true;
  } else {
    return false;
  }
}

void
PixelMapCollection::rebuildParameterVector() {
  // First assign new starting indices to every atomic map element

  // Restart the parameter counting:
  parameterCount = 0;
  // And map counting
  atomCount = 0;
  freeCount = 0;
  for (auto& i : mapElements) {
    //cerr << "rebuild mapElem" << endl;
    MapElement& map = i.second;
    // Only atomic map components go into the big parameter vector
    map.nParams = 0;
    map.number = -1;
    if (!(map.atom)) continue;
    atomCount++;
    int nMap = map.atom->nParams();
    if (nMap>0 && !map.isFixed) {
      // Have some parameters; append them to master list.
      map.nParams = nMap;
      map.startIndex = parameterCount;
      // Each atom with free parameters gets a serial number.
      map.number = freeCount;
      parameterCount += nMap;
      ++freeCount;
    }
  }

  // Then go through all extant SubMaps and update their versions of vectors & parameter counts.
  for (auto& i : mapElements) {
    MapElement& map = i.second;
    // Is there a SubMap here that we need to update?
    SubMap* sub = map.realization;
    if (!sub) continue;

    // Loop through all the map elements that this SubMap uses
    for (int iElement = 0; iElement < sub->vMaps.size(); iElement++) {
      string elementName = sub->vMaps[iElement]->getName();
      auto j = mapElements.find(elementName);
      if (j==mapElements.end()) 
	throw AstrometryError("PixelMapCollection::rebuildParameterVector could not find"
			      " map element with name " + elementName);
      if (!j->second.atom)
	throw AstrometryError("PixelMapCollection element " + elementName + " is not atomic"
			      " in rebuildParameterVector.");
      sub->vStartIndices[iElement] = j->second.startIndex;
      sub->vNSubParams[iElement] = j->second.nParams;
      sub->vMapNumbers[iElement] = j->second.number;
    }
    sub->countFreeParameters();
  }
}


//////////////////////////////////////////////////////////////
// Member maintenance
//////////////////////////////////////////////////////////////

vector<string>
PixelMapCollection::allMapNames() const {
  vector<string> output;
  for (auto& i : mapElements) {
    output.push_back(i.first);
  }
  return output;
}

vector<string>
PixelMapCollection::allWcsNames() const {
  vector<string> output;
  for (auto& i : wcsElements)
    output.push_back(i.first);
  return output;
}

void 
PixelMapCollection::learnMap(const PixelMap& pm, 
			     bool duplicateNamesAreExceptions,
			     bool rebuildIndices) {
  if (mapExists(pm.getName())) {
    if (duplicateNamesAreExceptions)
      throw AstrometryError("Duplicate map name in PixelMapCollection::learnMap at "
			    + pm.getName());
    // If not throwing an exception, we will just ignore this duplicate-named map.
    return;
  }
  const SubMap* sm = dynamic_cast<const SubMap*> (&pm);
  if (sm) {
    if (sm->vMaps.size()==1 && pm.getName()==sm->vMaps.front()->getName()) {
      // If this is a single-element SubMap that has same name as its dependent, then we'll
      // just learn the dependent:
      learnMap(*sm->vMaps.front(), duplicateNamesAreExceptions, rebuildIndices);
      return;
    } 
    // create a new compound map from this, learning each dependency:
    MapElement mel;
    for (int i = 0; i< sm->vMaps.size(); i++) {
      if (sm->getName() == sm->vMaps[i]->getName())
	throw AstrometryError("PixelMapCollection::learnMap encountered SubMap that has "
			      "a dependence on PixelMap with the same name: " + sm->getName());
      mel.subordinateMaps.push_back(sm->vMaps[i]->getName());
      learnMap(*sm->vMaps[i], duplicateNamesAreExceptions, rebuildIndices);
    }
    // Register this compound map
    mapElements.insert(std::pair<string,MapElement>(pm.getName(), mel));

  } else if ( const Wcs* wcsm = dynamic_cast<const Wcs*> (&pm)) {
    // A Wcs has been submitted as a PixelMap, which means it will be treated as
    // its pixel map followed by a reprojection.
    // Make sure we will not have a name conflict between this wcs and the PixelMap within it:
    if (pm.getName() == wcsm->getMap()->getName()) 
      throw AstrometryError("PixelMapCollection::learnMap received Wcs with same name as its "
			    " PixelMap: " + pm.getName());
    learnMap(*wcsm->getMap(), duplicateNamesAreExceptions, rebuildIndices);
    if (!wcsm->getTargetCoords()) 
      throw AstrometryError("PixelMapCollection::learnMap using a Wcs that has no target coord "
			    " system: " + wcsm->getName());
    string reprojectName = wcsm->getName() + "/reproject";
    ReprojectionMap repro(*wcsm->getNativeCoords(),
			  *wcsm->getTargetCoords(),
			  wcsm->getScale(),
			  reprojectName);
    learnMap(repro, duplicateNamesAreExceptions, rebuildIndices);
    MapElement mel;
    mel.subordinateMaps.push_back(wcsm->getMap()->getName());
    mel.subordinateMaps.push_back(reprojectName);
    // Register as compound map
    mapElements.insert(std::pair<string,MapElement>(pm.getName(), mel));
  } else {
    //  This should be an atomic PixelMap; simply add it to our element list:
    MapElement mel;
    mel.atom = pm.duplicate();
    mel.isFixed = mel.atom->nParams()==0;
    mapElements.insert(std::pair<string,MapElement>(pm.getName(), mel));
  }

  // Check that new map did not introduce circularity (though I do not know how it could)
  checkCircularDependence(pm.getName());

  // re-index everything
  if (rebuildIndices) rebuildParameterVector();
}

void
PixelMapCollection::learnWcs(const Wcs& wcs, 
			     bool duplicateNamesAreExceptions,
			     bool rebuildIndices) {
  // Check for duplicate Wcs name:
  if (wcsExists(wcs.getName())) {
    if (duplicateNamesAreExceptions) 
      throw AstrometryError("PixelMapCollection::learnWcs duplicates existing Wcs name " 
			    + wcs.getName());
    // We'll just ignore this if it's duplicated and we don't want to throw
    return;
  }
  // Learn the pixel map of the Wcs:
  learnMap(*wcs.getMap(), duplicateNamesAreExceptions, rebuildIndices);
  // Enter a new WcsElement into our tables:
  WcsElement wel;
  wel.mapName = wcs.getMap()->getName();
  wel.nativeCoords = wcs.getNativeCoords()->duplicate();
  wel.wScale = wcs.getScale();
  wcsElements.insert(std::pair<string,WcsElement>(wcs.getName(), wel));
}

void
PixelMapCollection::learn(PixelMapCollection& rhs, bool duplicateNamesAreExceptions) {
  if (&rhs == this) return;	// No need to learn self
  for (auto& iMap : rhs.mapElements) {
    const MapElement& incoming = iMap.second;
    if (mapExists(iMap.first)) {
      // incoming map duplicates and existing name.
      if (duplicateNamesAreExceptions) 
	throw AstrometryError("learn(PixelMapCollection) with duplicate map name "
			      + iMap.first);
      // Duplicate will just be ignored. 
    } else {
      // A new mapName for us.  Add its mapElement to our list.
      MapElement mel;
      if (iMap.second.atom) mel.atom = iMap.second.atom->duplicate();
      mel.isFixed = iMap.second.isFixed;
      mel.isDefaulted = iMap.second.isDefaulted;
      mel.subordinateMaps = iMap.second.subordinateMaps;
      mapElements.insert(std::pair<string,MapElement>(iMap.first, mel));
      // Note cannot introduce circularity if the top-level map is new.
    }
  } // end input map loop

  for (auto& iWcs : rhs.wcsElements) {
    WcsElement& incoming = iWcs.second;
    if (wcsExists(iWcs.first)) {
      // incoming wcs duplicates and existing name.
      if (duplicateNamesAreExceptions) 
	throw AstrometryError("learn(PixelMapCollection) with duplicate WCS name "
			      + iWcs.first);
      // Duplicate will just be ignored.  
    } else {
      // A new wcsName for us.  Add its info to our collection
      WcsElement wel;
      wel.mapName = iWcs.second.mapName;
      wel.nativeCoords = iWcs.second.nativeCoords->duplicate();
      wel.wScale = iWcs.second.wScale;
      wcsElements.insert(std::pair<string,WcsElement>(iWcs.first, wel));
    }
  } // end input wcs loop

  // And re-index everything
  rebuildParameterVector();
}

// Define a new pixelMap that is compounding of a list of other PixelMaps.  Order
// of the list is from pixel to world coordinates.
void 
PixelMapCollection::defineChain(string chainName, const list<string>& elements) {
  if (mapExists(chainName)) 
    throw AstrometryError("PixelMapCollection::defineChain with duplicate name: " 
			  + chainName);
  // Check that elements exist
  for (auto& i : elements) {
    if (!mapExists(i))
      throw AstrometryError("PixelMapCorrection::defineChain with unknown pixel map element: "
			    + i);
  }
  mapElements.insert(std::pair<string, MapElement>(chainName, MapElement()));
  mapElements[chainName].subordinateMaps = elements;
  // Note that adding a new chain does not change parameter vector assignments
  // Nor can it introduce circular dependence if chain name is new and all elements exist.
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
  if (!mapExists(mapName))
    throw AstrometryError("PixelMapCollection::issueMap requested for unknown PixelMap: "
			  +  mapName);
  MapElement& el = mapElements[mapName];

  if (!el.realization) {
    // Create a realization if one does not exist
    list<string> atomList = orderAtoms(mapName);
    vector<int> startIndices(atomList.size(), 0);
    vector<int> nParams(atomList.size(), 0);
    vector<int> mapNumbers(atomList.size(), -1);

    list<PixelMap*> atoms;
    int index=0;
    for (auto i = atomList.begin();
	 i != atomList.end();
	 ++i, ++index) {
      auto& mel = mapElements[*i];
      Assert(mel.atom);	// All elements should be atomic
      atoms.push_back(mel.atom);
      // fill in its indices into master vector:
      startIndices[index] = mel.startIndex;
      nParams[index] = mel.isFixed ? 0 : mel.nParams;
      mapNumbers[index] = mel.number;
    }
    SubMap* sm = new SubMap(atoms, mapName, true);
    sm->vStartIndices = startIndices;
    sm->vNSubParams = nParams;
    sm->vMapNumbers = mapNumbers;
    sm->countFreeParameters();
    el.realization = sm;
  }
  return el.realization;
}

PixelMap*
PixelMapCollection::cloneMap(string mapName) const {
  if (!mapExists(mapName))
    throw AstrometryError("PixelMapCollection::issueMap requested for unknown PixelMap: "
			  +  mapName);
  const MapElement& el = mapElements.find(mapName)->second;

  if (el.atom) 
    return el.atom->duplicate();

  // A composite one we will create as a SubMap:
  list<PixelMap*> vMaps;
  for (auto& i : el.subordinateMaps)
    vMaps.push_back(cloneMap(i));
  SubMap* retval = new SubMap(vMaps, mapName, false);
  // Clean up the stray clones since they've been duplicated by SubMap:
  for (auto i : vMaps)
    delete i;
  return retval;
}

// Return a pointer to a Wcs built from a SubMap and realizing the named coord system
Wcs* 
PixelMapCollection::issueWcs(string wcsName) {
  if (!wcsExists(wcsName))
    throw AstrometryError("PixelMapCollection::issueWcs requested for unknown Wcs: "
			  +  wcsName);
  WcsElement& el = wcsElements[wcsName];

  if (!el.realization) {
    // Create a realization if one does not exist.  Wcs set to share its PixelMap
    SubMap* sm = issueMap(el.mapName);
    el.realization = new Wcs(sm, *el.nativeCoords, wcsName, el.wScale, true);
  }
  return el.realization;
}

// Return a pointer to a Wcs built from a SubMap and realizing the named coord system
Wcs* 
PixelMapCollection::cloneWcs(string wcsName) const {
  if (!wcsExists(wcsName))
    throw AstrometryError("PixelMapCollection::cloneWcs requested for unknown Wcs: "
			  +  wcsName);
  const WcsElement& el = wcsElements.find(wcsName)->second;
  PixelMap* pm = cloneMap(el.mapName);
  Wcs* retval = new Wcs(pm, *el.nativeCoords, wcsName, el.wScale, false);
  delete pm;
  return retval;
}

set<string>
PixelMapCollection::dependencies(string target) const {
  // Find target MapElement and add target to output dependency list.  
  // Assume no circular dependence.
  auto iTarget = mapElements.find(target);
  if (iTarget == mapElements.end()) 
    throw AstrometryError("PixelMapCollection has no element " + target + " in dependencies()");
  set<string> output;
  output.insert(target);

  // Call this routine recursively for all the map elements that the target depends on
  for (auto& i : iTarget->second.subordinateMaps) {
    set<string> subs = dependencies(i);
    output.insert(subs.begin(), subs.end());
  }

  return output;
}

bool
PixelMapCollection::dependsOn(const string mapName, const string targetName) const {
  auto melptr = mapElements.find(mapName);
  if (melptr == mapElements.end()) 
    throw AstrometryError("PixelMapCollection has no element " + mapName + " in dependsOn()");

  if (mapName==targetName) return true;  // True if map is the target

  // Call this routine recursively for all the map elements that the target depends on
  for (auto& i : melptr->second.subordinateMaps)
    if (dependsOn(i, targetName)) return true;

  return false;
}

// Return a list of the atomic elements of the named map, in order of
// their application to data.  Assumes no circular dependence.
list<string> 
PixelMapCollection::orderAtoms(string mapName) const {
  if (!mapExists(mapName))
    throw AstrometryError("PixelMapCollection::orderAtoms requested for unknown map: "
			  + mapName);
  list<string> retval;
  const MapElement& m = mapElements.find(mapName)->second;
  if (m.atom) {
    // Simple atomic element just returns itself
    retval.push_back(mapName);
    return retval;
  }

  // Otherwise we have a compound map at work here.
  // Call recursively to all subordinateMaps:
  for (auto& i : m.subordinateMaps) {
    list<string> subList = orderAtoms(i);
    retval.splice(retval.end(), subList);
  }
  return retval;
}

void
PixelMapCollection::checkCircularDependence(string mapName,
					    const set<string>& ancestors) const {
  if (!mapExists(mapName))
    throw AstrometryError("Unknown mapName in PixelMapCollection::checkCircularDependency: "
			  + mapName);
  if (ancestors.count(mapName))
    throw AstrometryError("Circular dependency in PixelMapCollection at map "
			  + mapName);
  // Then call recursively to all subordinateMaps, adding self to ancestor list:
  set<string> ancestorsPlusSelf(ancestors);
  ancestorsPlusSelf.insert(mapName);
  const MapElement& m = mapElements.find(mapName)->second;
  for (auto& i : m.subordinateMaps) 
    checkCircularDependence(i, ancestorsPlusSelf);
}

void
PixelMapCollection::checkCompleteness() const {
  // Loop through all the (non-atomic) maps, throw exception if they refer to 
  // non-existent map elements or to themselves
  for (auto& i : mapElements) 
    checkCircularDependence(i.first);

  // Loop through all the Wcs's, making sure their pixelmaps exist.
  for (auto& i : wcsElements) {
    if (!mapExists(i.second.mapName))
      throw AstrometryError("PixelMapCollection Wcs <" + i.first + "> refers to unknown PixelMap <"
			    + i.second.mapName + ">");
  }
}

void
PixelMapCollection::invalidate(string wcsName) {
  auto it = wcsElements.find(wcsName);
  if (it==wcsElements.end()) return;
  it->second.isValid = false;
}

template <class T>
void
set_subtract(std::set<T>& s1, const std::set<T> s2) {
  auto i1 = s1.begin();
  auto i2 = s2.begin();
  while(i1!=s1.end() && i2!=s2.end()) {
    if (*i1 < *i2) {
      ++i1;
    } else if (*i2 < *i1) {
      ++i2;
    } else {
      // Equality: remove from s1
      i1 = s1.erase(i1);
      ++i2;
    }
  }
}
				       
void
PixelMapCollection::purgeInvalid() {
  // Collect names of all maps that invalid WCS's use,
  // and get rid of the WCS's themselves.
  set<string> unneeded;
  set<string> badWcs;
  for (auto& i : wcsElements) {
    if (i.second.isValid) continue;
    badWcs.insert(i.first);
    if (mapExists(i.first)) {
      // The WCS is also a map, build the dependency from there
      auto depend = dependencies(i.first);
      unneeded.insert(depend.begin(), depend.end());
    } else if (!i.second.mapName.empty() &&
	       i.second.mapName!="Identity") {
      // purge dependencies of the WCS's map
      auto depend = dependencies(i.second.mapName);
      unneeded.insert(depend.begin(), depend.end());
    }
  }
  // Now go through all PixelMaps.  If they are
  // not dependencies of the badWcs's, then
  // we keep everything that they need.
  for (auto& i : mapElements) {
    if (unneeded.count(i.first)) continue; // It's a dependency of our WCS.
    // If not part of a bad WCS, then it's part of a good one,
    // and we want to keep everything in it.
    set_subtract(unneeded, dependencies(i.first));
  }
  // Now kill all unneeded maps
  for (auto s : unneeded) {
    removeMap(s);
  }
  // And wcs's too
  for (auto s : badWcs) {
    removeWcs(s);
  }

  // Recalculate all parameter indices
  rebuildParameterVector();
}

//////////////////////////////////////////////////////////////
// YAML (De-) Serialization
//////////////////////////////////////////////////////////////

void
PixelMapCollection::writeSingleWcs(YAML::Emitter& os, const WcsElement& wel, string name) const {
  // This routine assumes we are already in the "WCS" node 
  os << YAML::Key << name
     << YAML::Value
     << YAML::BeginMap
     << YAML::Key << "Type" << YAML::Value << Wcs::type()
     << YAML::Key << "MapName" << YAML::Value << wel.mapName
     << YAML::Key << "Projection" << YAML::Value << *(wel.nativeCoords)
     << YAML::Key << "Scale" << YAML::Value << wel.wScale
     << YAML::EndMap;
}


void
PixelMapCollection::writeSingleMap(YAML::Emitter& os, const MapElement& mel, string name)  const {
  if (mel.atom) {
    // Atomic map, give all details:
    const ColorTerm* ct = dynamic_cast<const ColorTerm*> (mel.atom);
    if (ct) {
      // We have a color term which is wrapping another map.
      // Make a wrapper node of Type=Color and write the wrapped function
      os << YAML::Key << name
	 << YAML::Value
	 << YAML::BeginMap
	 << YAML::Key << "Type" << YAML::Value << ColorTerm::type()
	 << YAML::Key << "Reference" << YAML::Value << ct->refColor()
	 << YAML::Key << "Function" << YAML::Value;
      ct->map()->write(os);
      os << YAML::EndMap;
    } else {
      // A simple atomic node.  Let it write itself
      os << YAML::Key << name
	 << YAML::Value;
      mel.atom->write(os);
    }
  } else {
    // For compound map, just give names of constituents
    os << YAML::Key << name
       << YAML::Value
       << YAML::BeginMap
       << YAML::Key << "Type" << YAML::Value << SubMap::type()
       << YAML::Key << "Elements"
       << YAML::BeginSeq;
    for (auto subname : mel.subordinateMaps) {
      os << subname;
    }
    os << YAML::EndSeq
       << YAML::EndMap;
  }
}

void 
PixelMapCollection::write(ostream& os, string comment) const {
  // Create a YAML emitter.  Produce a map where one
  // key is the magic word, and the rest are names of maps.
  YAML::Emitter out;
  out << YAML::BeginMap;
  if (comment.size()>0)
    out << YAML::Comment(comment);
  out << YAML::Key << magicKey
      << YAML::Value << "This is a serialized PixelMapCollection";
  for (auto& melpair : mapElements) {
    writeSingleMap(out, melpair.second, melpair.first);
  }
  if (!wcsElements.empty()) {
    // Make a "WCS" node under the root
    out << YAML::Key << "WCS"
	<< YAML::Value
	<< YAML::BeginMap;
    for (auto& wcspair : wcsElements) {
      writeSingleWcs(out, wcspair.second, wcspair.first);
    }
    out << YAML::EndMap;  // Close the WCS node
  }
  out << YAML::EndMap;
  // Send the YAML to the output stream
  os << out.c_str() << endl;
}

void 
PixelMapCollection::writeMap(ostream& os, string name, string comment) const {
  if (!mapExists(name))
    throw AstrometryError("PixelMapCollection::writeMap() for unknown map: " + name);
  set<string> allmaps = dependencies(name);
  YAML::Emitter out;
  out << YAML::BeginMap;
  if (comment.size()>0)
    out << YAML::Comment(comment);
  out << YAML::Key << magicKey
      << YAML::Value << "This is a serialized PixelMapCollection";
  // Write the map elements
  for (auto mapname : allmaps) {
    auto j = mapElements.find(mapname);
    writeSingleMap(out, j->second, mapname);
  }
  out << YAML::EndMap;
  // Send the YAML to the output stream
  os << out.c_str() << endl;
}

void 
PixelMapCollection::writeWcs(ostream& os, string name, string comment) const {
  if (!wcsExists(name))
    throw AstrometryError("PixelMapCollection::writeWcs() for unknown Wcs: " + name);
  auto j = wcsElements.find(name);
  YAML::Emitter out;
  out << YAML::BeginMap;
  if (comment.size()>0)
    out << YAML::Comment(comment);
  out << YAML::Key << magicKey
      << YAML::Value << "This is a serialized PixelMapCollection";
  // Write the Wcs information:
  out << YAML::Key << "WCS"
      << YAML::Value
      << YAML::BeginMap;
  writeSingleWcs(out, j->second, name);
  out << YAML::EndMap; // close out the WCS node.
  // Then all of the dependent map elements
  set<string> allmaps = dependencies(j->second.mapName);
  // Write the map elements
  for (auto mapname : allmaps) {
    auto ii = mapElements.find(mapname);
    writeSingleMap(out, ii->second, mapname);
  }
  out << YAML::EndMap;
  // Send the YAML to the output stream
  os << out.c_str() << endl;
}

bool
PixelMapCollection::read(istream& is, string namePrefix) {
  string buffer;
  //try {
  {

    YAML::Node root = YAML::Load(is);
    if (!root.IsMap() || !root[magicKey]) {
      // A valid YAML document that does not have the structure 
      // and key that our files should have.
      return false;
    }

    // Every entry in the root node map is a map to enter
    for (YAML::const_iterator i = root.begin();
	 i != root.end();
	 ++i) {
      // Skip the map key that is our magic word
      string keyName = i->first.as<string>();
      if (keyName == magicKey) continue;

      string name = namePrefix + keyName;
      
      const YAML::Node& node=i->second;
      // If there is a node named "WCS" it holds all of our WCS definitions
      if (keyName=="WCS") {
	if (!node.IsMap()) 
	  throw AstrometryError("YAML node 'WCS' is not a map");
	// Loop through its entries
	for (auto j : node) {
	  string wName = namePrefix + j.first.as<string>();
	  auto& wNode=j.second;
	  // All nodes here should be of type WCS
	  if (!wNode.IsMap() || !wNode["Type"] ||
	      !stringstuff::nocaseEqual(wNode["Type"].as<string>(), Wcs::type())) { 
	    throw AstrometryError("YAML node " + wName + " in WCS node is not of type WCS");
	  }
	  if (!wNode["MapName"] || !wNode["Projection"] || !wNode["Scale"])
      throw AstrometryError("Missing YAML keys for deserializing Wcs " + wName);
	  // Enter a new WcsElement into our tables:
	  WcsElement wel;
	  wel.mapName = wNode["MapName"].as<string>();
    //wel.wScale = wNode["Scale"].as<double>();
	  wel.wScale = stod(wNode["Scale"].as<string>());
	  wel.nativeCoords = SphericalCoords::deserialize(wNode["Projection"]);
	  wcsElements.insert(std::pair<string,WcsElement>(wName, wel));
	} // Done parsing a WCS node
	continue;  // Move on to the next subnode of root.
      } // Done parsing the whole WCS map.
      
      // All further nodes should be maps:
      if (!node.IsMap() || !node["Type"]) 
	throw AstrometryError("Non-map or missing Type key in PixelMapCollection YAML entry " +
			      name);
      string mapType = node["Type"].as<string>();

      if (stringstuff::nocaseEqual(mapType, Wcs::type())) {
	// All the WCS's should have been in the WCS node.
	throw AstrometryError("Node " + name + " of type WCS found outside the WCS library");
      }	else if (stringstuff::nocaseEqual(mapType,ColorTerm::type())) {
	// A color term has the function it modifies in a child node:
  if (!node["Function"])
	throw AstrometryError("Missing Function key in ColorTerm YAML entry " +
			      name);
	double reference = 0.; // Default reference color
	if (node["Reference"])
	  reference = node["Reference"].as<double>();
	bool defaulted;
	PixelMap* pm = createAtomFromNode(node["Function"], defaulted, "");
	defaulted = false;  // *** Always going to assume default values ok for color terms**/

	// create the MapElement
	MapElement mel;
	mel.atom = new ColorTerm(pm, reference, name);
	mel.isDefaulted = defaulted;
	mel.isFixed = pm->nParams()==0;
	mapElements.insert(std::pair<string,MapElement> (name, mel));
	
      } else if (stringstuff::nocaseEqual(mapType,SubMap::type())) {
	// Composite maps: all names stored in a child node
	if (!node["Elements"])
	  throw AstrometryError("Missing Elements key in SubMap YAML entry " + name);
	list<string> submaps = node["Elements"].as<list<string> >();

	// Create the MapElement
	MapElement mel;
	mel.subordinateMaps = submaps;
	mapElements.insert(std::pair<string,MapElement>(name, mel));

      } else {
	// Should have some kind of atomic map.  Call a function to parse the node.
	// This will throw exception if we have an unknown Type of map
	  
	// create the MapElement
	MapElement mel;
	bool defaulted;
	mel.atom = createAtomFromNode(node,defaulted, name);
	mel.isDefaulted = defaulted;
	mel.isFixed = mel.atom->nParams()==0;  // Set fixed if no free params
	mapElements.insert(std::pair<string,MapElement> (name, mel));
      }

    } // end loop through root node map members

  }// catch (YAML::Exception& e) {
  //  /**/cerr << "PixelMapCollection::read() caught: " << e.what() << endl;
  //  // Return false if not a valid YAML file
  //  return false;
  //}

  // Check all maps for circular dependence and completeness.
  checkCompleteness();

  // Recalculate indices
  rebuildParameterVector();
  return true;
}

PixelMap*
PixelMapCollection::createAtomFromNode(const YAML::Node& node,
				       bool& defaulted,
				       string name) const {
  // Access the static arrays of map type names and creators
  if (!node["Type"])
    throw AstrometryError("Missing Type key in createAtomFromNode at " + name);

  string mapType = node["Type"].as<string>();
  auto it = mapCreators.find(mapType);
  if (it==mapCreators.end())
    throw AstrometryError("PixelMapCollection does not recognize PixelMap Type " + mapType);
  PixelMap* out = (*(it->second))(node, defaulted, name);
  return out;
}

string
PixelMapCollection::atomHavingParameter(int parameterIndex) const {
  // Find the atomic map that this parameter number influences
  for (auto& melpair : mapElements) {
    if (!melpair.second.atom) continue;
    if ( parameterIndex >= melpair.second.startIndex &&
	 parameterIndex < melpair.second.startIndex + melpair.second.nParams)
      return melpair.first;
  }
  return "";	// Nothing in range, return empty string.
}

void
PixelMapCollection::parameterIndicesOf(string mapname,
				       int& startIndex,
				       int& nParams) const {
  if (!mapExists(mapname))
    throw AstrometryError("PixelMapCollection::parameterIndicesof() unknown map: "
			  + mapname);
  startIndex = 0;
  nParams = 0;
  auto& mel = mapElements.at(mapname);
  if (!mel.atom || mel.isFixed) return;
  startIndex = mel.startIndex;
  nParams = mel.nParams;
  return;
}

#endif // USE_YAML
