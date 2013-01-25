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
  int mapNumber = 0;
  for (MapIter i = mapElements.begin(); i!=mapElements.end(); ++i, ++mapNumber) {
    MapElement& map = i->second;
    // Only atomic map components go into the big parameter vector
    if (!(map.atom)) continue;
    int nMap = map.atom->nParams();
    if (map.isFixed || nMap==0) {
      // No free parameters here.
      map.nParams = 0;
      map.number = -1;
    } else {
      // Have some parameters; append them to master list.
      map.nParams = nMap;
      map.startIndex = parameterCount;
      map.number = mapNumber;
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
PixelMapCollection::learnMap(const PixelMap& pm, bool duplicateNamesAreExceptions) {
  if (mapExists(pm.getName())) {
    if (duplicateNamesAreExceptions)
      throw AstrometryError("Duplicate map name in PixelMapCollection::learnMap at "
			    + pm.getName());
    // If not throwing an exception, we will just ignore this duplicate-named map.
    return;
  }
  const SubMap* sm = dynamic_cast<const SubMap*> (&pm);
  const CompoundPixelMap* cpm = dynamic_cast<const CompoundPixelMap*> (&pm);
  if (sm) {
    if (sm->vMaps.size()==1 && pm.getName()==sm->vMaps.front()->getName()) {
      // If this is a single-element SubMap that has same name as its dependent, then we'll
      // just learn the dependent:
      learnMap(*sm->vMaps.front(), duplicateNamesAreExceptions);
      return;
    } 
    // create a new compound map from this, learning each dependency:
    MapElement mel;
    for (int i = 0; i< sm->vMaps.size(); i++) {
      if (sm->getName() == sm->vMaps[i]->getName())
	throw AstrometryError("PixelMapCollection::learnMap encountered SubMap that has "
			      "a dependence on PixelMap with the same name: " + sm->getName());
      mel.subordinateMaps.push_back(sm->vMaps[i]->getName());
      learnMap(*sm->vMaps[i], duplicateNamesAreExceptions);
    }
    // Register this compound map
    mapElements.insert(std::pair<string,MapElement>(pm.getName(), mel));

  } else if (cpm) {
    // learn each dependency
    MapElement mel;
    const list<PixelMap*>& pmlist = cpm->getList();
    for (list<PixelMap*>::const_iterator i = pmlist.begin();
	 i != pmlist.end();
	 ++i) {
      if (cpm->getName() == (*i)->getName())
	throw AstrometryError("PixelMapCollection::learnMap encountered CompoundPixelMap that "
			      "has a dependence on PixelMap with the same name: "
			      + cpm->getName());
      mel.subordinateMaps.push_back((*i)->getName());
      learnMap(**i, duplicateNamesAreExceptions);
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
    learnMap(*wcsm->getMap(), duplicateNamesAreExceptions);
    if (!wcsm->getTargetCoords()) 
      throw AstrometryError("PixelMapCollection::learnMap using a Wcs that has no target coord "
			    " system: " + wcsm->getName());
    string reprojectName = wcsm->getName() + "_reproject";
    ReprojectionMap repro(*wcsm->getNativeCoords(),
			  *wcsm->getTargetCoords(),
			  wcsm->getScale(),
			  reprojectName);
    learnMap(repro, duplicateNamesAreExceptions);
    MapElement mel;
    mel.subordinateMaps.push_back(wcsm->getMap()->getName());
    mel.subordinateMaps.push_back(reprojectName);
    // Register as compound map
    mapElements.insert(std::pair<string,MapElement>(pm.getName(), mel));
  } else {
    //  This should be an atomic PixelMap; simply add it to our element list:
    MapElement mel;
    mel.atom = pm.duplicate();
    mapElements.insert(std::pair<string,MapElement>(pm.getName(), mel));
  }

  // Check that new map did not introduce circularity (though I do not know how it could)
  checkCircularDependence(pm.getName());

  // re-index everything
  rebuildParameterVector();
}

void
PixelMapCollection::learnWcs(const Wcs& wcs, bool duplicateNamesAreExceptions) {
  // Check for duplicate Wcs name:
  if (wcsExists(wcs.getName())) {
    if (duplicateNamesAreExceptions) 
      throw AstrometryError("PixelMapCollection::learnWcs duplicates existing Wcs name " 
			    + wcs.getName());
    // We'll just ignore this if it's duplicated and we don't want to throw
    return;
  }
  // Learn the pixel map of the Wcs:
  learnMap(*wcs.getMap());
  // Enter a new WcsElement into our tables:
  WcsElement wel;
  wel.mapName = wcs.getMap()->getName();
  wel.nativeCoords = wcs.getNativeCoords()->duplicate();
  wel.wScale = wcs.getScale();
  wcsElements.insert(std::pair<string,WcsElement>(wcs.getName(), wel));
}

void
PixelMapCollection::learn(PixelMapCollection& rhs, bool duplicateNamesAreExceptions) {
  for (ConstMapIter iMap = rhs.mapElements.begin();
       iMap != rhs.mapElements.end();
       ++iMap) {
    const MapElement& incoming = iMap->second;
    if (mapExists(iMap->first)) {
      // incoming map duplicates and existing name.
      if (duplicateNamesAreExceptions) 
	throw AstrometryError("learn(PixelMapCollection) with duplicate map name "
			      + iMap->first);
      // Duplicate will just be ignored. 
    } else {
      // A new mapName for us.  Add its mapElement to our list.
      MapElement mel;
      mel.atom = iMap->second.atom->duplicate();
      mel.isFixed = iMap->second.isFixed;
      mel.subordinateMaps = iMap->second.subordinateMaps;
      mapElements.insert(std::pair<string,MapElement>(iMap->first, mel));
      // Note cannot introduce circularity if the top-level map is new.
    }
  } // end input map loop

  for (WcsIter iWcs = rhs.wcsElements.begin();
       iWcs != rhs.wcsElements.end();
       ++iWcs) {
    WcsElement& incoming = iWcs->second;
    if (wcsExists(iWcs->first)) {
      // incoming wcs duplicates and existing name.
      if (duplicateNamesAreExceptions) 
	throw AstrometryError("learn(PixelMapCollection) with duplicate WCS name "
			      + iWcs->first);
      // Duplicate will just be ignored.  
    } else {
      // A new wcsName for us.  Add its info to our collection
      WcsElement wel;
      wel.mapName = iWcs->second.mapName;
      wel.nativeCoords = iWcs->second.nativeCoords->duplicate();
      wel.wScale = iWcs->second.wScale;
      wcsElements.insert(std::pair<string,WcsElement>(iWcs->first, wel));
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
  for (list<string>::const_iterator i = elements.begin();
       i != elements.end();
       ++i) 
    if (!mapExists(*i))
      throw AstrometryError("PixelMapCorrection::defineChain with unknown pixel map element: "
			    + *i);
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
    for (list<string>::const_iterator i = atomList.begin();
	 i != atomList.end();
	 ++i, ++index) {
      Assert(mapElements[*i].atom);	// All elements should be atomic
      atoms.push_back(mapElements[*i].atom);
      // fill in its indices into master vector:
      startIndices[index] = mapElements[*i].startIndex;
      nParams[index] = mapElements[*i].isFixed ? 0 : mapElements[*i].nParams;
      mapNumbers[index] = mapElements[*i].number;
    }
    SubMap* sm = new SubMap(atoms, mapName);
    sm->vStartIndices = startIndices;
    sm->vNSubParams = nParams;
    sm->vMapNumbers = mapNumbers;
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
    // Create a realization if one does not exist.  Wcs set to share its PixelMap
    SubMap* sm = issueMap(el.mapName);
    el.realization = new Wcs(sm, *el.nativeCoords, wcsName, el.wScale, true);
  }
  return el.realization;
}

set<string>
PixelMapCollection::dependencies(string target) const {
  // Find target MapElement and add target to output dependency list.  
  // Assume no circular dependence.
  ConstMapIter iTarget = mapElements.find(target);
  if (iTarget == mapElements.end()) 
    throw AstrometryError("PixelMapCollection has no element " + target + " in dependencies()");
  set<string> output;
  output.insert(target);

  // Call this routine recursively for all the map elements that the target depends on
  for (list<string>::const_iterator i = iTarget->second.subordinateMaps.begin();
       i != iTarget->second.subordinateMaps.end();
       ++i) {
    set<string> subs = dependencies(*i);
    output.insert(subs.begin(), subs.end());
  }

  return output;
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
  for (list<string>::const_iterator i = m.subordinateMaps.begin();
       i != m.subordinateMaps.end();
       ++i) {
    list<string> subList = orderAtoms(*i);
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
  for (list<string>::const_iterator i = m.subordinateMaps.begin();
       i != m.subordinateMaps.end();
       ++i) 
    checkCircularDependence(*i, ancestorsPlusSelf);
}


//////////////////////////////////////////////////////////////
// (De-) Serialization
//////////////////////////////////////////////////////////////

const string magicHeader = "PixelMapCollection";
void 
PixelMapCollection::read(istream& is, string namePrefix) {
  // ???
}

void
PixelMapCollection::writeSingleWcs(const WcsElement& wel, string name, ostream& os) const {
  os << Wcs::mapType()
     << " " << name
     << " " << wel.mapName
     << endl;
  const Gnomonic* gn = dynamic_cast<const Gnomonic*> (wel.nativeCoords);
  if (!gn)
    throw AstrometryError("PixelMapCollection can only serialize Gnomonic projections at Wcs name "
			  + name);
  os << *gn->getOrient() << endl;
  os << wel.wScale << endl;
}

void
PixelMapCollection::writeSingleMap(const MapElement& mel, string name, ostream& os)  const {
  if (mel.atom) {
    // Atomic map, give all details:
    os << mel.atom->getType() 
       << " " << name
       << " " << mel.atom->getPixelStep()
       << endl;
    mel.atom->write(os);
  } else {
    // For compound map, just give names of constituents
    os << SubMap::mapType()
       << " " << name
       << " " << mel.subordinateMaps.size()
       << endl;
    const int maxOnLine = 5;
    int iOnLine = 0;
    for (list<string>::const_iterator i=mel.subordinateMaps.begin();
	 i != mel.subordinateMaps.end();
	 ++i) {
      if (iOnLine >= maxOnLine) {
	os << endl;
	iOnLine = 0;
      }
      os << *i << " ";
      iOnLine++;
    }
    os << endl;
  }
}

void 
PixelMapCollection::write(ostream& os) const {
  os << magicHeader << endl;
  // Write the map elements
  for (ConstMapIter i = mapElements.begin();
       i != mapElements.end();
       ++i) {
    writeSingleMap(i->second, i->first, os);
  }
  // Then write the Wcs's:
  for (ConstWcsIter i = wcsElements.begin();
       i != wcsElements.end();
       ++i) {
    writeSingleWcs(i->second, i->first, os);
  }
}

void 
PixelMapCollection::writeMap(ostream& os, string name) const {
  // ???
}
void 
PixelMapCollection::writeWcs(ostream& os, string name) const {
  // ???
}
