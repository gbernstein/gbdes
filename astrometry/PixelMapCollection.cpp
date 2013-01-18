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
PixelMapCollection::deleteMap(PixelMap* pm, const set<PixelMap*>& ancestors) const {
  if (ancestors.find(pm))
    throw AstrometryError("Circular PixelMap dependence detected in "
			  "PixelMapCollection::deleteMap at " + pm->getName());
  if ( (SubMap* sm = dynamic_cast<SubMap*> (pm))!=0) {
    // Should not be deleting a SubMap that is owned by another PixelMapCollection:
    if (sm->isOwned)
      throw AstrometryError("PixelMapCollection::deleteMap called on owned SubMap "
			    + pm->getName());
    // Call recursively on all elements of SubMap:
    set<PixelMap*> ancestorsPlusSelf(ancestors);
    ancestorsPlusSelf.insert(pm);
    for (int i = 0; i< sm->vMaps.size(); i++)
      deleteMap(sm->vMaps[i], ancestorsPlusSelf);
  } else if ( (Wcs* wcsm = dynamic_cast<wcs*> (pm))!=0) {
    // If the PixelMap is a Wcs, call on its associated pixelMap
    deleteMap(wcsm->getMap(), ancestors);
  }
  // See if this pm is referenced by any part of this PixelMapCollection
  bool deleteIt = true;
  for (ConstMapIter i = mapElements.begin();
       i != mapElements.end();
       ++i) 
    if (i->second.atom == pm) {
      deleteIt = false;
      break;
    }
  // If not in use, delete it.
  if (deleteIt) delete pm;
}

void 
PixelMapCollection::absorbMap(PixelMap* pm, bool duplicateNamesAreExceptions) {
  if (mapExists(pm->getName())) {
    if (duplicateNamesAreExceptions)
      throw AstrometryError("Duplicate map name in PixelMapCollection::absorbMap at "
			    + pm->getName());
    // We will just ignore this duplicate-named map.
    // Delete it and its dependencies if they are not in this PixelMapCollection:
    deleteMap(pm);
    return;
  }
  if ( (SubMap* sm = dynamic_cast<SubMap*> (pm))!=0) {
    // Should not be absorbing a SubMap that is owned by another PixelMapCollection:
    if (sm->isOwned)
      throw AstrometryError("PixelMapCollection::absorbMap called on already-owned SubMap "
			    + pm->getName());
    // If this is a single-element SubMap that has same name as its dependent, then we'll
    // just absorb the dependent:
    if (sm->vMaps.size()==1 && pm->getName()==sm->vMaps.front()->getName()) {
      absorbMap(sm->vMaps.front(), duplicateNamesAreExceptions);
      // Delete this SubMap and be done
      delete sm;
    } else {
      // create a new compound map from this, absorbing each dependency:
      MapElement mel;
      for (int i = 0; i< sm->vMaps.size(); i++) {
	if (sm->getName() == sm->vMaps[i]->getName())
	  throw AstrometryError("PixelMapCollection::absorbMap encountered SubMap that has "
				"a dependence on PixelMap with the same name: " + sm->getName());
	mel.subordinateMaps.push_back(sm->vMaps[i]->getName());
	absorbMap(sm->vMaps[i], duplicateNamesAreExceptions);
      }
      // Register this compound map
      mapElements.insert(std::pair<string,MapElement>(pm->getName(), mel));
      // Delete the original SubMap
      delete sm;
    }
  } else   if ( (WcsMap* wcsm = dynamic_cast<WcsMap*> (pm))!=0) {
    // A Wcs has been submitted as a PixelMap, which is just too confusing to deal with
    throw AstrometryError("Passed the Wcs <" + pm->getName() "> as a PixelMap to"
			  " PixelMapCollection::absorbMap.\n"
			  "This is not supported, use defineWcs and/or ReprojectionMaps.");
  } else if ( (CompoundPixelMap* cpm = dynamic_cast<CompoundPixelMap*> (pm))!=0) {
    // absorb each dependency:
    MapElement mel;
    for (list<PixelMap*>::iterator i = cpm->pmlist.begin();
	 i != cpm->pmlist.end();
	 ++i) {
      if (cpm->getName() == (*i)->getName())
	throw AstrometryError("PixelMapCollection::absorbMap encountered CompoundPixelMap that "
			      "has a dependence on PixelMap with the same name: "
			      + cpm->getName());
      mel.subordinateMaps.push_back((*i)->getName());
      absorbMap(*i, duplicateNamesAreExceptions);
    }
    // Register this compound map
    mapElements.insert(std::pair<string,MapElement>(pm->getName(), mel));
    // Delete the original CompoundPixelMap
    delete cpm;
  } else {
    //  This should be an atomic PixelMap; simply add it to our element list:
    MapElement mel;
    mel.atomic = pm;
    mapElements.insert(std::pair<string,MapElement>(pm->getName(), mel));
  }

  // re-index everything
  rebuildParameterVector();
}

void
PixelMapCollection::absorb(PixelMapCollection& rhs, bool duplicateNamesAreExceptions) {
  set<PixelMap*> atomsToKill;
  for (MapIter iMap = rhs.mapElements.begin();
       iMap != rhs.mapElements.end();
       ++iMap) {
    MapElement& incoming = iMap->second;
    // Kill any existing realization of incoming map
    if (incoming.realization) {
      delete incoming.realization;
      incoming.realization = 0;
    }

    if (mapExists(iMap->first)) {
      // incoming map duplicates and existing name.
      if (duplicateNamesAreExceptions) 
	throw AstrometryError("PixelMapCollection::absorb with duplicate map name "
			      + iMap->first);
      // Duplicate will just be ignored.  We want to delete its
      // atom, will check later if it's needed elsewhere though.
      atomsToKill.insert(incoming.atom);
    } else {
      // A new mapName for us.  Add its mapElement to our list.
      mapElements.insert(*iMap);
    }
  }

  // Remove from the kill list any duplicate-name atomic PixelMaps that are still in use
  for (MapIter iMap = mapElements.begin();
       iMap != mapElements.end();
       ++iMap)
    if (iMap->second.atom)
      atomsToKill.erase(iMap->second.atom);
  // Then delete all the duplicaate atomic PixelMaps that are no longer in use
  for (set<PixelMap*>::iterator i = atomsToKill.begin();
       i != atomsToKill.end();
       ++i)
    delete *i;


  for (WcsIter iWcs = rhs.wcsElements.begin();
       iWcs != rhs.wcsElements.end();
       ++iWcs) {
    WcsElement& incoming = iWcs->second;
    // Kill any existing realization of incoming wcs
    if (incoming.realization) {
      delete incoming.realization;
      incoming.realization = 0;
    }

    if (wcsExists(iWcs->first)) {
      // incoming wcs duplicates and existing name.
      if (duplicateNamesAreExceptions) 
	throw AstrometryError("PixelMapCollection::absorb with duplicate wcs name "
			      + iWcs->first);
      // Duplicate will just be ignored.  
    } else {
      // A new wcsName for us.  Add its wcsElement to our list.
      wcsElements.insert(*iWcs);
    }
  }

  // And re-index everything
  rebuildParameterVector();

  // Empty out the incoming collection:
  rhs.mapElements.clear();
  rhs.wcsElements.clear();
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
    sm->wasIssued = true;	// Mark this SubMap as being owned by someone.
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

// Return a list of the atomic elements of the named map, in order of
// their application to data.  Checking for circular dependece.
list<string> 
PixelMapCollection::orderAtoms(string mapName,
			       const set<string>& ancestors) const {
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

  // Otherwise we have a compound map at work here.  First check for circularity:
  if (ancestors.find(mapName) != ancestors.end())
    throw AstrometryError("Circular dependence in PixelMapCollection::orderAtoms at map "
			  + mapName);
  // Then call recursively to all subordinateMaps, adding self to ancestor list:
  set<string> ancestorsPlusSelf(ancestors);
  ancestorsPlusSelf.insert(mapName);
  for (list<string>::const_iterator i = m.subordinateMaps.begin();
       i != m.subordinateMaps.end();
       ++i) {
    list<string> subList = orderAtoms(*i, ancestorsPlusSelf);
    retval.splice(retval.end(), subList);
  }
  return retval;
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
