// Define all the PhotoMapCollection routines.

#include "PhotoMapCollection.h"
#include "StringStuff.h"

using namespace photometry;

// The word that we expect all serialized collections to have on their first line
const string magicHeader = "PhotoMapCollection";

//////////////////////////////////////////////////////////////
// Static data
//////////////////////////////////////////////////////////////
vector<string>
PhotoMapCollection::mapTypeNames;

vector<PhotoMapCollection::Creator>
PhotoMapCollection::mapCreators;

bool
PhotoMapCollection::creatorsInitialized=false;

void
PhotoMapCollection::PhotoMapTypeInitialize() {
  if (creatorsInitialized) return;
  creatorsInitialized = true;
  registerMapType<IdentityMap>();
  registerMapType<PolyMap>();
}

//////////////////////////////////////////////////////////////
// Construction / destruction
//////////////////////////////////////////////////////////////

PhotoMapCollection::PhotoMapCollection(): parameterCount(0) {
  PhotoMapTypeInitialize();
}

PhotoMapCollection::~PhotoMapCollection() {
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
}

//////////////////////////////////////////////////////////////
// Parameter manipluation
//////////////////////////////////////////////////////////////

void
PhotoMapCollection::setParams(const DVector& p) {
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
PhotoMapCollection::getParams() const {
  DVector p(parameterCount, 888.);
  for (ConstMapIter i = mapElements.begin(); i!=mapElements.end(); ++i) {
    const MapElement& map = i->second;
    if (map.isFixed) continue;
    int nSub = map.nParams;
    if (nSub<=0) continue;
    if (!map.atom) cerr << "mapElement is not atomic: " << i->first << " params: " << nSub << endl;
    Assert(map.atom);
    p.subVector(map.startIndex, map.startIndex+nSub) = 
      map.atom->getParams().subVector(0,nSub);
  }
  return p;
}

void 
PhotoMapCollection::setFixed(list<string> nameList, bool isFixed) {
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
PhotoMapCollection::getFixed(string name) const {
  ConstMapIter i = mapElements.find(name);
  if (i==mapElements.end())
    throw PhotometryError("No member of PhotoMapCollection named " + name);
  return i->second.isFixed;
}

void
PhotoMapCollection::rebuildParameterVector() {
  // First assign new starting indices to every atomic map element

  // Restart the parameter counting:
  parameterCount = 0;
  // And map counting
  atomCount = 0;
  freeCount = 0;
  for (MapIter i = mapElements.begin(); i!=mapElements.end(); ++i) {
    MapElement& map = i->second;
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
	throw PhotometryError("PhotoMapCollection::rebuildParameterVector could not find"
			      " map element with name " + elementName);
      if (!j->second.atom)
	throw PhotometryError("PhotoMapCollection element " + elementName + " is not atomic"
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
PhotoMapCollection::allMapNames() const {
  vector<string> output;
  for (ConstMapIter i = mapElements.begin();
       i != mapElements.end();
       ++i) 
    output.push_back(i->first);
  return output;
}

void 
PhotoMapCollection::learnMap(const PhotoMap& pm, bool duplicateNamesAreExceptions) {
  if (mapExists(pm.getName())) {
    if (duplicateNamesAreExceptions)
      throw PhotometryError("Duplicate map name in PhotoMapCollection::learnMap at "
			    + pm.getName());
    // If not throwing an exception, we will just ignore this duplicate-named map.
    return;
  }
  const SubMap* sm = dynamic_cast<const SubMap*> (&pm);
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
	throw PhotometryError("PhotoMapCollection::learnMap encountered SubMap that has "
			      "a dependence on PhotoMap with the same name: " + sm->getName());
      mel.subordinateMaps.push_back(sm->vMaps[i]->getName());
      learnMap(*sm->vMaps[i], duplicateNamesAreExceptions);
    }
    // Register this compound map
    mapElements.insert(std::pair<string,MapElement>(pm.getName(), mel));

  } else {
    //  This should be an atomic PhotoMap; simply add it to our element list:
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
PhotoMapCollection::learn(PhotoMapCollection& rhs, bool duplicateNamesAreExceptions) {
  for (ConstMapIter iMap = rhs.mapElements.begin();
       iMap != rhs.mapElements.end();
       ++iMap) {
    const MapElement& incoming = iMap->second;
    if (mapExists(iMap->first)) {
      // incoming map duplicates and existing name.
      if (duplicateNamesAreExceptions) 
	throw PhotometryError("learn(PhotoMapCollection) with duplicate map name "
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

  // And re-index everything
  rebuildParameterVector();
}

// Define a new PhotoMap that is compounding of a list of other PhotoMaps.  Order
// of the list is order of transform from mapIn to mapOut
void 
PhotoMapCollection::defineChain(string chainName, const list<string>& elements) {
  if (mapExists(chainName)) 
    throw PhotometryError("PhotoMapCollection::defineChain with duplicate name: " 
			  + chainName);
  // Check that elements exist
  for (list<string>::const_iterator i = elements.begin();
       i != elements.end();
       ++i) 
    if (!mapExists(*i))
      throw PhotometryError("PhotoMapCorrection::defineChain with unknown photo map element: "
			    + *i);
  mapElements.insert(std::pair<string, MapElement>(chainName, MapElement()));
  mapElements[chainName].subordinateMaps = elements;
  // Note that adding a new chain does not change parameter vector assignments
  // Nor can it introduce circular dependence if chain name is new and all elements exist.
}

// Return pointer to a SubMap realizing the named magnitude transformation
SubMap* 
PhotoMapCollection::issueMap(string mapName) {
  if (!mapExists(mapName))
    throw PhotometryError("PhotoMapCollection::issueMap requested for unknown PhotoMap: "
			  +  mapName);
  MapElement& el = mapElements[mapName];

  if (!el.realization) {
    // Create a realization if one does not exist
    list<string> atomList = orderAtoms(mapName);
    vector<int> startIndices(atomList.size(), 0);
    vector<int> nParams(atomList.size(), 0);
    vector<int> mapNumbers(atomList.size(), -1);

    list<PhotoMap*> atoms;
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
    SubMap* sm = new SubMap(atoms, mapName, true);
    sm->vStartIndices = startIndices;
    sm->vNSubParams = nParams;
    sm->vMapNumbers = mapNumbers;
    sm->countFreeParameters();
    el.realization = sm;
  }
  return el.realization;
}

PhotoMap*
PhotoMapCollection::cloneMap(string mapName) const {
  if (!mapExists(mapName))
    throw PhotometryError("PhotoMapCollection::issueMap requested for unknown PhotoMap: "
			  +  mapName);
  const MapElement& el = mapElements.find(mapName)->second;

  if (el.atom) 
    return el.atom->duplicate();

  // A composite one we will create as a SubMap:
  list<PhotoMap*> vMaps;
  for (list<string>::const_iterator i = el.subordinateMaps.begin();
	 i != el.subordinateMaps.end();
	 ++i) 
    vMaps.push_back(cloneMap(*i));
  return new SubMap(vMaps, mapName, false);
}

set<string>
PhotoMapCollection::dependencies(string target) const {
  // Find target MapElement and add target to output dependency list.  
  // Assume no circular dependence.
  ConstMapIter iTarget = mapElements.find(target);
  if (iTarget == mapElements.end()) 
    throw PhotometryError("PhotoMapCollection has no element " + target + " in dependencies()");
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
PhotoMapCollection::orderAtoms(string mapName) const {
  if (!mapExists(mapName))
    throw PhotometryError("PhotoMapCollection::orderAtoms requested for unknown map: "
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
PhotoMapCollection::checkCircularDependence(string mapName,
					    const set<string>& ancestors) const {
  if (!mapExists(mapName))
    throw PhotometryError("Unknown mapName in PhotoMapCollection::checkCircularDependency: "
			  + mapName);
  if (ancestors.count(mapName))
    throw PhotometryError("Circular dependency in PhotoMapCollection at map "
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

void
PhotoMapCollection::checkCompleteness() const {
  // Loop through all the (non-atomic) maps, throw exception if the refer to 
  // non-existent map elements or to themselves
  for (ConstMapIter i = mapElements.begin();
       i != mapElements.end();
       ++i)
    checkCircularDependence(i->first);
}


//////////////////////////////////////////////////////////////
// (De-) Serialization
//////////////////////////////////////////////////////////////

void
PhotoMapCollection::writeSingleMap(ostream& os, const MapElement& mel, string name)  const {
  if (mel.atom) {
    // Atomic map, give all details:
    os << mel.atom->getType() 
       << " " << name
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
PhotoMapCollection::write(ostream& os) const {
  os << magicHeader << endl;
  // Write the map elements
  for (ConstMapIter i = mapElements.begin();
       i != mapElements.end();
       ++i) {
    os << "#" << endl;
    writeSingleMap(os, i->second, i->first);
  }
}

void 
PhotoMapCollection::writeMap(ostream& os, string name) const {
  if (!mapExists(name))
    throw PhotometryError("PhotoMapCollection::writeMap() for unknown map: " + name);
  set<string> allmaps = dependencies(name);
  os << magicHeader << endl;
  // Write the map elements
  for (set<string>::const_iterator i = allmaps.begin();
       i != allmaps.end();
       ++i) {
    ConstMapIter j = mapElements.find(*i);
    writeSingleMap(os, j->second, *i);
  }
}

bool
PhotoMapCollection::read(istream& is, string namePrefix) {
  string buffer;
  // Check for magic word
  if (!stringstuff::getlineNoComment(is, buffer)) 
    throw PhotometryError("PhotoMapCollection::read() found no information for " + namePrefix);
  istringstream iss(buffer);
  {
    string word;
    iss >> word;
    if (!stringstuff::nocaseEqual(word, magicHeader))
      return false;
  }
  while (stringstuff::getlineNoComment(is, buffer)) {
    iss.str(buffer);
    iss.clear();
    // Read first word to get type of next element and its name
    string otype;
    string name;
    if (!(iss >> otype >> name)) 
      throw PhotometryError("PhotoMapCollection::read() format error at line: " + buffer);
    name = namePrefix + name;

    if (stringstuff::nocaseEqual(otype, SubMap::mapType())) {
      int nmaps;
      if (!(iss >> nmaps)) 
	throw PhotometryError("PhotoMapCollection::read() missing number of SubMap components at line "
			      + buffer);
      // Read names of components
      list<string> submaps;
      while (submaps.size() < nmaps) {
	if (!(stringstuff::getlineNoComment(is, buffer)))
	  throw PhotometryError("PhotoMapCollection::read() out of input for SubMap: " + name);
	string n;
	iss.str(buffer);
	iss.clear();
	do {
	  iss >> n;
	  if (!iss.fail()) {
	    submaps.push_back(n);
	  } else if (iss.eof()) {
	    break;
	  } else {
	    // Failure other than eof:
	    throw PhotometryError("PhotoMapCollection::read() error getting SubMaps at <" 
				  + buffer + ">");
	  }
	} while (submaps.size() < nmaps);
      }

      // Create the MapElement
      MapElement mel;
      mel.subordinateMaps = submaps;
      mapElements.insert(std::pair<string,MapElement>(name, mel));

    } else {
      // Access the static arrays of map type names and creators
      Assert (mapTypeNames.size() == mapCreators.size());
      PhotoMap* atom = 0;
      // Loop through registered map types
      for (int iType=0; iType < mapTypeNames.size(); iType++) {
	if (stringstuff::nocaseEqual(otype, mapTypeNames[iType])) {
	  // call the create() for the one that matches
	  atom = (*mapCreators[iType])(is, name);
	  break;
	}
      }
      // throw for unknown map type
      if (!atom)
	throw PhotometryError("PhotoMapCollection does not recognize PhotoMap type " + otype);

      // create the MapElement
      MapElement mel;
      mel.atom = atom;
      mapElements.insert(std::pair<string,MapElement> (name, mel));
    }
  } // end input line loop

  // Check all maps for circular dependence and completeness.
  checkCompleteness();

  // Recalculate indices
  rebuildParameterVector();

  return true;
}

