// Define all the PhotoMapCollection routines.

#include "PhotoMapCollection.h"
#include "StringStuff.h"

using namespace photometry;

// The word that we expect all serialized collections to have on their first line
const string
PhotoMapCollection::magicKey = "PhotoMapCollection";

//////////////////////////////////////////////////////////////
// Static data
//////////////////////////////////////////////////////////////
map<string,PhotoMapCollection::Creator>
PhotoMapCollection::mapCreators;

bool
PhotoMapCollection::creatorsInitialized=false;

void
PhotoMapCollection::PhotoMapTypeInitialize() {
  if (creatorsInitialized) return;
  creatorsInitialized = true;
  registerMapType<IdentityMap>();
  registerMapType<ConstantMap>();
  registerMapType<PolyMap>();
}

//////////////////////////////////////////////////////////////
// Construction / destruction
//////////////////////////////////////////////////////////////

PhotoMapCollection::PhotoMapCollection() {
  PhotoMapTypeInitialize();
  rebuildParameterVector();
}

PhotoMapCollection::~PhotoMapCollection() {
  // Destroy all of the compounded maps and submaps we've made:
  for (auto& i : mapElements) {
      delete i.second.atom;
      delete i.second.realization;
  }
}

//////////////////////////////////////////////////////////////
// Parameter manipluation
//////////////////////////////////////////////////////////////

void
PhotoMapCollection::setParams(const DVector& p) {
  Assert(p.size()==parameterCount);
  for (auto& i : mapElements) {
    MapElement& map = i.second;
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
  DVector p(parameterCount, NODATA);
  for (auto& melpair : mapElements) {
    const MapElement& map = melpair.second;
    if (map.isFixed) continue;
    int nSub = map.nParams;
    if (nSub<=0) continue;
    if (!map.atom) 
      cerr << "mapElement is not atomic: " << melpair.first << " params: " << nSub << endl;
    Assert(map.atom);
    p.subVector(map.startIndex, map.startIndex+nSub) = map.atom->getParams();
  }
  return p;
}

void
PhotoMapCollection::copyParamsFrom(const PhotoMap& pm) {
  auto melpair = mapElements.find(pm.getName());
  if (melpair == mapElements.end())
    throw PhotometryError("Attempt to copyParamsFrom non-existent map element " + pm.getName());
  if (!melpair->second.atom)
    throw PhotometryError("Attempt to copyParamsFrom non-atomic map element " + pm.getName());
  melpair->second.atom->setParams(pm.getParams());
}

void 
PhotoMapCollection::setFixed(string name, bool isFixed) {
  set<string> s;
  s.insert(name);
  setFixed(s, isFixed);
}

void 
PhotoMapCollection::setFixed(set<string> nameList, bool isFixed) {
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
PhotoMapCollection::getFixed(string name) const {
  auto i = mapElements.find(name);
  if (i==mapElements.end())
    throw PhotometryError("No member of PhotoMapCollection named " + name);
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
PhotoMapCollection::isAtomic(string name) const {
  auto i = mapElements.find(name);
  if (i==mapElements.end())
    throw PhotometryError("No member of PhotoMapCollection named " + name);
  if (i->second.atom) {
    return true;
  } else {
    return false;
  }
}

void
PhotoMapCollection::rebuildParameterVector() {
  // First assign new starting indices to every atomic map element

  // Restart the parameter counting:
  parameterCount = 0;
  // And map counting
  atomCount = 0;
  freeCount = 0;
  for (auto& i : mapElements) {
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
  for (auto& i : mapElements)
    output.push_back(i.first);
  return output;
}

void 
PhotoMapCollection::learnMap(const PhotoMap& pm,
			     bool duplicateNamesAreExceptions,
			     bool rebuildIndices) {
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
      learnMap(*sm->vMaps.front(), duplicateNamesAreExceptions, rebuildIndices);
      return;
    } 
    // create a new compound map from this, learning each dependency:
    MapElement mel;
    for (auto mptr : sm->vMaps) {
      if (sm->getName() == mptr->getName())
	throw PhotometryError("PhotoMapCollection::learnMap encountered SubMap that has "
			      "a dependence on PhotoMap with the same name: " + sm->getName());
      mel.subordinateMaps.push_back(mptr->getName());
      learnMap(*mptr, duplicateNamesAreExceptions, rebuildIndices);
    }
    // Register this compound map
    mapElements.insert(std::pair<string,MapElement>(pm.getName(), mel));

  } else {
    //  This should be an atomic PhotoMap; simply add it to our element list:
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
PhotoMapCollection::learn(PhotoMapCollection& rhs, bool duplicateNamesAreExceptions) {
  if (&rhs == this) return;	// No need to learn self
  for (auto& iMap : rhs.mapElements) {
    const MapElement& incoming = iMap.second;
    if (mapExists(iMap.first)) {
      // incoming map duplicates and existing name.
      if (duplicateNamesAreExceptions) 
	throw PhotometryError("learn(PhotoMapCollection) with duplicate map name "
			      + iMap.first);
      // Duplicate will just be ignored. 
    } else {
      // A new mapName for us.  Add its mapElement to our list.
      MapElement mel;
      if (iMap.second.atom) mel.atom = iMap.second.atom->duplicate();
      mel.isFixed = iMap.second.isFixed;
      mel.subordinateMaps = iMap.second.subordinateMaps;
      mapElements.insert(std::pair<string,MapElement>(iMap.first, mel));
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
  for (auto& i : elements)
    if (!mapExists(i))
      throw PhotometryError("PhotoMapCorrection::defineChain with unknown photo map element: "
			    + i);
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
  for (auto& i : el.subordinateMaps)
    vMaps.push_back(cloneMap(i));
  SubMap* retval = new SubMap(vMaps, mapName, false);
  // Clean up the stray clones since they've been duplicated by SubMap:
  for (auto i : vMaps)
    delete i;
  return retval;
}

set<string>
PhotoMapCollection::dependencies(string target) const {
  // Find target MapElement and add target to output dependency list.  
  // Assume no circular dependence.
  auto iTarget = mapElements.find(target);
  if (iTarget == mapElements.end()) 
    throw PhotometryError("PhotoMapCollection has no element " + target + " in dependencies()");
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
PhotoMapCollection::dependsOn(const string mapName, const string targetName) const {
  auto melptr = mapElements.find(mapName);
  if (melptr == mapElements.end()) 
    throw PhotometryError("PhotoMapCollection has no element " + mapName + " in dependsOn()");

  if (mapName==targetName) return true;  // True if map is the target

  // Call this routine recursively for all the map elements that the target depends on
  for (auto& i : melptr->second.subordinateMaps)
    if (dependsOn(i, targetName)) return true;

  return false;
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
  for (auto& i : m.subordinateMaps) {
    list<string> subList = orderAtoms(i);
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
  for (auto& i : m.subordinateMaps) 
    checkCircularDependence(i, ancestorsPlusSelf);
}

void
PhotoMapCollection::checkCompleteness() const {
  // Loop through all the (non-atomic) maps, throw exception if the refer to 
  // non-existent map elements or to themselves
  for (auto& i : mapElements) 
    checkCircularDependence(i.first);
}


void
PhotoMapCollection::invalidate(string mapName) {
  auto it = mapElements.find(mapName);
  if (it==mapElements.end()) return;
  it->second.isValid = false;
}

template <class T>
void
set_subtract(std::set<T>& s1, const std::set<T> s2) {
  auto i1 = s1.begin();
  auto i2 = s2.begin();
  while(i1!=s1.end() && i2!=s2.end()) {
    if (*i1 < *i2) {
      i1++;
    } else if (*i2 < *i1) {
      i2++;
    } else {
      // Equality: remove from s1
      i1 = s1.erase(i1);
      i2++;
    }
  }
}
				       
void
PhotoMapCollection::purgeInvalid() {
  // Collect names of all maps that invalid maps use,
  // and get rid of them.
  set<string> unneeded;
  for (auto& i : mapElements) {
    if (i.second.isValid) continue;
    auto depend = dependencies(i.first);
    unneeded.insert(depend.begin(), depend.end());
  }
  // Now go through all PhotoMaps.  If they are
  // not dependencies of the bad ones, then
  // we keep everything that they need.
  for (auto& i : mapElements) {
    if (unneeded.count(i.first)) continue; // It's a dependency of bad ones.
    // If not, then it's part of a good one,
    // and we want to keep everything it needs.
    set_subtract(unneeded, dependencies(i.first));
  }
  // Now kill all unneeded maps
  for (auto s : unneeded)
    removeMap(s);

  // Recalculate all parameter indices
  rebuildParameterVector();
}

void
PhotoMapCollection::removeMap(string mapName) {
  auto it = mapElements.find(mapName);
  if (it == mapElements.end())
    return; // Do nothing if there is no such map
  if (it->second.atom) {
    delete it->second.atom;
    it->second.atom = nullptr;
  }
  if (it->second.realization) {
    delete it->second.realization;
    it->second.realization = nullptr;
  }
  mapElements.erase(it);
}
//////////////////////////////////////////////////////////////
// YAML (De-) Serialization
//////////////////////////////////////////////////////////////

void
PhotoMapCollection::writeSingleMap(YAML::Emitter& os,
				   const MapElement& mel, string name)  const {
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
       << YAML::Key << "Type"
       << YAML::Value << SubMap::type()
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
PhotoMapCollection::write(ostream& os, string comment) const {
  // Create a YAML emitter.  Produce a map where one
  // key is the magic word, and the rest are names of maps.
  YAML::Emitter out;
  out << YAML::BeginMap;
  if (comment.size()>0)
    out << YAML::Comment(comment);
  out << YAML::Key << magicKey
      << YAML::Value << "This is a serialized PhotoMapCollection";
  for (auto& melpair : mapElements) {
    writeSingleMap(out, melpair.second, melpair.first);
  }
  out << YAML::EndMap;
  // Send the YAML to the output stream
  os << out.c_str() << endl;
}

void 
PhotoMapCollection::writeMap(ostream& os, string name, string comment) const {
  if (!mapExists(name))
    throw PhotometryError("PhotoMapCollection::writeMap() for unknown map: " + name);
  set<string> allmaps = dependencies(name);
  YAML::Emitter out;
  out << YAML::BeginMap;
  if (comment.size()>0)
    out << YAML::Comment(comment);
  out << YAML::Key << magicKey
      << YAML::Value << "This is a serialized PhotoMapCollection";
  // Write the map elements
  for (auto mapname : allmaps) {
    auto j = mapElements.find(mapname);
    writeSingleMap(out, j->second, mapname);
  }
  out << YAML::EndMap;
  // Send the YAML to the output stream
  os << out.c_str() << endl;
}

bool
PhotoMapCollection::read(istream& is, string namePrefix) {
  string buffer;
  try {
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

      if (!node.IsMap() || !node["Type"]) 
	throw PhotometryError("Non-map or missing Type key in PhotoMapCollection YAML entry " +
			      name);
      string mapType = node["Type"].as<string>();
      if (stringstuff::nocaseEqual(mapType,ColorTerm::type())) {
	// A color term has the function it modifies in a child node:
	if (!node["Function"])
	throw PhotometryError("Missing Function key in ColorTerm YAML entry " +
			      name);
	double reference = 0.; // Default reference color
	if (node["Reference"])
	  reference = node["Reference"].as<double>();
	PhotoMap* pm = createAtomFromNode(node["Function"], "");

	// create the MapElement
	MapElement mel;
	mel.atom = new ColorTerm(pm, reference, name);
	mapElements.insert(std::pair<string,MapElement> (name, mel));
	
      } else if (stringstuff::nocaseEqual(mapType,SubMap::type())) {
	// Composite maps: all names stored in a child node
	if (!node["Elements"])
	  throw PhotometryError("Missing Elements key in SubMap YAML entry " + name);
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
	mel.atom = createAtomFromNode(node,name);
	mel.isFixed = mel.atom->nParams()==0;  // Set fixed if no free params
	mapElements.insert(std::pair<string,MapElement> (name, mel));
      }

    } // end loop through root node map members

  } catch (YAML::Exception& e) {
    /**/cerr << "PhotoMapCollection::read() caught: " << e.what() << endl;
    // Return false if not a valid YAML file
    return false;
  }

  // Check all maps for circular dependence and completeness.
  checkCompleteness();

  // Recalculate indices
  rebuildParameterVector();

  return true;
}

PhotoMap*
PhotoMapCollection::createAtomFromNode(const YAML::Node& node, string name) const {
  // Access the static arrays of map type names and creators
  if (!node["Type"])
    throw PhotometryError("Missing Type key in createAtomFromNode at " + name);

  string mapType = node["Type"].as<string>();
  auto it = mapCreators.find(mapType);
  if (it==mapCreators.end())
    throw PhotometryError("PhotoMapCollection does not recognize PhotoMap Type " + mapType);
  return (*(it->second))(node, name);
}

//////////////////////////////////////////////////////////////////////////////

string
PhotoMapCollection::atomHavingParameter(int parameterIndex) const {
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
PhotoMapCollection::parameterIndicesOf(string mapname,
				       int& startIndex,
				       int& nParams) const {
  if (!mapExists(mapname))
    throw PhotometryError("PhotoMapCollection::parameterIndicesof() unknown map: "
			  + mapname);
  startIndex = 0;
  nParams = 0;
  auto& mel = mapElements.at(mapname);
  if (!mel.atom || mel.isFixed) return;
  startIndex = mel.startIndex;
  nParams = mel.nParams;
  return;
}
