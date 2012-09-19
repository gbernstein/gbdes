// $Id: PixelMapCollection.cpp,v 1.1 2012/02/02 23:34:25 garyb Exp $
#include "PixelMapCollection.h"
using namespace astrometry;

void
PixelMapCollection::setParams(const DVector& p) {
  Assert(p.size()==parameterCount);
  for (int i=0; i<members.size(); i++) {
    int nSub = members[i]->nParams();
    if (nSub<=0) continue;
    DVector subp(nSub, 0.);
    subp = p.subVector(startIndices[i],
		       startIndices[i]+nSub);
    members[i]->setParams(subp);
  }
}

DVector
PixelMapCollection::getParams() const {
  DVector p(parameterCount, 888.);
  for (int i=0; i<members.size(); i++) {
    int nSub = members[i]->nParams();
    if (nSub<=0) continue;
    p.subVector(startIndices[i],
    		startIndices[i]+nSub) = members[i]->getParams().subVector(0, nSub);
  }
  return p;
}

PixelMapCollection::~PixelMapCollection() {
  // Destroy all of the compounded maps and submaps we've made:
  for (list<SubMap*>::iterator i = createdSubMaps.begin();
       i != createdSubMaps.end();
       ++i) {
    if (*i) delete *i;
    *i = 0;
  }
  for (list<CompoundPixelMap*>::iterator i = createdCompoundMaps.begin();
       i != createdCompoundMaps.end();
       ++i) {
    if (*i) delete *i;
    *i = 0;
  }
  // And all of the map elements contributed to the collection:
  for (vector<PixelMap*>::iterator i = members.begin();
       i != members.end();
       ++i) {
    if (*i) delete *i;
    *i = 0;
  }
}

PixelMapKey
PixelMapCollection::add(PixelMap* pm, string name) {
  members.push_back(pm);
  startIndices.push_back(parameterCount);
  parameterCount += pm->nParams();
  PixelMapKey m = members.size()-1;
  string useName=name;
  if (useName.empty()) {
    ostringstream oss;
    oss << "Element_" << m ;
    useName = oss.str();
  }
  if (names.has(useName))
    FormatAndThrow<AstrometryError>() << "Duplicate PixelMapCollection name: " << useName;
  PixelMapKey m2 = names.append(useName);
  Assert(m2 == m);
  return m;
}

SubMap* 
PixelMapCollection::issue(PixelMapChain& chain) {
  Assert(!chain.empty());
  // Start compound map
  PixelMapChain::iterator i = chain.begin();
  PixelMapKey index = *i;
  CompoundPixelMap* cpm = new CompoundPixelMap(members[index]);
  // And SubMap
  SubMap* sm = new SubMap(cpm);
  sm->startIndices.push_back(startIndices[index]);
  sm->nSubParams.push_back(members[index]->nParams());

  // add the rest of the elements of the map chain
  for (++i; i!=chain.end(); ++i) {
    index = *i;
    cpm->append(members[index]);
    sm->startIndices.push_back(startIndices[index]);
    sm->nSubParams.push_back(members[index]->nParams());
  }
  sm->chain = chain;

  // Record these and return the SubMap pointer
  createdCompoundMaps.push_back(cpm);
  createdSubMaps.push_back(sm);
  return sm;
}

