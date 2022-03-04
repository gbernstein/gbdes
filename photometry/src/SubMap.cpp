// SubMap implementation: a chain of maps, which can be referred back to 
// a master parameter vector (usually held by PhotoMapCollection).

// The vStartIndices, vNSubParams are assumed to have been set up by PhotoMapCollection.
// If a component of this map chain has its parameters fixed, the corresponding
// vNSubParams is set to zero.

#include "PhotoMapCollection.h"

using namespace photometry;

SubMap::SubMap(const list<PhotoMap*>& photoMaps, 
	       string name,
	       bool shareMaps): PhotoMap(name), 
				totalFreeParameters(0),
				ownMaps(!shareMaps)
{
  // set up parameter vectors, making all PhotoMap parameters free and consectuve by default
  for (auto i : photoMaps) {
    if (shareMaps) {
      vMaps.push_back(i);
    } else {
      vMaps.push_back(i->duplicate());
    }
  }
  vNSubParams.clear();
  vStartIndices.clear();
  vMapNumbers.clear();
  anyColor = false;
  for (int i=0; i<nMaps(); i++) {
    vStartIndices.push_back(totalFreeParameters);
    vNSubParams.push_back(vMaps[i]->nParams());
    vMapNumbers.push_back(i);
    totalFreeParameters += vMaps[i]->nParams();
    if (vMaps[i]->needsColor()) anyColor = true;
  }
}

SubMap::~SubMap() {
  if (ownMaps) 
    for (auto i : vMaps)
      delete i;
}

PhotoMap*
SubMap::duplicate() const {
  list <PhotoMap*> pmlist;
  for (auto i : vMaps)
    pmlist.push_back(i);
  return new SubMap(pmlist, getName(), !ownMaps);
}

void
SubMap::countFreeParameters() {
  totalFreeParameters = 0;
  for (int i=0; i<vNSubParams.size(); i++)
    totalFreeParameters += vNSubParams[i];
}

// When get/setParams are called from here directly, we just concatenate the free
// parameters of our components in order:
DVector
SubMap::getParams() const {
  if (nParams()==0) return DVector(0);
  int nm = nMaps();
  if (nm==1) return vMaps.front()->getParams();
  DVector p(nParams(),0.);
  int outIndex = 0;
  for (int iMap = 0; iMap < nm; ++iMap) {
    int nsub = nSubParams(iMap);
    if (nsub==0) continue;
    Assert(nsub == vMaps[iMap]->nParams());
    DVector pp = vMaps[iMap]->getParams();
    for (int j=0; j<nsub; j++, outIndex++) p[outIndex]=pp[j];
  }
  return p;
}

void
SubMap::setParams(const DVector& p) {
  if (nParams()==0) return;
  int nm = nMaps();
  if (nm==1) vMaps.front()->setParams(p);
  Assert(p.size() == nParams());
  int index = 0;
  for (int iMap = 0; iMap < nm; iMap++) {
    int nsub = nSubParams(iMap);
    if (nsub==0) continue;
    Assert(nsub == vMaps[iMap]->nParams());
    DVector pp(nsub);
    for (int j=0; j<nsub; j++, index++) pp[j]=p[index];
    vMaps[iMap]->setParams(pp);
  }
}

double
SubMap::forward(double magIn, const PhotoArguments& args) const {
  // Make behavior with zero elements be identity matrix:
  double magOut = magIn;
  for (int iMap = 0; iMap < nMaps(); iMap++) {
    magIn = magOut;
    magOut = vMaps[iMap]->forward(magIn, args);
  }
  return magOut;
}

double
SubMap::inverse(double magOut, const PhotoArguments& args) const {
  // Make behavior with zero elements be identity matrix:
  double magIn = magOut;
  for (int iMap = nMaps()-1; iMap >= 0; iMap--) {
    magOut = magIn;
    magIn = vMaps[iMap]->inverse(magOut, args);
  }
  return magIn;
}


// In derivatives vectors, parameters of map components are concatenated, skipping
// those whose params are fixed:
double
SubMap::forwardDerivs( double magIn, const PhotoArguments& args,
		       DVector& derivs) const {
  Assert(derivs.size()==nParams());
  if (nParams()==0) 
    return forward(magIn, args);
  int nm = nMaps();
  if (nm==1) 
    return vMaps.front()->forwardDerivs(magIn, args, derivs);

  derivs.setZero();

  // Start with identity transformation:
  double magOut = magIn;

  int index=0;	// starting index of parameters for this map
  for (int iMap = 0; iMap < nm; iMap++) {
    magIn = magOut;
    if (index>0) {
      // transform derivs of previous maps into new coords:
      derivs.subVector(0,index) *= vMaps[iMap]->derivative(magIn, args);
    }
    int nNext = nSubParams(iMap);  
    if (nNext>0) {
      Assert(nNext == vMaps[iMap]->nParams());
      DVector dd(nNext);
      magOut = vMaps[iMap]->forwardDerivs(magIn, args, dd);
      derivs.subVector(index, index+nNext) += dd; 
      index += nNext;
    } else {
      // No parameters to worry about for this map:
      magOut = vMaps[iMap]->forward(magIn, args);
    }
  }
  return magOut;
}

// Calculate dOut / dIn, with chain rule:
double
SubMap::derivative(double magIn, const PhotoArguments& args) const {
  if (nMaps()==0) return 1.;
  double deriv = 1.;
  double magOut = magIn;
  // Chain all but the last map together:
  for (int iMap = 0; iMap < nMaps()-1; iMap++) {
    magIn = magOut;
    magOut = vMaps[iMap]->forward(magIn, args);
    deriv *= vMaps[iMap]->derivative(magIn, args);
  }
  magIn = magOut;
  deriv *= vMaps.back()->derivative(magIn, args);
  return deriv;
}

// Serialization (here for completeness only; current chains are
// serialized through PhotoMapCollection::writeSingleMap
void
SubMap::write(YAML::Emitter& os) const {
  os << YAML::BeginMap
     << YAML::Key << "Type" << YAML::Value << SubMap::type()
     << YAML::Key << "Elements"
     << YAML::BeginSeq;
  for (auto i=vMaps.begin(); i != vMaps.end(); ++i) {
    os << (*i)->getName();
  }
  os << YAML::EndSeq;
  os << YAML::EndMap;
  return;
}
