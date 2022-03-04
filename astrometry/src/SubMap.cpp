// SubMap implementation: a chain of maps, which can be referred back to 
// a master parameter vector (usually held by PixelMapCollection).

// The vStartIndices, vNSubParams are assumed to have been set up by PixelMapCollection.
// If a component of this map chain has its parameters fixed, the corresponding
// vNSubParams is set to zero.

#include "PixelMapCollection.h"

using namespace astrometry;

SubMap::SubMap(const list<PixelMap*>& pixelMaps, 
	       string name,
	       bool shareMaps): PixelMap(name), 
				totalFreeParameters(0),
				ownMaps(!shareMaps),
				anyColor(false)
{
  // set up parameter vectors, making all PixelMap parameters free and consectuve by default
  for (auto& i : pixelMaps) {
    if (shareMaps) {
      vMaps.push_back(i);
    } else {
      vMaps.push_back(i->duplicate());
    }
  }
  vNSubParams.clear();
  vStartIndices.clear();
  vMapNumbers.clear();
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

PixelMap*
SubMap::duplicate() const {
  list <PixelMap*> pmlist;
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
// parameters of our components in order.
// If nSubParams[i] has been set to zero, then parameters of map i are considered fixed
// and do not appear in parameter vectors.
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

void
SubMap::toWorld(double xpix, double ypix,
		double& xworld, double& yworld,
		double color) const {
  // Make behavior with zero elements be identity matrix:
  xworld = xpix;
  yworld = ypix;
  for (int iMap = 0; iMap < nMaps(); iMap++) {
    vMaps[iMap]->toWorld(xworld, yworld, xworld, yworld, color);
  }
}

void
SubMap::toPix(double xworld, double yworld,
	      double& xpix, double& ypix,
	      double color) const {
  const int nm = nMaps();
  if (nm==1) {
    vMaps.front()->toPix(xworld, yworld, xpix, ypix, color);
    return;
  }

  // Behave as identity for empty map vector:
  if (nm==0) {
    xpix = xworld;
    ypix = yworld;
    return;
  }

  // Could consider doing the whole thing by Newton iteration ????
  // ??? A problem with current method is that tolerances on inverses
  // will build up to larger errors than perhaps desired.

  // xpix and ypix on input are supposed to be guesses at the
  // answer.  Propagate them forward to create a list of guesses
  // of "pixel" coords for all the intermediate steps.
  vector<double> xguess(nm);
  vector<double> yguess(nm);
  xguess[0] = xpix;
  yguess[0] = ypix;
  for (int iMap=0; iMap<nm-1; iMap++) {
    vMaps[iMap]->toWorld(xpix,ypix,xpix,ypix,color);
    xguess[iMap+1] = xpix;
    yguess[iMap+1] = ypix;
  }

  // Now propagate solution backwards:
  xpix = xworld;
  ypix = yworld;
  for (int iMap = nm-1; iMap>=0; iMap--) {
    xworld = xpix;
    yworld = ypix;
    xpix = xguess[iMap];
    ypix = yguess[iMap];
    vMaps[iMap]->toPix(xworld,yworld,xpix,ypix,color);
  }
}

// In derivatives vectors, parameters of map components are concatenated, skipping
// those whose params are fixed.
void 
SubMap::toWorldDerivs( double xpix, double ypix,
		       double &xworld, double &yworld,
		       DMatrix& derivs,
		       double color) const {
  if (nParams()==0) {
    toWorld(xpix, ypix, xworld, yworld, color);
    return;
  }
  int nm = nMaps();
  if (nm==1) {
    vMaps.front()->toWorldDerivs(xpix, ypix, xworld, yworld, derivs, color);
    return;
  }

  Assert(derivs.rows()==2 && derivs.cols()==nParams());
  derivs.setZero();

  // Start with identity transformation:
  xworld=xpix;
  yworld=ypix;

  int index=0;	// starting index of parameters for this map
  for (int iMap = 0; iMap < nm; iMap++) {
    // Assume coords from previous step:
    xpix = xworld;
    ypix = yworld;
    if (index>0) {
      // transform derivs of previous maps into new coords:
      DMatrix tmp = vMaps[iMap]->dWorlddPix(xpix,ypix,color) * derivs.subMatrix(0,2,0,index);
      derivs.subMatrix(0,2,0,index) = tmp; 
    }
    int nNext = nSubParams(iMap);  
    if (nNext>0) {
      Assert(nNext == vMaps[iMap]->nParams());
      DMatrix dd(2,nNext);
      vMaps[iMap]->toWorldDerivs(xpix,ypix,xworld,yworld,dd, color);
      derivs.subMatrix(0,2,index, index+nNext) = dd; 
      index += nNext;
    } else {
      // No parameters to worry about for this map:
      vMaps[iMap]->toWorld(xpix,ypix,xworld,yworld,color);
    }
  }
}

void 
SubMap::toPixDerivs( double xworld, double yworld,
		     double &xpix, double &ypix,
		     DMatrix& derivs,
		     double color) const {
  if (nParams()==0) {
    toPix(xworld, yworld, xpix, ypix, color);
    return;
  }
  if (nMaps()==1) {
    vMaps.front()->toPixDerivs(xworld, yworld, xpix, ypix, derivs, color);
    return;
  }

  // Do the inverse map then call the forward routine to find dWorld / dParams:
  toPix(xworld,yworld,xpix,ypix,color);
  toWorldDerivs(xpix, ypix, xworld, yworld, derivs, color);
}

// Calculate dWorld / dPix, we'll revert to the base class implementation of finite differences
// if we have multiple map components, guessing this will be faster than propagating
// derivatives that may or may not already be from finite differences.
Matrix22
SubMap::dWorlddPix(double xpix, double ypix, double color) const {
  if (nMaps()==0) return Matrix22().setToIdentity();
  else if (nMaps()==1) return vMaps.front()->dWorlddPix(xpix,ypix,color);
  else return PixelMap::dWorlddPix(xpix,ypix,color);
}

// Also revert to base class implementation (inverting dWorlddPix) for multiple components:
Matrix22
SubMap::dPixdWorld(double xworld, double yworld, double color) const {
  if (nMaps()==0) return Matrix22().setToIdentity();
  else if (nMaps()==1) return vMaps.front()->dPixdWorld(xworld,yworld, color);
  else return PixelMap::dPixdWorld(xworld,yworld, color);
}

void
SubMap::setPixelStep(double ps) {
  if (!vMaps.empty()) vMaps.front()->setPixelStep(ps);
}

double
SubMap::getPixelStep() const {
  return vMaps.empty() ? 1. : vMaps.front()->getPixelStep();
}

  
#ifdef USE_YAML
// Serialization (here for completeness only; current chains are
// serialized through PixelMapCollection::writeSingleMap
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
#endif
