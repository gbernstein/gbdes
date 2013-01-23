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
				ownMaps(!shareMaps)
{
  // set up parameter vectors, making all PixelMap parameters free and consectuve by default
  for (list<PixelMap*>::const_iterator i = pixelMaps.begin();
       i != pixelMaps.end();
       ++i) {
    if (shareMaps) {
      vMaps.push_back(*i);
    } else {
      vMaps.push_back((*i)->duplicate());
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
  }
}

SubMap::~SubMap() {
  if (ownMaps) 
    for (int i=0; i<vMaps.size(); i++)
      delete vMaps[i];
}

PixelMap*
SubMap::duplicate() const {
  list <PixelMap*> pmlist;
  for (int i=0; i<nMaps(); i++)
    pmlist.push_back(vMaps[i]);
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

void
SubMap::toWorld(double xpix, double ypix,
		double& xworld, double& yworld) const {
  // Make behavior with zero elements be identity matrix:
  xworld = xpix;
  yworld = ypix;
  for (int iMap = 0; iMap < nMaps(); iMap++) {
    vMaps[iMap]->toWorld(xworld, yworld, xworld, yworld);
  }
}

void
SubMap::toPix(double xworld, double yworld,
	      double& xpix, double& ypix) const {
  const int nm = nMaps();
  if (nm==1) {
    vMaps.front()->toPix(xworld, yworld, xpix, ypix);
    return;
  }

  // Behave as identity for empty map vector:
  if (nm==0) return;

  // xpix and ypix on input are supposed to be guesses at the
  // answer.  Propagate them forward to create a list of guesses
  // of "pixel" coords for all the intermediate steps.
  vector<double> xguess(nm);
  vector<double> yguess(nm);
  for (int iMap=0; iMap<nm-1; iMap++) {
    xguess[iMap] = xpix;
    yguess[iMap] = ypix;
    vMaps[iMap]->toWorld(xpix,ypix,xpix,ypix);
  }

  xpix = xworld;
  ypix = yworld;

  // Now propagate solution backwards:
  for (int iMap = nm-1; iMap>=0; iMap--) {
    xworld = xpix;
    yworld = ypix;
    xpix = xguess[iMap];
    ypix = yguess[iMap];
    vMaps[iMap]->toPix(xworld,yworld,xpix,ypix);
    cerr << iMap << " " << xworld << "," << yworld << " " << xpix << "," << ypix << endl;
  }
}

// In derivatives vectors, parameters of map components are concatenated, skipping
// those whose params are fixed:
void 
SubMap::toWorldDerivs( double xpix, double ypix,
		       double &xworld, double &yworld,
		       DMatrix& derivs) const {
  if (nParams()==0) {
    toWorld(xpix, ypix, xworld, yworld);
    return;
  }
  int nm = nMaps();
  if (nm==1) vMaps.front()->toWorldDerivs(xpix, ypix, xworld, yworld, derivs);

  Assert(derivs.nrows()==2 && derivs.ncols()==nParams());
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
      // transform derives of previous maps into new coords:
      DMatrix tmp = vMaps[iMap]->dWorlddPix(xpix,ypix) * derivs.colRange(0,index);
      derivs.colRange(0,index) = tmp; 
    }
    int nNext = nSubParams(iMap);  
    if (nNext>0) {
      Assert(nNext == vMaps[iMap]->nParams());
      DMatrix dd(2,nNext);
      vMaps[iMap]->toWorldDerivs(xpix,ypix,xworld,yworld,dd);
      derivs.colRange(index, index+nNext) += dd;  // Could just set =, not += here??
      index += nNext;
    } else {
      // No parameters to worry about for this map:
      vMaps[iMap]->toWorld(xpix,ypix,xworld,yworld);
    }
  }
}

void 
SubMap::toPixDerivs( double xworld, double yworld,
		     double &xpix, double &ypix,
		     DMatrix& derivs) const {
  if (nParams()==0) {
    toPix(xworld, yworld, xpix, ypix);
    return;
  }
  if (nMaps()==1) vMaps.front()->toPixDerivs(xworld, yworld, xpix, ypix, derivs);

  // Do the inverse map then call the forward routine to find dWorld / dParams:
  toPix(xworld,yworld,xpix,ypix);
  toWorldDerivs(xpix, ypix, xworld, yworld, derivs);
}

// Calculate dWorld / dPix, we'll revert to the base class implementation of finite differences
// if we have multiple map components, guessing this will be faster than propagating
// derivatives that may or may not already be from finite differences.
Matrix22
SubMap::dWorlddPix(double xpix, double ypix) const {
  if (nMaps()==0) return Matrix22().setToIdentity();
  else if (nMaps()==1) return vMaps.front()->dWorlddPix(xpix,ypix);
  else return PixelMap::dWorlddPix(xpix,ypix);
}

// Also revert to base class implementation (inverting dWorlddPix) for multiple components:
Matrix22
SubMap::dPixdWorld(double xworld, double yworld) const {
  if (nMaps()==0) return Matrix22().setToIdentity();
  else if (nMaps()==1) return vMaps.front()->dPixdWorld(xworld,yworld);
  else return PixelMap::dPixdWorld(xworld,yworld);
}

void
SubMap::setPixelStep(double ps) {
  if (!vMaps.empty()) vMaps.front()->setPixelStep(ps);
}

double
SubMap::getPixelStep() const {
  return vMaps.empty() ? 1. : vMaps.front()->getPixelStep();
}

  
