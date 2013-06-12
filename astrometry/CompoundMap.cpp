// Code for DEPRECATED CompoundMap class

#include <cctype>
#include "CompoundMap.h"
#include "StringStuff.h"

using namespace astrometry;

// Deep copy
PixelMap*
CompoundPixelMap::duplicate() const {
  // Note possibility of endless recursion from circular dependence here!
  CompoundPixelMap* retval = new CompoundPixelMap(pmlist.front()->duplicate(), getName());
  list<PixelMap*>::const_iterator i = pmlist.begin();
  for ( ++i; i != pmlist.end(); ++i)
    retval->append((*i)->duplicate());
  return retval;
}

void
CompoundPixelMap::toWorld(double xpix, double ypix,
			  double& xworld, double& yworld) const {
  double x=xpix;
  double y=ypix;
  for (list<PixelMap*>::const_iterator i=pmlist.begin();
       i != pmlist.end();
       ++i) {
    (*i)->toWorld(x,y,x,y);
  }
  xworld = x;
  yworld = y;
}

void
CompoundPixelMap::toPix(double xworld, double yworld,
			double& xpix, double& ypix) const {
  double x=xworld;
  double y=yworld;
  // xpix and ypix on input are supposed to be guesses at the
  // answer.  Propagate them forward to create a list of guesses
  // of "pixel" coords for all the intermediate steps.
  list<double> xguess;
  list<double> yguess;
  for (list<PixelMap*>::const_iterator i=pmlist.begin();
       i != pmlist.end(); i++) {
    xguess.push_back(xpix);
    yguess.push_back(ypix);
    (*i)->toWorld(xpix,ypix,xpix,ypix);
  }
  list<double>::const_reverse_iterator ix=xguess.rbegin();
  list<double>::const_reverse_iterator iy=yguess.rbegin();

  for (list<PixelMap*>::const_reverse_iterator i=pmlist.rbegin();
      i != pmlist.rend();
       ++i, ++ix, ++iy) {
    double xg = *ix;
    double yg = *iy;
    (*i)->toPix(x,y,xg,yg);
    x = xg;
    y = yg;
  }
  xpix = x;
  ypix = y;
}

int
CompoundPixelMap::nParams() const {
  int n=0;
  for (list<PixelMap*>::const_iterator i=pmlist.begin();
       i != pmlist.end();
       ++i) n+= (*i)->nParams();
  return n;
}

DVector
CompoundPixelMap::getParams() const {
  int n=0;
  DVector p(nParams());
  int index = 0;
  for (list<PixelMap*>::const_iterator i=pmlist.begin();
       i != pmlist.end();
       ++i) {
    DVector pp = (*i)->getParams();
    for (int j=0; j<(*i)->nParams(); j++, index++) p[index]=pp[j];
  }
  return p;
}

void
CompoundPixelMap::setParams(const DVector& p) {
  Assert(p.size() == nParams());
  int n=0;
  int index = 0;
  for (list<PixelMap*>::const_iterator i=pmlist.begin();
       i != pmlist.end();
       ++i) {
    DVector pp((*i)->nParams());
    for (int j=0; j<pp.size(); j++, index++) pp[j]=p[index];
    (*i)->setParams(pp);
  }
}

void 
CompoundPixelMap::toWorldDerivs( double xpix, double ypix,
				 double &xworld, double &yworld,
				 DMatrix& derivs) const {
  double x=xpix;
  double y=ypix;
  Assert(derivs.nrows()==2 && derivs.ncols()==nParams());
  derivs.setZero();

  int index=0;	// starting index of parameters for this map
  for (list<PixelMap*>::const_iterator i=pmlist.begin();
       i != pmlist.end();
       ++i) {
    if (index>0) {
      DMatrix tmp = (*i)->dWorlddPix(x,y) * derivs.colRange(0,index);
      derivs.colRange(0,index) = tmp; 
    }
    int nNext = (*i)->nParams();
    if (nNext>0) {
      DMatrix dd(2,nNext);
      (*i)->toWorldDerivs(x,y,x,y,dd);
      derivs.colRange(index, index+nNext) += dd;
      index += nNext;
    } else {
      (*i)->toWorld(x,y,x,y);
    }
  }
  xworld = x; yworld = y;
}

void 
CompoundPixelMap::toPixDerivs( double xworld, double yworld,
			       double &xpix, double &ypix,
			       DMatrix& derivs) const {
  toPix(xworld,yworld,xpix,ypix);
  double x=xpix;
  double y=ypix;
  
  Assert(derivs.nrows()==2 && derivs.ncols()==nParams());
  derivs.setZero();

  int index=0;	// starting index of parameters for this map
  for (list<PixelMap*>::const_iterator i=pmlist.begin();
       i != pmlist.end();
       ++i) {
    if (index>0) {
      //      derivs.colRange(0,index) = (*i)->dWorlddPix(x,y) * derivs.colRange(0,index);
      DMatrix tmp = (*i)->dWorlddPix(x,y) * derivs.colRange(0,index);
      derivs.colRange(0,index) = tmp; 
    }
    int nNext = (*i)->nParams();
    if (nNext>0) {
      DMatrix dd(2,nNext);
      (*i)->toWorldDerivs(x,y,x,y,dd);
      derivs.colRange(index, index+nNext) += dd;
      index += nNext;
    } else {
      (*i)->toWorld(x,y,x,y);
    }
  }
}
