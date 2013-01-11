// $Id: PixelMap.cpp,v 1.12 2012/02/02 02:23:31 garyb Exp $

#include "PixelMap.h"
#include <cctype>

using namespace astrometry;

int 
PixelMap::anonymousCounter=0;

PixelMap::PixelMap(string name_): name(name_), pixStep(1.) {
  if (name.empty()) {
    std::ostringstream oss;
    oss << "map_" << anonymousCounter++;
    name = oss.str();
  } else {
    // Replace white space with underscores in PixelMap names:
    for (string::iterator i=name.begin(); i!=name.end(); ++i)
      if (std::isspace(*i)) *i='_';
  }
}

double
PixelMap::pixelArea(double xpix, double ypix) const {
  double xw, yw, dxdx, dxdy, dydx, dydy;
  toWorld(xpix, ypix, xw, yw);
  toWorld(xpix+getPixelStep(), ypix, dxdx, dydx);
  dxdx -= xw; dydx-=yw;
  toWorld(xpix, ypix+getPixelStep(), dxdy, dydy);
  dxdy -= xw; dydy-=yw;
  return abs(dxdx*dydy-dxdy*dydx);
}

Matrix22
PixelMap::dPixdWorld(double xworld, double yworld) const {
  double xpix, ypix;
  toPix(xworld,yworld,xpix,ypix);
  return dWorlddPix(xpix,ypix).inverse();
}

Matrix22
PixelMap::dWorlddPix(double xpix, double ypix) const {
  // Use 1 pixel as a step for finite derivatives
  const double step=getPixelStep();
  double xw0, yw0, xw1, yw1;
  toWorld(xpix, ypix, xw0, yw0);
  Matrix22 dd;
  toWorld(xpix+step, ypix, xw1, yw1);
  dd(0,0) = (xw1-xw0)/step;
  dd(1,0) = (yw1-yw0)/step;

  toWorld(xpix, ypix+step, xw1, yw1);
  dd(0,1) = (xw1-xw0)/step;
  dd(1,1) = (yw1-yw0)/step;
  return dd;
}

void
PixelMap::NewtonInverse(double xworld, double yworld,
			double& xpix, double& ypix,
			double worldTolerance) const {
  // Iterate Newton's method to find pix coords given world
  // Start with input xpix/ypix as trial
  int nIter=0;
  const int maxIterations=10;
  while (nIter < maxIterations) {
    double xout, yout;
    toWorld(xpix,ypix,xout,yout);
    xout -= xworld;
    yout -= yworld;
    /**cerr << nIter << " world " << xworld << " " << yworld
	     << " pix " << xpix << " " << ypix 
	     << " resid " << xout << " " << yout
	     << endl; //**/
    if ( hypot(xout,yout)<worldTolerance) return;
    Matrix22 dPdW = dWorlddPix(xpix,ypix).inverse();
    xpix -= xout * dPdW(0,0) + yout * dPdW(0,1);
    ypix -= xout * dPdW(1,0) + yout * dPdW(1,1);
    nIter++;
  }
  throw AstrometryError("NewtonInverse did not converge");
}
    
void 
IdentityMap::toWorld(double xpix, double ypix,
		     double& xworld, double& yworld) const {
  xworld=xpix; yworld=ypix;
}
void 
IdentityMap::toPix( double xworld, double yworld,
		    double &xpix, double &ypix) const {
  xpix=xworld; ypix=yworld;
}

Matrix22
IdentityMap::dWorlddPix(double xpix, double ypix) const {
  return Matrix22().setToIdentity();
}
Matrix22
IdentityMap::dPixdWorld(double xworld, double yworld) const {
  return Matrix22().setToIdentity();
}

void 
ReprojectionMap::toWorld(double xpix, double ypix,
			 double& xworld, double& yworld) const {
  SphericalCoords* in = pix->duplicate();
  in->setLonLat(xpix*scaleFactor, ypix*scaleFactor);
  SphericalCoords* out = world->duplicate();
  out->convertFrom(*in);
  Vector2 vv = out->getLonLat() / scaleFactor;
  xworld=vv[0];
  yworld=vv[1];
  delete in; delete out;
}
void 
ReprojectionMap::toPix( double xworld, double yworld,
			double &xpix, double &ypix) const {
  SphericalCoords* in = world->duplicate();
  SphericalCoords* out = pix->duplicate();
  in->setLonLat(xworld*scaleFactor, yworld*scaleFactor);
  out->convertFrom(*in);
  Vector2 vv = out->getLonLat() / scaleFactor;
  xpix=vv[0];
  ypix=vv[1];
  delete in; delete out;
}
Matrix22 
ReprojectionMap::dWorlddPix(double xpix, double ypix) const {
  Matrix22 m;
  SphericalCoords* in = pix->duplicate();
  SphericalCoords* out = world->duplicate();
  in->setLonLat(xpix*scaleFactor, ypix*scaleFactor);
  out->convertFrom(*in,m);
  delete in; delete out;
  return m;
}
Matrix22 
ReprojectionMap::dPixdWorld(double xworld, double yworld) const {
  Matrix22 m;
  SphericalCoords* in = world->duplicate();
  SphericalCoords* out = pix->duplicate();
  in->setLonLat(xworld*scaleFactor, yworld*scaleFactor);
  out->convertFrom(*in, m);
  delete in; delete out;
  return m;
}


// Chain together some maps:
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
