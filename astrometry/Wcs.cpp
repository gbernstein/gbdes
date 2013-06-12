// Implementation of Wcs class
#include "Wcs.h"

using namespace astrometry;

Wcs::Wcs(PixelMap* pm_, const SphericalCoords& nativeCoords_, string name, 
	 double wScale_, bool shareMap_): PixelMap(name), wScale(wScale_), 
					 nativeCoords(0),
					 targetCoords(0),
					 shareMap(shareMap_)
{
  nativeCoords = nativeCoords_.duplicate();
  if (shareMap) {
    pm = pm_;
  } else {
    pm = pm_->duplicate();
  }
}

Wcs::~Wcs() {
  if (nativeCoords) delete nativeCoords;
  if (targetCoords) delete targetCoords;
  if (!shareMap) delete pm;
}

Wcs*
Wcs::duplicate() const {
  Wcs* retval = new Wcs(pm, *nativeCoords, getName(), wScale, shareMap);
  retval->targetCoords = targetCoords->duplicate();
}

SphericalICRS 
Wcs::toSky(double xpix, double ypix) const {
  double xw,yw;
  pm->toWorld(xpix, ypix, xw, yw);
  // Create duplicate coordinates for thread safety:
  SphericalCoords* nc = nativeCoords->duplicate();
  nc->setLonLat(xw*wScale, yw*wScale);
  SphericalICRS retval(*nativeCoords);
  delete nc;
  return retval;
}

void 
Wcs::fromSky(const SphericalCoords& sky, double& xpix, double& ypix) const {
  // Create duplicate coordinates for thread safety:
  SphericalCoords* nc = nativeCoords->duplicate();
  nc->convertFrom(sky);
  double xw, yw;
  nc->getLonLat(xw, yw);
  xw /= wScale;
  yw /= wScale;
  pm->toPix(xw, yw, xpix, ypix);
  delete nc;
}

void
Wcs::reprojectTo(const SphericalCoords& targetCoords_) {
  useNativeProjection(); // clear out old reprojection info
  // Will not bother checking to see if the new projection matches old one.
  targetCoords = targetCoords_.duplicate();
}

void
Wcs::useNativeProjection() {
  if (targetCoords) {
    delete targetCoords;
    targetCoords = 0;
  }
}  

// Now we implement the PixelMap interface.  In each case, we check targetCoords to see
// if we need to reproject the underlying PixelMap's "world" coordinates from the native 
// to the target projection before passing them out as the "world" coordinates of this class.

void 
Wcs::toWorld(double xpix, double ypix,
	     double& xworld, double& yworld) const {
  pm->toWorld(xpix, ypix, xworld, yworld);
  if (targetCoords) {
    // Create duplicates for thread safety:
    SphericalCoords* nc = nativeCoords->duplicate();
    SphericalCoords* tc = targetCoords->duplicate();
    nc->setLonLat(xworld * wScale, yworld * wScale);
    tc->convertFrom(*nc);
    tc->getLonLat(xworld, yworld);
    xworld /= wScale;
    yworld /= wScale;
    delete tc;
    delete nc;
  }
}

void 
Wcs::toPix( double xworld, double yworld,
	    double &xpix, double &ypix) const {
  if (targetCoords) {
    // Create duplicates for thread safety:
    SphericalCoords* nc = nativeCoords->duplicate();
    SphericalCoords* tc = targetCoords->duplicate();
    tc->setLonLat(xworld*wScale, yworld*wScale);
    nc->convertFrom(*nc);
    nc->getLonLat(xworld, yworld);
    xworld /= wScale;
    yworld /= wScale;
    delete tc;
    delete nc;
  }
  pm->toPix(xworld, yworld, xpix, ypix);
}

Matrix22 
Wcs::dWorlddPix(double xpix, double ypix) const {
  Matrix22 m1 = pm->dWorlddPix(xpix, ypix);
  if (!targetCoords) return m1;

  // Create duplicates for thread safety:
  SphericalCoords* nc = nativeCoords->duplicate();
  SphericalCoords* tc = targetCoords->duplicate();
  Matrix22 mReproject;
  double xw, yw;
  pm->toWorld(xpix, ypix, xw, yw);
  nc->setLonLat(xw * wScale, yw * wScale);
  tc->convertFrom(*nc, mReproject);
  delete tc;
  delete nc;
  return mReproject*m1;
}

Matrix22 
Wcs::dPixdWorld(double xworld, double yworld) const {
  if (!targetCoords) return pm->dPixdWorld(xworld, yworld);

  // Create duplicates for thread safety:
  SphericalCoords* nc = nativeCoords->duplicate();
  SphericalCoords* tc = targetCoords->duplicate();
  Matrix22 mReproject;
  tc->setLonLat(xworld*wScale, yworld*wScale);
  nc->convertFrom(*tc, mReproject);
  nc->getLonLat(xworld, yworld);
  delete tc;
  delete nc;
  return pm->dPixdWorld(xworld/wScale, yworld/wScale) * mReproject;
}

void 
Wcs::toWorldDerivs(double xpix, double ypix,
		   double& xworld, double& yworld,
		   DMatrix& derivs) const {
  pm->toWorldDerivs(xpix, ypix, xworld, yworld, derivs);
  if (targetCoords) {
    // Create duplicates for thread safety:
    SphericalCoords* nc = nativeCoords->duplicate();
    SphericalCoords* tc = targetCoords->duplicate();
    Matrix22 mReproject;
    nc->setLonLat(xworld * wScale, yworld * wScale);
    tc->convertFrom(*nc, mReproject);
    tc->getLonLat(xworld, yworld);
    xworld /= wScale;
    yworld /= wScale;
    if (pm->nParams()>0) {
      // Avoid possible aliasing:
      DMatrix tmp = mReproject * derivs;
      derivs = tmp;  
    }
    delete tc;
    delete nc;
  }
}

void 
Wcs::toPixDerivs( double xworld, double yworld,
		  double &xpix, double &ypix,
		  DMatrix& derivs) const {
  if (!targetCoords) pm->toPixDerivs(xworld, yworld, xpix, ypix, derivs);

  // Create duplicates for thread safety:
  SphericalCoords* nc = nativeCoords->duplicate();
  SphericalCoords* tc = targetCoords->duplicate();
  tc->setLonLat(xworld*wScale, yworld*wScale);
  nc->convertFrom(*tc);
  nc->getLonLat(xworld, yworld);
  pm->toPixDerivs(xworld, yworld, xpix, ypix, derivs);
  if (pm->nParams()>0) {
    Matrix22 mReproject;
    // We need the derivative dW/dP for the reprojection here:
    tc->convertFrom(*nc, mReproject);
    // Avoid possible aliasing:
    DMatrix tmp = mReproject * derivs;
    derivs = tmp;  
  }
  delete tc;
  delete nc;
}
