// Implementation of Wcs class
#include "Wcs.h"

using namespace astrometry;

Wcs::Wcs(PixelMap* pm_, const SphericalCoords& nativeCoords_, string name, 
	 double wScale_, bool shareMap): PixelMap(name), wScale(wScale), 
					 nativeCoords(0),
					 targetCoords(0),
					 ownMap(!shareMap)
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
  if (ownMap) delete pm;
}

Wcs*
Wcs::duplicate() const {
  Wcs* retval = new Wcs(pm, *nativeCoords, getName(), wScale, !ownMap);
  retval->targetCoords = targetCoords->duplicate();
}

SphericalICRS 
Wcs::toSky(double xpix, double ypix) const {
  double xw,yw;
  pm->toWorld(xpix, ypix, xw, yw);
  nativeCoords->setLonLat(xw*wScale, yw*wScale);
  return SphericalICRS( *nativeCoords);
}

void 
Wcs::fromSky(const SphericalCoords& sky, double& xpix, double& ypix) const {
  nativeCoords->convertFrom(sky);
  double xw, yw;
  nativeCoords->getLonLat(xw, yw);
  xw /= wScale;
  yw /= wScale;
  pm->toPix(xw, yw, xpix, ypix);
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
    nativeCoords->setLonLat(xworld * wScale, yworld * wScale);
    targetCoords->convertFrom(*nativeCoords);
    targetCoords->getLonLat(xworld, yworld);
    xworld /= wScale;
    yworld /= wScale;
  }
}

void 
Wcs::toPix( double xworld, double yworld,
	    double &xpix, double &ypix) const {
  if (targetCoords) {
    targetCoords->setLonLat(xworld*wScale, yworld*wScale);
    nativeCoords->convertFrom(*targetCoords);
    nativeCoords->getLonLat(xworld, yworld);
    xworld /= wScale;
    yworld /= wScale;
  }
  pm->toPix(xworld, yworld, xpix, ypix);
}

Matrix22 
Wcs::dWorlddPix(double xpix, double ypix) const {
  Matrix22 m1 = pm->dWorlddPix(xpix, ypix);
  if (!targetCoords) return m1;
  Matrix22 mReproject;
  double xw, yw;
  pm->toWorld(xpix, ypix, xw, yw);
  nativeCoords->setLonLat(xw * wScale, yw * wScale);
  targetCoords->convertFrom(*nativeCoords, mReproject);
  return mReproject*m1;
}

Matrix22 
Wcs::dPixdWorld(double xworld, double yworld) const {
  if (!targetCoords) return pm->dPixdWorld(xworld, yworld);
  Matrix22 mReproject;
  targetCoords->setLonLat(xworld*wScale, yworld*wScale);
  nativeCoords->convertFrom(*targetCoords, mReproject);
  nativeCoords->getLonLat(xworld, yworld);
  return pm->dPixdWorld(xworld/wScale, yworld/wScale) * mReproject;
}

void 
Wcs::toWorldDerivs(double xpix, double ypix,
		   double& xworld, double& yworld,
		   DMatrix& derivs) const {
  pm->toWorldDerivs(xpix, ypix, xworld, yworld, derivs);
  if (targetCoords) {
    Matrix22 mReproject;
    nativeCoords->setLonLat(xworld * wScale, yworld * wScale);
    targetCoords->convertFrom(*nativeCoords, mReproject);
    targetCoords->getLonLat(xworld, yworld);
    xworld /= wScale;
    yworld /= wScale;
    if (pm->nParams()>0) {
      // Avoid possible aliasing:
      DMatrix tmp = mReproject * derivs;
      derivs = tmp;  
    }
  }
}

void 
Wcs::toPixDerivs( double xworld, double yworld,
		  double &xpix, double &ypix,
		  DMatrix& derivs) const {
  if (!targetCoords) pm->toPixDerivs(xworld, yworld, xpix, ypix, derivs);
  targetCoords->setLonLat(xworld*wScale, yworld*wScale);
  nativeCoords->convertFrom(*targetCoords);
  nativeCoords->getLonLat(xworld, yworld);
  pm->toPixDerivs(xworld, yworld, xpix, ypix, derivs);
  if (pm->nParams()>0) {
    Matrix22 mReproject;
    // We need the derivative dW/dP for the reprojection here:
    targetCoords->convertFrom(*nativeCoords, mReproject);
    // Avoid possible aliasing:
    DMatrix tmp = mReproject * derivs;
    derivs = tmp;  
  }
}

/// ???? Read and write
