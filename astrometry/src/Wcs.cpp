// Implementation of Wcs class
#include "Wcs.h"

using namespace astrometry;

Wcs::Wcs(PixelMap* pm_, const SphericalCoords& nativeCoords_, string name, 
	 double wScale_, bool shareMap_): PixelMap(name), wScale(wScale_), 
					 nativeCoords(nullptr),
					 targetCoords(nullptr),
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
  return retval;
}

SphericalICRS 
Wcs::toSky(double xpix, double ypix, double color) const {
  double xw,yw;
  pm->toWorld(xpix, ypix, xw, yw,color);
  // Create duplicate coordinates for thread safety:
  SphericalCoords* nc = nativeCoords->duplicate();
  nc->setLonLat(xw*wScale, yw*wScale);
  SphericalICRS retval(*nc);
  delete nc;
  return retval;
}

void 
Wcs::fromSky(const SphericalCoords& sky, 
	     double& xpix, double& ypix,
	     double color) const {
  // Create duplicate coordinates for thread safety:
  SphericalCoords* nc = nativeCoords->duplicate();
  nc->convertFrom(sky);
  double xw, yw;
  nc->getLonLat(xw, yw);
  xw /= wScale;
  yw /= wScale;
  pm->toPix(xw, yw, xpix, ypix, color);
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
    targetCoords = nullptr;
  }
}  

// Now we implement the PixelMap interface.  In each case, we check targetCoords to see
// if we need to reproject the underlying PixelMap's "world" coordinates from the native 
// to the target projection before passing them out as the "world" coordinates of this class.

void 
Wcs::toWorld(double xpix, double ypix,
	     double& xworld, double& yworld,
	     double color) const {
  pm->toWorld(xpix, ypix, xworld, yworld, color);
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
	    double &xpix, double &ypix,
	    double color) const {
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
  pm->toPix(xworld, yworld, xpix, ypix, color);
}

Matrix22 
Wcs::dWorlddPix(double xpix, double ypix, 
		double color) const {
  Matrix22 m1 = pm->dWorlddPix(xpix, ypix, color);
  if (!targetCoords) return m1;

  // Create duplicates for thread safety:
  SphericalCoords* nc = nativeCoords->duplicate();
  SphericalCoords* tc = targetCoords->duplicate();
  Matrix22 mReproject;
  double xw, yw;
  pm->toWorld(xpix, ypix, xw, yw, color);
  nc->setLonLat(xw * wScale, yw * wScale);
  tc->convertFrom(*nc, mReproject);
  delete tc;
  delete nc;
  return mReproject*m1;
}

Matrix22 
Wcs::dPixdWorld(double xworld, double yworld,
		double color) const {
  if (!targetCoords) return pm->dPixdWorld(xworld, yworld, color);

  // Create duplicates for thread safety:
  SphericalCoords* nc = nativeCoords->duplicate();
  SphericalCoords* tc = targetCoords->duplicate();
  Matrix22 mReproject;
  tc->setLonLat(xworld*wScale, yworld*wScale);
  nc->convertFrom(*tc, mReproject);
  nc->getLonLat(xworld, yworld);
  delete tc;
  delete nc;
  return pm->dPixdWorld(xworld/wScale, yworld/wScale, color) * mReproject;
}

void 
Wcs::toWorldDerivs(double xpix, double ypix,
		   double& xworld, double& yworld,
		   DMatrix& derivs,
		   double color) const {
  pm->toWorldDerivs(xpix, ypix, xworld, yworld, derivs, color);
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
		  DMatrix& derivs,
		  double color) const {
  if (!targetCoords) pm->toPixDerivs(xworld, yworld, xpix, ypix, derivs, color);

  // Create duplicates for thread safety:
  SphericalCoords* nc = nativeCoords->duplicate();
  SphericalCoords* tc = targetCoords->duplicate();
  tc->setLonLat(xworld*wScale, yworld*wScale);
  nc->convertFrom(*tc);
  nc->getLonLat(xworld, yworld);
  pm->toPixDerivs(xworld, yworld, xpix, ypix, derivs, color);
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

#ifdef USE_YAML
void
Wcs::write(YAML::Emitter& os) const {
  os << YAML::BeginMap
     << YAML::Key << "Type" << YAML::Value << type()
     << YAML::Key << "MapName" << YAML::Value << pm->getName()
     << YAML::Key << "Projection" << YAML::Value << *nativeCoords
     << YAML::Key << "Scale" << YAML::Value << wScale;
  if (targetCoords)
    os << YAML::Key << "Target_Projection" << YAML::Value << *targetCoords;
  os << YAML::EndMap;
  return;
}
#endif
