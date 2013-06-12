#include <cctype>
#include "PixelMap.h"
#include "StringStuff.h"
#include "SerializeProjection.h"

using namespace astrometry;

int 
PixelMap::anonymousCounter=0;

PixelMap::PixelMap(string name_): name(name_), pixStep(1.) {
  if (name.empty()) {
    std::ostringstream oss;
    oss << "map_" << anonymousCounter++;
    name = oss.str();
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
  // Default implementation: use pixelStep for numerical derivatives
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

ReprojectionMap::ReprojectionMap(const ReprojectionMap& rhs): PixelMap(rhs),
							      pix(rhs.pix->duplicate()), 
							      world(rhs.world->duplicate()), 
							      scaleFactor(rhs.scaleFactor)  {
  setPixelStep(rhs.getPixelStep());
}

ReprojectionMap&
ReprojectionMap::operator=(const ReprojectionMap& rhs) {
  if (this == &rhs) return *this;  // self-assignment
  delete pix;
  delete world;
  PixelMap::operator=(rhs);  // Sets name & pixel step
  pix = rhs.pix->duplicate(); 
  world = rhs.world->duplicate(); 
  scaleFactor = rhs.scaleFactor;
}

PixelMap*
ReprojectionMap::duplicate() const {
  return new ReprojectionMap(*this);
}
  
void 
ReprojectionMap::toWorld(double xpix, double ypix,
			 double& xworld, double& yworld) const {
  // Creation of duplicates of SphericalCoords allows this to be const and thread-safe.
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
  // Creation of duplicates of SphericalCoords allows this to be const and thread-safe.
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



PixelMap*
ReprojectionMap::create(std::istream& is, string name) {
  string buffer;
  vector<SphericalCoords*> coords(2);  // Will read input and output projections
  for (int i=0; i<coords.size(); i++) {
    if (!stringstuff::getlineNoComment(is, buffer)) 
      throw AstrometryError("Stream input failure in ReprojectionMap::create() for name " + name);
    coords[i] = deserializeProjection(buffer);
  }

  // Now read the scale factor:
  if (!stringstuff::getlineNoComment(is, buffer)) 
    throw AstrometryError("Stream input failure in ReprojectionMap::create() for name " + name);
  istringstream iss(buffer);
  double scale;
  if (!(iss >> scale))
    throw AstrometryError("Bad scale input to ReprojectionMap::create(): " + buffer);

  ReprojectionMap* out = new ReprojectionMap(*coords[0], *coords[1], scale, name);
  delete coords[0];
  delete coords[1];
  return out;
}

void
ReprojectionMap::write(std::ostream& os) const {
  os << serializeProjection(pix) << endl;
  os << serializeProjection(world) << endl;
  os << scaleFactor << endl;
}

