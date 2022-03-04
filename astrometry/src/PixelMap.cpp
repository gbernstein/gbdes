#include <cctype>
#include "PixelMap.h"
#include "StringStuff.h"

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
PixelMap::pixelArea(double xpix, double ypix, double color) const {
  double xw, yw, dxdx, dxdy, dydx, dydy;
  toWorld(xpix, ypix, xw, yw, color);
  toWorld(xpix+getPixelStep(), ypix, dxdx, dydx, color);
  dxdx -= xw; dydx-=yw;
  toWorld(xpix, ypix+getPixelStep(), dxdy, dydy, color);
  dxdy -= xw; dydy-=yw;
  return abs(dxdx*dydy-dxdy*dydx);
}

Matrix22
PixelMap::dPixdWorld(double xworld, double yworld, double color) const {
  double xpix=500., ypix=500.;
  toPix(xworld,yworld,xpix,ypix, color);
  return dWorlddPix(xpix,ypix,color).inverse();
}

Matrix22
PixelMap::dWorlddPix(double xpix, double ypix, double color) const {
  // Default implementation: use pixelStep for numerical derivatives
  const double step=getPixelStep();
  double xw0, yw0, xw1, yw1;
  toWorld(xpix, ypix, xw0, yw0, color);
  Matrix22 dd;
  toWorld(xpix+step, ypix, xw1, yw1, color);
  dd(0,0) = (xw1-xw0)/step;
  dd(1,0) = (yw1-yw0)/step;

  toWorld(xpix, ypix+step, xw1, yw1, color);
  dd(0,1) = (xw1-xw0)/step;
  dd(1,1) = (yw1-yw0)/step;
  return dd;
}

void
PixelMap::NewtonInverse(double xworld, double yworld,
			double& xpix, double& ypix,
			double worldTolerance,
			double color) const {
  // Iterate Newton's method to find pix coords given world
  // Start with input xpix/ypix as trial
  int nIter=0;
  const int maxIterations=10;
  while (nIter < maxIterations) {
    double xout, yout;
    toWorld(xpix,ypix,xout,yout,color);
    xout -= xworld;
    yout -= yworld;
    /**cerr << nIter << " world " << xworld << " " << yworld
	     << " pix " << xpix << " " << ypix 
	     << " resid " << xout << " " << yout
	     << endl; //**/
    if ( hypot(xout,yout)<worldTolerance) return;
    Matrix22 dPdW = dWorlddPix(xpix,ypix,color).inverse();
    xpix -= xout * dPdW(0,0) + yout * dPdW(0,1);
    ypix -= xout * dPdW(1,0) + yout * dPdW(1,1);
    nIter++;
  }
  throw AstrometryError("NewtonInverse did not converge");
}
    
void 
IdentityMap::toWorld(double xpix, double ypix,
		     double& xworld, double& yworld,
		     double color) const {
  xworld=xpix; yworld=ypix;
}
void 
IdentityMap::toPix( double xworld, double yworld,
		    double &xpix, double &ypix,
		    double color) const {
  xpix=xworld; ypix=yworld;
}

Matrix22
IdentityMap::dWorlddPix(double xpix, double ypix, double color) const {
  return Matrix22().setToIdentity();
}
Matrix22
IdentityMap::dPixdWorld(double xworld, double yworld, double color) const {
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
  return *this;
}

PixelMap*
ReprojectionMap::duplicate() const {
  return new ReprojectionMap(*this);
}
  
void 
ReprojectionMap::toWorld(double xpix, double ypix,
			 double& xworld, double& yworld,
			 double color) const {
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
			double &xpix, double &ypix,
			double color) const {
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
ReprojectionMap::dWorlddPix(double xpix, double ypix, double color) const {
  Matrix22 m;
  SphericalCoords* in = pix->duplicate();
  SphericalCoords* out = world->duplicate();
  in->setLonLat(xpix*scaleFactor, ypix*scaleFactor);
  out->convertFrom(*in,m);
  delete in; delete out;
  return m;
}
Matrix22 
ReprojectionMap::dPixdWorld(double xworld, double yworld, double color) const {
  Matrix22 m;
  SphericalCoords* in = world->duplicate();
  SphericalCoords* out = pix->duplicate();
  in->setLonLat(xworld*scaleFactor, yworld*scaleFactor);
  out->convertFrom(*in, m);
  delete in; delete out;
  return m;
}

#ifdef USE_YAML
PixelMap*
ReprojectionMap::create(const YAML::Node& node,
			bool& defaulted,
			string name) {
  string buffer;
  if (!node.IsMap() || 
      !node["Type"] || node["Type"].as<string>() != type() ||
      !node["InProjection"] ||
      !node["OutProjection"] ||
      !node["ScaleFactor"] )
    throw AstrometryError("ReprojectionMap::create() is missing YAML keys for name " + name);

  vector<SphericalCoords*> coords(2);  // Will read input and output projections
  coords[0] = SphericalCoords::deserialize(node["InProjection"]);
  coords[1] = SphericalCoords::deserialize(node["OutProjection"]);
  double scale = node["ScaleFactor"].as<double>();
  ReprojectionMap* out = new ReprojectionMap(*coords[0], *coords[1], scale, name);
  defaulted = false;
  delete coords[0];
  delete coords[1];
  return out;
}

void
ReprojectionMap::write(YAML::Emitter& os) const {
  os << YAML::BeginMap
     << YAML::Key << "Type" << YAML::Value << type()
     << YAML::Key << "InProjection"
     << YAML::Value << *pix
     << YAML::Key << "OutProjection"
     << YAML::Value << *world
     << YAML::Key << "ScaleFactor"
     << YAML::Value << scaleFactor
     << YAML::EndMap;
}
#endif

///////////////////////////////////////////////////////////////////
// Color term calculations
///////////////////////////////////////////////////////////////////
void
ColorTerm::toPix( double xworld, double yworld,
		  double &xpix, double &ypix,
		  double color) const {
  // A fast iterative solution assumes that the daughter map has xw-xp varying
  // very slowly as a function of xp across a range the size of the typical color
  // position shift.
  checkColor(color);
  xpix = xworld;
  ypix = yworld;
  double c = color - reference;
  for (int i=0; i<2; i++) {
    double xw1, yw1;
    pm->toWorld(xpix,ypix,xw1,yw1);
    xpix = xworld - c*(xw1-xpix);
    ypix = yworld - c*(yw1-ypix);
  }
}
Matrix22
ColorTerm::dPixdWorld(double xworld, double yworld, double color) const {
  double xpix, ypix;
  toPix(xworld,yworld,xpix,ypix,color);
  return dWorlddPix(xpix,ypix,color).inverse();
}

Matrix22
ColorTerm::dWorlddPix(double xpix, double ypix, 
		      double color) const {
  checkColor(color);
  double c = (color-reference); 
  Matrix22 m = c * pm->dWorlddPix(xpix,ypix);
  m(0,0) += 1.-c;
  m(1,1) += 1.-c;
  return m;
}

void
ColorTerm::toWorldDerivs(double xpix, double ypix,
			 double& xworld, double& yworld,
			 DMatrix& derivs,
			 double color) const {
  checkColor(color);
  pm->toWorldDerivs(xpix,ypix,xworld,yworld,derivs);
  double c = color - reference;
  xworld = xpix + c*(xworld-xpix);
  yworld = ypix + c*(yworld-ypix);
  derivs *= c;
}

void
ColorTerm::toPixDerivs( double xworld, double yworld,
			double &xpix, double &ypix,
			DMatrix& derivs,
			double color) const {
  // Do the quick iteration assuming small displacements and weak derivs
  checkColor(color);
  xpix = xworld;
  ypix = yworld;
  double c = color - reference;
  double xw1, yw1;
  pm->toWorld(xpix,ypix,xw1,yw1);
  xpix = xworld - c*(xw1-xpix);
  ypix = yworld - c*(yw1-ypix);
  pm->toWorldDerivs(xpix,ypix,xw1,yw1,derivs);
  xpix = xworld - c*(xw1-xpix);
  ypix = yworld - c*(yw1-ypix);
  derivs *= c;
}

#ifdef USE_YAML
void
ColorTerm::write(YAML::Emitter& os) const {
  os << YAML::BeginMap
     << YAML::Key << "Type" << YAML::Value << ColorTerm::type()
     << YAML::Key << "Reference" << YAML::Value << refColor()
     << YAML::Key << "Function" << YAML::Value;
  pm->write(os);
  os << YAML::EndMap;
}

PixelMap*
ConstantMap::create(const YAML::Node& node,
		    bool& defaulted,
		    string name) {
  if (!node.IsMap() || 
      !node["Type"] || node["Type"].as<string>() != type())
    throw AstrometryError("ConstantMap has missing <Type> fields at YAML node " + name);

  if (node["Parameters"]) {
    vector<double> vv = node["Parameters"].as<vector<double> >();
    if (vv.size()!=2)
      throw AstrometryError("ConstantMap has wrong # of parameters at YAML node " + name);
    defaulted = false;
    return new ConstantMap(name, vv[0], vv[1]);
  } else {
    // Build default
    defaulted = true;
    return new ConstantMap(name);
  }
}
#endif
