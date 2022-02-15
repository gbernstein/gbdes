// Definitions for 2d polynomial PixelMap

#include "PolyMap.h"
#include "StringStuff.h"
#include <sstream>

using namespace astrometry;

void
PolyMap::setDomain(Bounds<double> domain) {
  // Will be mapping the bounds into (-1,1) intervals
  if (domain.area()==0.)
    FormatAndThrow<AstrometryError>() << "Degenerate domain for PolyMap: " << domain;
  xshift = 0.5*(domain.getXMax() + domain.getXMin());
  xscale = 2. /(domain.getXMax() - domain.getXMin());
  yshift = 0.5*(domain.getYMax() + domain.getYMin());
  yscale = 2. /(domain.getYMax() - domain.getYMin());
}

Bounds<double>
PolyMap::getDomain() const {
  return Bounds<double>( xshift - 1./xscale,
			 xshift + 1./xscale,
			 yshift - 1./yscale,
			 yshift + 1./yscale);
}

void 
PolyMap::toWorld(double xpix, double ypix,
		 double& xworld, double& yworld,
		 double color) const {
  rescale(xpix, ypix);
  xworld = xpoly(xpix,ypix);
  yworld = ypoly(xpix,ypix);
}

void 
PolyMap::toPix( double xworld, double yworld,
		double &xpix, double &ypix,
		double color) const {
  // Note assuming that xpix and ypix are a decent guess on input
  NewtonInverse(xworld, yworld, xpix, ypix, worldTolerance, color);
}

Matrix22 
PolyMap::dWorlddPix(double xpix, double ypix, 
		    double color) const {
  Matrix22 d;
  rescale(xpix, ypix);
  d(0,0) = xpoly.derivx(xpix,ypix) * xscale;
  d(0,1) = xpoly.derivy(xpix,ypix) * yscale;
  d(1,0) = ypoly.derivx(xpix,ypix) * xscale;
  d(1,1) = ypoly.derivy(xpix,ypix) * yscale;
  return d;
}

void
PolyMap::setParams(const DVector& p) {
  Assert(p.size()==nParams());
  DVector vx=p.subVector(0,xpoly.nCoeffs());
  xpoly.setC(vx);
  DVector vy=p.subVector(xpoly.nCoeffs(),xpoly.nCoeffs()+ypoly.nCoeffs());
  ypoly.setC(vy);
}

DVector
PolyMap::getParams() const {
  DVector p(nParams(), 0.);
  p.subVector(0,xpoly.nCoeffs()) = xpoly.getC();
  p.subVector(xpoly.nCoeffs(),xpoly.nCoeffs()+ypoly.nCoeffs()) = ypoly.getC();
  return p;
}

void
PolyMap::setToIdentity() {
  if (xpoly.getOrderX()<1 || ypoly.getOrderY()<1)
    throw AstrometryError("PolyMap " + getName() + " has insufficient order to setToIdentity");
  DVector p = getParams();
  p.setZero();
  p[xpoly.vectorIndex(1,0)] = 1.;
  p[xpoly.nCoeffs() + ypoly.vectorIndex(0,1)] = 1.;
  setParams(p);
}

void 
PolyMap::toPixDerivs( double xworld, double yworld,
		      double &xpix, double &ypix,
		      DMatrix& derivs,
		      double color) const {
  toPix(xworld, yworld, xpix, ypix);
  double xx=xpix;
  double yy=ypix;
  rescale(xx,yy);
  DVector dx = xpoly.derivC(xx,yy);
  DVector dy = ypoly.derivC(xx,yy);
  Assert(derivs.cols()==2 && derivs.rows()==(dx.size()+dy.size()));
  derivs.setZero();
#ifdef USE_TMV
  derivs.row(0,0,dx.size()) = dx;
  derivs.row(1,dx.size(), dx.size()+dy.size()) = dy;
#elif defined USE_EIGEN
  derivs.block(0,0,1,dx.size()) = dx.transpose();
  derivs.block(1,dx.size(),1,dy.size()) = dy.transpose();
#endif
}

void 
PolyMap::toWorldDerivs(double xpix, double ypix,
		       double& xworld, double& yworld,
		       DMatrix& derivs,
		       double color) const {
  toWorld(xpix,ypix,xworld, yworld);
  rescale(xpix,ypix);
  DVector dx = xpoly.derivC(xpix, ypix);
  DVector dy = ypoly.derivC(xpix, ypix);
  Assert(derivs.rows()==2 && derivs.cols()==(dx.size()+dy.size()));
  derivs.setZero();
#ifdef USE_TMV
  derivs.row(0,0,dx.size()) = dx;
  derivs.row(1,dx.size(), dx.size()+dy.size()) = dy;
#elif defined USE_EIGEN
  derivs.block(0,0,1,dx.size()) = dx.transpose();
  derivs.block(1,dx.size(),1,dy.size()) = dy.transpose();
#endif
}

// A linear map too
void
LinearMap::toWorld(double xpix, double ypix,
		   double& xworld, double& yworld,
		   double color) const {
  xworld = v[0] + v[1]*xpix + v[2]*ypix;
  yworld = v[3] + v[4]*xpix + v[5]*ypix;
}
void 
LinearMap::toPix( double xworld, double yworld,
		  double &xpix, double &ypix,
		  double color) const {
  xpix = vinv[0] + vinv[1]*xworld + vinv[2]*yworld;
  ypix = vinv[3] + vinv[4]*xworld + vinv[5]*yworld;
}

Matrix22 
LinearMap::dWorlddPix(double xpix, double ypix, 
		      double color) const {
  Matrix22 m;
  m(0,0)=v[1];
  m(0,1)=v[2];
  m(1,0)=v[4];
  m(1,1)=v[5];
  return m;
  }
Matrix22 
LinearMap::dPixdWorld(double xworld, double yworld, 
		      double color) const {
  Matrix22 m;
  m(0,0)=vinv[1];
  m(0,1)=vinv[2];
  m(1,0)=vinv[4];
  m(1,1)=vinv[5];
  return m;
}

void 
LinearMap::toPixDerivs( double xworld, double yworld,
			double &xpix, double &ypix,
			DMatrix& derivs,
			double color) const {
  toPix(xworld, yworld, xpix, ypix);
  Assert(derivs.cols()==DIM && derivs.rows()==2);
  derivs.setZero();
  derivs(0,0) = 1.;
  derivs(0,1) = xpix;
  derivs(0,2) = ypix;
  derivs(1,3) = 1.;
  derivs(1,4) = xpix;
  derivs(1,5) = ypix;
}

void 
LinearMap::toWorldDerivs(double xpix, double ypix,
			 double& xworld, double& yworld,
			 DMatrix& derivs,
			 double color) const {
  toWorld(xpix, ypix, xworld, yworld);
  Assert(derivs.cols()==DIM && derivs.rows()==2);
  derivs.setZero();
  derivs(0,0) = 1.;
  derivs(0,1) = xpix;
  derivs(0,2) = ypix;
  derivs(1,3) = 1.;
  derivs(1,4) = xpix;
  derivs(1,5) = ypix;
}

void
LinearMap::makeInv() {
  double det = v[1]*v[5]-v[2]*v[4];
  vinv[1] = v[5]/det;
  vinv[2] = -v[2]/det;
  vinv[4] = -v[4]/det;
  vinv[5] = v[1]/det;
  vinv[0] = -(vinv[1]*v[0]+vinv[2]*v[3]);
  vinv[3] = -(vinv[4]*v[0]+vinv[5]*v[3]);
}

////////////////////////////////////////////////////////////////////////////////////
// YAML (de-)serializations
////////////////////////////////////////////////////////////////////////////////////
#ifdef USE_YAML

void
LinearMap::write(YAML::Emitter& os) const {
  vector<double> vv(DIM);
  for (int i=0; i<vv.size(); i++) vv[i] = v[i];
  os << YAML::BeginMap
     << YAML::Key << "Type" << YAML::Value << type()
     << YAML::Key << "Coefficients" << YAML::Flow << YAML::Value << vv
     << YAML::EndMap;
}

PixelMap*
LinearMap::create(const YAML::Node& node,
		  bool& defaulted,
		  string name) {
  if (!node.IsMap() || 
      !node["Type"] || node["Type"].as<string>() != type())
    throw AstrometryError("LinearMap has missing <Type> fields at YAML node " + name);

  if (node["Coefficients"]) {
    vector<string> vvStr = node["Coefficients"].as<vector<string> >();
    vector<double> vv(vvStr.size());
    transform(vvStr.begin(), vvStr.end(), vv.begin(), [](const string& val)
    {
    return std::stod(val);
    });
    if (vv.size()!=DIM)
      throw AstrometryError("LinearMap has wrong # of coefficients at YAML node " + name);
    DVector v(DIM);
    for (int i=0; i<v.size(); i++) v[i]=vv[i];
    defaulted = false;
    return new LinearMap(v, name);
  } else {
    // Build default
    defaulted = true;
    return new LinearMap(name);
  }
}

void
PolyMap::write(YAML::Emitter& os) const {
  Bounds<double> b = getDomain();
  os << YAML::BeginMap
     << YAML::Key << "Type" << YAML::Value << type()
     << YAML::Key << "XMin" << YAML::Value << b.getXMin()
     << YAML::Key << "XMax" << YAML::Value << b.getXMax()
     << YAML::Key << "YMin" << YAML::Value << b.getYMin()
     << YAML::Key << "YMax" << YAML::Value << b.getYMax()
     << YAML::Key << "XPoly" << YAML::Value;
  xpoly.write(os);
  os << YAML::Key << "YPoly" << YAML::Value;
  ypoly.write(os);
  os << YAML::Key << "Tolerance" << YAML::Value << worldTolerance
     << YAML::EndMap;
}
 
PixelMap*
PolyMap::create(const YAML::Node& node,
		bool& defaulted,
		string name) {
  if (!node.IsMap() || 
    !node["Type"] || node["Type"].as<string>() != type())
   throw AstrometryError("Missing PolyMap key <Type> at YAML node " + name);
  if (!node["XPoly"])
   throw AstrometryError("Missing PolyMap key <XPoly> at YAML node " + name);
  if (!node["YPoly"])
   throw AstrometryError("Missing PolyMap key <YPoly> at YAML node " + name);

  // Set up domain, starting with default
  Bounds<double> b(-1.,1.,-1.,1.);
  //if (node["XMin"]) b.setXMax(node["XMin"].as<double>());
  if (node["XMin"]) {
    // TODO: fix whatever is going on with inputYAML so that this string-to-double isn't required:
    b.setXMin(stod(node["XMin"].as<string>()));
  }
  
  //if (node["XMax"]) b.setXMax(node["XMax"].as<double>());
  if (node["XMax"]) b.setXMax(stod(node["XMax"].as<string>()));
  //if (node["YMin"]) b.setYMin(node["YMin"].as<double>());
  if (node["YMin"]) b.setYMin(stod(node["YMin"].as<string>()));
  //if (node["YMax"]) b.setYMax(node["YMax"].as<double>());
  if (node["YMax"]) b.setYMax(stod(node["YMax"].as<string>()));
  
  bool defaultX, defaultY;
  poly2d::Poly2d* px = poly2d::Poly2d::create(node["XPoly"], defaultX);
  poly2d::Poly2d* py = poly2d::Poly2d::create(node["YPoly"], defaultY);
  defaulted = defaultX || defaultY;
  double tol = 0.;  // zero value will indicate default
  if (node["Tolerance"])
    tol = stod(node["Tolerance"].as<string>());

  PolyMap* pm = 0;
  if (tol>0.)
    pm = new PolyMap(*px, *py, name, b, tol);
  else {
    // Take default on tolerance
    pm = new PolyMap(*px, *py, name, b);
  }

  // If the linear terms are zero, we will assume that we have a map
  // that was initialized with all zeros.  This means we should in fact
  // be initializing to default identity polynomials
  int orderx = px->getOrderX();
  int ordery = px->getOrderY();
  bool bad = false;
  if (orderx>0 && ordery>0) {
    // One of the two linear terms needs to be nonzero to be non-default
    bad =  px->getC()[px->vectorIndex(1,0)]==0. && px->getC()[px->vectorIndex(0,1)]==0.;
  } else if (orderx > 0) {
    // The x linear term needs to be nonzero
    bad =  px->getC()[px->vectorIndex(1,0)]==0.;
  } else if (ordery > 0) {
    // The y linear term needs to be nonzero
    bad =  px->getC()[px->vectorIndex(0,1)]==0.;
  } else {
    // The polynomial is a constant, which is no good...
    throw AstrometryError("Cannot use constant polynomial for PolyMap from YAML node "+name);
  }
  // Now check the y polynomial too
  orderx = py->getOrderX();
  ordery = py->getOrderY();
  if (orderx>0 && ordery>0) {
    // One of the two linear terms needs to be nonzero to be non-default
    bad = bad || ( py->getC()[py->vectorIndex(1,0)]==0. && 
		   py->getC()[py->vectorIndex(0,1)]==0.);
  } else if (orderx > 0) {
    // The x linear term needs to be nonzero
    bad = bad || py->getC()[py->vectorIndex(1,0)]==0.;
  } else if (ordery > 0) {
    // The y linear term needs to be nonzero
    bad = bad || py->getC()[py->vectorIndex(0,1)]==0.;
  } else {
    // The polynomial is a constant, which is no good...
    throw AstrometryError("Cannot use constant polynomial for PolyMap from YAML node "+name);
  }

  delete px;
  delete py;

  if (bad)
    pm->setToIdentity();
  
  return pm;
}

#endif // end YAML requirement
