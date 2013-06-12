// Definitions for 2d polynomial PixelMap

#include "PolyMap.h"
#include "StringStuff.h"
#include <sstream>

using namespace astrometry;

void 
PolyMap::toWorld(double xpix, double ypix,
		 double& xworld, double& yworld) const {
  xworld = xpoly(xpix,ypix);
  yworld = ypoly(xpix,ypix);
}

void 
PolyMap::toPix( double xworld, double yworld,
		double &xpix, double &ypix) const {
  // Note assuming that xpix and ypix are a decent guess on input
  NewtonInverse(xworld, yworld, xpix, ypix, worldTolerance);
}

Matrix22 
PolyMap::dWorlddPix(double xpix, double ypix) const {
  Matrix22 d;
  d(0,0) = xpoly.derivx(xpix,ypix);
  d(0,1) = xpoly.derivy(xpix,ypix);
  d(1,0) = ypoly.derivx(xpix,ypix);
  d(1,1) = ypoly.derivy(xpix,ypix);
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
		      DMatrix& derivs) const {
  toPix(xworld, yworld, xpix, ypix);
  DVector dx = xpoly.derivC(xpix, ypix);
  DVector dy = ypoly.derivC(xpix, ypix);
  Assert(derivs.ncols()==2 && derivs.nrows()==(dx.size()+dy.size()));
  derivs.setZero();
  derivs.row(0,0,dx.size()) = dx;
  derivs.row(1,dx.size(), dx.size()+dy.size()) = dy;
}

void 
PolyMap::toWorldDerivs(double xpix, double ypix,
		       double& xworld, double& yworld,
		       DMatrix& derivs) const {
  toWorld(xpix,ypix,xworld, yworld);
  DVector dx = xpoly.derivC(xpix, ypix);
  DVector dy = ypoly.derivC(xpix, ypix);
  Assert(derivs.nrows()==2 && derivs.ncols()==(dx.size()+dy.size()));
  derivs.setZero();
  derivs.row(0,0,dx.size()) = dx;
  derivs.row(1,dx.size(), dx.size()+dy.size()) = dy;
}

// A linear map too
void
LinearMap::toWorld(double xpix, double ypix,
		   double& xworld, double& yworld) const {
  xworld = v[0] + v[1]*xpix + v[2]*ypix;
  yworld = v[3] + v[4]*xpix + v[5]*ypix;
}
void 
LinearMap::toPix( double xworld, double yworld,
		  double &xpix, double &ypix) const {
  xpix = vinv[0] + vinv[1]*xworld + vinv[2]*yworld;
  ypix = vinv[3] + vinv[4]*xworld + vinv[5]*yworld;
}
Matrix22 
LinearMap::dWorlddPix(double xpix, double ypix) const {
  Matrix22 m;
  m(0,0)=v[1];
  m(0,1)=v[2];
  m(1,0)=v[4];
  m(1,1)=v[5];
  return m;
}
Matrix22 
LinearMap::dPixdWorld(double xworld, double yworld) const {
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
			DMatrix& derivs) const {
  toPix(xworld, yworld, xpix, ypix);
  Assert(derivs.ncols()==6 && derivs.nrows()==2);
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
			 DMatrix& derivs) const {
  toWorld(xpix, ypix, xworld, yworld);
  Assert(derivs.ncols()==6 && derivs.nrows()==2);
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

PixelMap*
LinearMap::create(std::istream& is, string name) {
  DVector v(6);
  string buffer;
  getlineNoComment(is, buffer);
  istringstream iss(buffer);
  if (!(iss >> v[0] >> v[1] >> v[2]))
    throw AstrometryError("LinearMap::create has bad x coefficients: " + buffer);
  getlineNoComment(is, buffer);
  iss.str(buffer);
  iss.clear();
  if (!(iss >> v[3] >> v[4] >> v[5])) 
    throw AstrometryError("LinearMap::create has bad y coefficients: " + buffer);
  return new LinearMap(v);
}

void
LinearMap::write(std::ostream& os, int precision) const {
  stringstuff::StreamSaver ss(os);
  os.precision(precision);
  os.setf( ios_base::showpos | ios_base::scientific);
  os << v[0] << " " << v[1] << " " << v[2] << endl;
  os << v[3] << " " << v[4] << " " << v[5] << endl;
}

void
PolyMap::write(std::ostream& os, int precision) const {
  xpoly.write(os, precision);
  ypoly.write(os, precision);
  // ??? Worry about precision on the tolerance line
  os << worldTolerance << endl;
}

PixelMap*
PolyMap::create(std::istream& is, string name) {
  poly2d::Poly2d* px = poly2d::Poly2d::create(is);
  poly2d::Poly2d* py = poly2d::Poly2d::create(is);
  double tol;
  string buffer;
  if (!getlineNoComment(is, buffer)) 
    throw AstrometryError("PolyMap::create() is missing tolerance value at name " + name);
  istringstream iss(buffer);
  if (!(iss >> tol))
    throw AstrometryError("PolyMap::create() has bad tolerance value: " + buffer);
  PolyMap* pm =  new PolyMap(*px, *py, name, tol);
  delete px;
  delete py;
  return pm;
}

