// $Id: PolyMap.h,v 1.3 2010/07/03 23:03:26 garyb Exp $ 
// PixelMap that is affine transformation followed by 
// polynomial distortion, encoded into FITS headers
// as per SCAMP and the never-adopted standard for
// WCS polynomial distortions.
// Note that the "World" system is defined to be a gnomonic
// (tangent-plane) projection about the CRVAL[12] location,
// which is available via the getOrientation() method.

#ifndef POLYMAP_H
#define POLYMAP_H

#include "PixelMap.h"
#include "Poly2d.h"

namespace astrometry {

  class PolyMap: public PixelMap {
  public:
    PolyMap(const poly2d::Poly2d& px, const poly2d::Poly2d& py,
	    double tol_=0.001/3600.):
      xpoly(px), ypoly(py), worldTolerance(tol_) {}
    PolyMap(int orderx, int ordery, double tol_=0.001/3600.):
      xpoly(orderx,ordery), ypoly(orderx,ordery), worldTolerance(tol_) {}
    PolyMap(int order, double tol_=0.001/3600.):
      xpoly(order), ypoly(order), worldTolerance(tol_) {}
    // Note that default tolerance is set to be 1 mas if world units are degrees.

    ~PolyMap() {}

    // Implement all the virtual PixelMap calls that need to be overridden
    void toWorld(double xpix, double ypix,
		 double& xworld, double& yworld) const;
    // inverse map:
    void toPix( double xworld, double yworld,
		double &xpix, double &ypix) const;
    Matrix22 dWorlddPix(double xpix, double ypix) const;
    void toPix( double xworld, double yworld,
		double &xpix, double &ypix,
		DMatrix& derivs) const;
    void toWorld(double xpix, double ypix,
		 double& xworld, double& yworld,
		 DMatrix& derivs) const;

    void setParams(const DVector& p);
    DVector getParams() const;
    int nParams() const {return xpoly.nCoeffs() + ypoly.nCoeffs();}


    // Access routines for this derived class:
    poly2d::Poly2d getXPoly() const {return xpoly;}
    poly2d::Poly2d getYPoly() const {return ypoly;}

    // Set tolerance in world coords for soln of inverse
    void setWorldTolerance(double wt) {worldTolerance=wt;}

  private:
    poly2d::Poly2d xpoly;
    poly2d::Poly2d ypoly;
    double worldTolerance;
  };

  class LinearMap: public PixelMap {
  public:
    LinearMap(const DVector& v_): v(v_), vinv(6) {Assert(v.size()==6); makeInv();}
    LinearMap(): v(6,0.), vinv(6, 0.) {}
    ~LinearMap() {}

    // Implement all the virtual PixelMap calls that need to be overridden
    void toWorld(double xpix, double ypix,
		 double& xworld, double& yworld) const;
    // inverse map:
    void toPix( double xworld, double yworld,
		double &xpix, double &ypix) const;
    Matrix22 dWorlddPix(double xpix, double ypix) const;
    Matrix22 dPixdWorld(double xpix, double ypix) const;
    void toPix( double xworld, double yworld,
		double &xpix, double &ypix,
		DMatrix& derivs) const;
    void toWorld(double xpix, double ypix,
		 double& xworld, double& yworld,
		 DMatrix& derivs) const;

    void setParams(const DVector& p) {Assert(p.size()==6); v=p; makeInv();}
    DVector getParams() const {return v;}
    int nParams() const {return 6;}

  private:
    // Forward and inverse transformations
    DVector v;
    DVector vinv;
    void makeInv();
  };

} // namespace astrometry

#endif // POLYMAP_H
