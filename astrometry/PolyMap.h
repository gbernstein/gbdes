// Implementations of PixelMap that are Linear and higher-order 2d polynomial
// functions of coordinates.

#ifndef POLYMAP_H
#define POLYMAP_H

#include "PixelMap.h"
#include "Poly2d.h"

namespace astrometry {

  class PolyMap: public PixelMap {
  public:
    // Constructor with 1 order has terms with sum of x and y powers up to this order.
    // Constructor with 2 orders has all terms w/ powers of x up to orderx, y to ordery
    PolyMap(const poly2d::Poly2d& px, const poly2d::Poly2d& py, 
	    string name="",
	    double tol_=0.001/3600.):
    PixelMap(name), xpoly(px), ypoly(py), worldTolerance(tol_) {}
    PolyMap(int orderx, int ordery, string name="", double tol_=0.001/3600.):
      PixelMap(name), xpoly(orderx,ordery), ypoly(orderx,ordery), worldTolerance(tol_) {}
    PolyMap(int order, string name="", double tol_=0.001/3600.):
      PixelMap(name), xpoly(order), ypoly(order), worldTolerance(tol_) {}
    // Note that default tolerance is set to be 1 mas if world units are degrees.
    virtual PixelMap* duplicate() const {return new PolyMap(*this);}
      
    ~PolyMap() {}

    // Implement all the virtual PixelMap calls that need to be overridden
    void toWorld(double xpix, double ypix,
		 double& xworld, double& yworld) const;
    // inverse map:
    void toPix( double xworld, double yworld,
		double &xpix, double &ypix) const;
    Matrix22 dWorlddPix(double xpix, double ypix) const;
    void toPixDerivs( double xworld, double yworld,
		      double &xpix, double &ypix,
		      DMatrix& derivs) const;
    void toWorldDerivs(double xpix, double ypix,
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

    static string mapType() {return "Poly";}
    virtual string getType() const {return mapType();}
    static PixelMap* create(std::istream& is, string name="");
    void write(std::ostream& os) const;

  private:
    poly2d::Poly2d xpoly;
    poly2d::Poly2d ypoly;
    double worldTolerance;
  };

  class LinearMap: public PixelMap {
  public:
    LinearMap(const DVector& v_, string name=""): 
      PixelMap(name), v(v_), vinv(6) {Assert(v.size()==6); makeInv();}
    LinearMap(string name=""): PixelMap(name), v(6,0.), vinv(6, 0.) {}
    virtual PixelMap* duplicate() const {return new LinearMap(*this);}
    ~LinearMap() {}

    // Implement all the virtual PixelMap calls that need to be overridden
    void toWorld(double xpix, double ypix,
		 double& xworld, double& yworld) const;
    // inverse map:
    void toPix( double xworld, double yworld,
		double &xpix, double &ypix) const;
    Matrix22 dWorlddPix(double xpix, double ypix) const;
    Matrix22 dPixdWorld(double xpix, double ypix) const;
    void toPixDerivs( double xworld, double yworld,
		      double &xpix, double &ypix,
		      DMatrix& derivs) const;
    void toWorldDerivs(double xpix, double ypix,
		       double& xworld, double& yworld,
		       DMatrix& derivs) const;

    void setParams(const DVector& p) {Assert(p.size()==6); v=p; makeInv();}
    DVector getParams() const {return v;}
    int nParams() const {return 6;}

    static string mapType() {return "Linear";}
    virtual string getType() const {return mapType();}
    static PixelMap* create(std::istream& is, string name="");
    void write(std::ostream& os) const;

  private:
    // Forward and inverse transformations
    DVector v;
    DVector vinv;
    void makeInv();
  };

} // namespace astrometry

#endif // POLYMAP_H
