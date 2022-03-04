// Implementations of PixelMap that are Linear and higher-order 2d polynomial
// functions of coordinates.

#ifndef POLYMAP_H
#define POLYMAP_H

#include "PixelMap.h"
#include "Poly2d.h"
#include "Bounds.h"
#ifdef USE_YAML
#include "yaml-cpp/yaml.h"
#endif

namespace astrometry {

  class PolyMap: public PixelMap {
  public:
    // Constructor with 1 order has terms with sum of x and y powers up to this order.
    // Constructor with 2 orders has all terms w/ powers of x up to orderx, y to ordery
    // Maps set to identity on construction unless polynomial coefficients are given
    PolyMap(const poly2d::Poly2d& px, const poly2d::Poly2d& py, 
	    string name="",
	    Bounds<double> domain=Bounds<double>(-1.,1.,-1.,1.),
	    double tol_=0.001/3600.):
      PixelMap(name), xpoly(px), ypoly(py), worldTolerance(tol_) {setDomain(domain);}
      
    PolyMap(int orderx, int ordery,
	    string name="",
	    Bounds<double> domain=Bounds<double>(-1.,1.,-1.,1.),
	    double tol_=0.001/3600.):
      PixelMap(name), xpoly(orderx,ordery), ypoly(orderx,ordery), worldTolerance(tol_) {
	setDomain(domain);
	setToIdentity();
      }
      
    PolyMap(int order,
	    string name="",
	    Bounds<double> domain=Bounds<double>(-1.,1.,-1.,1.),
	    double tol_=0.001/3600.):
      PixelMap(name), xpoly(order), ypoly(order), worldTolerance(tol_) {
	setDomain(domain);
	setToIdentity();
      }
    // Note that default tolerance is set to be 1 mas if world units are degrees.
    virtual PixelMap* duplicate() const {return new PolyMap(*this);}
      
    ~PolyMap() {}

    // Implement all the virtual PixelMap calls that need to be overridden
    void toWorld(double xpix, double ypix,
		 double& xworld, double& yworld,
		 double color=astrometry::NODATA) const;
    // inverse map:
    void toPix( double xworld, double yworld,
		double &xpix, double &ypix,
		double color=astrometry::NODATA) const;
    Matrix22 dWorlddPix(double xpix, double ypix, 
			double color=astrometry::NODATA) const;
    void toPixDerivs( double xworld, double yworld,
		      double &xpix, double &ypix,
		      DMatrix& derivs,
		      double color=astrometry::NODATA) const;
    void toWorldDerivs(double xpix, double ypix,
		       double& xworld, double& yworld,
		       DMatrix& derivs,
		       double color=astrometry::NODATA) const;

    void setParams(const DVector& p);
    DVector getParams() const;
    int nParams() const {return xpoly.nCoeffs() + ypoly.nCoeffs();}


    // Access routines for this derived class:
    poly2d::Poly2d getXPoly() const {return xpoly;}
    poly2d::Poly2d getYPoly() const {return ypoly;}

    // Set tolerance in world coords for soln of inverse
    void setWorldTolerance(double wt) {worldTolerance=wt;}
    // Set coefficients to give identity transformation:
    void setToIdentity();
    // Set up domain scaling
    void setDomain(Bounds<double> domain);
    // Return rectangular bounds that rescale to (-1,1)
    Bounds<double> getDomain() const;
		   
    static string type() {return "Poly";}
    virtual string getType() const {return type();}
#ifdef USE_YAML
    static PixelMap* create(const YAML::Node& node,
			    bool& defaulted,
			    string name="");
    void write(YAML::Emitter& os) const;
#endif
    
  private:
    poly2d::Poly2d xpoly;
    poly2d::Poly2d ypoly;
    // scale input data so that x->xscale*(x-shift), nominally (-1,1.)
    double xshift, xscale;
    double yshift, yscale;
    void rescale(double& xpix, double& ypix) const {
      // Apply domain scaling
      xpix = xscale * (xpix-xshift);
      ypix = yscale * (ypix-yshift);
    }
    double worldTolerance;
  };

  // ??? rescale linearMap xpix,ypix too?
  
  class LinearMap: public PixelMap {
  public:
    LinearMap(const DVector& v_, string name=""): 
      PixelMap(name), v(v_), vinv(DIM) {Assert(v.size()==DIM); makeInv();}
    // The default constructor sets the transformation to identity
    LinearMap(string name=""): PixelMap(name), v(DIM,0.), vinv(DIM, 0.) {setToIdentity();}
    virtual PixelMap* duplicate() const {return new LinearMap(*this);}
    ~LinearMap() {}

    // Implement all the virtual PixelMap calls that need to be overridden
    void toWorld(double xpix, double ypix,
		 double& xworld, double& yworld,
		 double color=astrometry::NODATA) const;
    // inverse map:
    void toPix( double xworld, double yworld,
		double &xpix, double &ypix,
		double color=astrometry::NODATA) const;
    Matrix22 dWorlddPix(double xpix, double ypix,
			double color=astrometry::NODATA) const;
    Matrix22 dPixdWorld(double xworld, double yworld, 
			double color=astrometry::NODATA) const;
    void toPixDerivs( double xworld, double yworld,
		      double &xpix, double &ypix,
		      DMatrix& derivs,
		      double color=astrometry::NODATA) const;
    void toWorldDerivs(double xpix, double ypix,
		       double& xworld, double& yworld,
		       DMatrix& derivs,
		       double color=astrometry::NODATA) const;

    void setParams(const DVector& p) {Assert(p.size()==DIM); v=p; makeInv();}
    DVector getParams() const {return v;}
    int nParams() const {return DIM;}

    void setToIdentity() {v.setZero(); v[1]=1.; v[5]=1.; vinv.setZero(); vinv[1]=1.; vinv[5]=1.;}

    static string type() {return "Linear";}
    virtual string getType() const {return type();}
#ifdef USE_YAML
    static PixelMap* create(const YAML::Node& node,
			    bool& defaulted,
			    string name="");
    void write(YAML::Emitter& os) const;
#endif
  private:
    // Forward and inverse transformations are always kept in sync
    const static int DIM=6;
    DVector v;
    DVector vinv;
    void makeInv();
  };

} // namespace astrometry

#endif // POLYMAP_H
