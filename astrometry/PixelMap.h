// $Id: PixelMap.h,v 1.10 2012/02/02 23:34:25 garyb Exp $
// Match several images together and derive position (mag?) map tweaks
// to obtain optimal agreement.
#ifndef PIXELMAP_H
#define PIXELMAP_H

#include <string>
#include <list>
#include "Std.h"
#include "UseTMV.h"
#include "Astrometry.h"

// Note that Astrometry.h has typdefs for the 2- and 3-element small vectors & matrices

namespace astrometry {

  class PixelMap {
  public:
    PixelMap(): pixStep(1.) {}
    virtual ~PixelMap() {}
    // pixel coords to world coords map
    virtual void toWorld(double xpix, double ypix,
			 double& xworld, double& yworld) const=0;
    // inverse map:
    virtual void toPix( double xworld, double yworld,
			double &xpix, double &ypix) const=0;
    virtual Matrix22 dPixdWorld(double xworld, double yworld) const;
    virtual Matrix22 dWorlddPix(double xpix, double ypix) const;
    // Solid angle of a pixel (sr)
    double pixelArea(double xpix, double ypix) const;

    // Get parameters of the mapping.  Note default behavior of derived classes
    // is to ignore anything to do with map parameters.
    virtual void setParams(const DVector& p) {}
    // Note here that a zero-length vector is returned as parameters by default
    virtual DVector getParams() const {return DVector(0);}
    virtual int nParams() const {return 0;}
    // Here the derivatives are of world coords w.r.t. map parameters
    virtual void toPix( double xworld, double yworld,
			double &xpix, double &ypix,
			DMatrix& derivs) const {
      toPix(xworld,yworld,xpix,ypix);}
    // And derivs of world w.r.t. parameters of map.  Resizing of derivs
    // done if needed.
    virtual void toWorld(double xpix, double ypix,
			 double& xworld, double& yworld,
			 DMatrix& derivs) const {
      toWorld(xpix,ypix,xworld,yworld);}
    // Step size for derivatives in "pixel" space:
    virtual double getPixelStep() const {return pixStep;}
    virtual void setPixelStep(double ps) {pixStep=ps;}
  private:
    double pixStep;
  protected:
    // Invert the pix->world map using Newton iteration.
    // Likely that many derived PixelMap classes will want this.
    void NewtonInverse(double xworld, double yworld, 
		       double& xpix, double& pix,
		       double worldTolerance) const;
  };

  class CompoundPixelMap: public PixelMap {
  private:
    // PixelMap whose pix - to - world map follows series of transformations in order of list
    // Keeps pointers to elements of compound map; destruction is up to user.
    list<PixelMap*> pmlist;
  protected:
    CompoundPixelMap() {}
  public:
    CompoundPixelMap(PixelMap* pm1) {pmlist.push_back(pm1);}
    ~CompoundPixelMap() {}
    void append(PixelMap* pmnew) {pmlist.push_back(pmnew);}
    void prepend(PixelMap* pmnew) {pmlist.push_front(pmnew);}
    // pixel coords to world coords map
    virtual void toWorld(double xpix, double ypix,
			 double& xworld, double& yworld) const;
    // inverse map:
    virtual void toPix( double xworld, double yworld,
			double &xpix, double &ypix) const;

    // Get parameters of the mapping.  Note default behavior of derived classes
    // is to ignore anything to do with map parameters.
    virtual void setParams(const DVector& p);
    virtual DVector getParams() const;
    virtual int nParams() const;

    // Here the derivatives are of world coords w.r.t. map parameters
    virtual void toPix( double xworld, double yworld,
			double &xpix, double &ypix,
			DMatrix& derivs) const;
    // And derivs of world w.r.t. parameters of map.  Resizing of derivs
    // done if needed.
    virtual void toWorld(double xpix, double ypix,
			 double& xworld, double& yworld,
			 DMatrix& derivs) const;
    virtual double getPixelStep() const {return pmlist.front()->getPixelStep();}
    virtual void setPixelStep(double p) {return pmlist.front()->setPixelStep(p);}
    // ??? should set all steps by propagating derivatives??
  };

  class IdentityMap: public PixelMap {
  public:
    void toWorld(double xpix, double ypix,
		 double& xworld, double& yworld) const;
    void toPix( double xworld, double yworld,
		double &xpix, double &ypix) const;
    Matrix22 dPixdWorld(double xworld, double yworld) const;    
    Matrix22 dWorlddPix(double xpix, double ypix) const;    
    // Call base class for these other signatures of toPix and toWorld
    // (otherwise these calls are hidden for the derived class, weird C++ thing)
    virtual void toPix( double xworld, double yworld,
			double &xpix, double &ypix,
			DMatrix& derivs) const {PixelMap::toPix(xworld,yworld,xpix,ypix,derivs);}
    virtual void toWorld(double xpix, double ypix,
			 double& xworld, double& yworld,
			 DMatrix& derivs) const {PixelMap::toPix(xpix,ypix,xworld,yworld,derivs);}
  };

  // PixelMap that takes coordinates in one projection
  // into another.  No adjustable parameters
  class ReprojectionMap: public PixelMap {
  public:
    // On creation, give reference to a set of coordinates in the
    // source (pix) and destination (world) systems.  
    // scaleFactor will be applied to the lat/lon in both systems
    // to put them into the native radian units of this class.
    // For instance enter AstronomicalConstants::DEGREE to
    // have coordinates in degrees.
    // setPixelScale to 1 arcsecond for any numerical derivatives of this projection,
    // although derivatives are done analytically.
    ReprojectionMap(const SphericalCoords& pixCoords,
		    const SphericalCoords& worldCoords,
		    double scale_=1.):
    pix(pixCoords.duplicate()), world(worldCoords.duplicate()), scaleFactor(scale_) {
      setPixelStep(ARCSEC/scale_);
    }
    ~ReprojectionMap() {delete pix; delete world;}
    void toWorld(double xpix, double ypix,
		 double& xworld, double& yworld) const;
    void toPix( double xworld, double yworld,
		double &xpix, double &ypix) const;
    Matrix22 dPixdWorld(double xworld, double yworld) const;    
    Matrix22 dWorlddPix(double xpix, double ypix) const;    
    // Call base class for these other signatures of toPix and toWorld
    // (otherwise these calls are hidden for the derived class, weird C++ thing)
    virtual void toPix( double xworld, double yworld,
			double &xpix, double &ypix,
			DMatrix& derivs) const {PixelMap::toPix(xworld,yworld,xpix,ypix,derivs);}
    virtual void toWorld(double xpix, double ypix,
			 double& xworld, double& yworld,
			 DMatrix& derivs) const {PixelMap::toPix(xpix,ypix,xworld,yworld,derivs);}

  private:
    SphericalCoords* pix;
    SphericalCoords* world;
    double scaleFactor;
  };

} // namespace astrometry
#endif //PIXMAP_H
