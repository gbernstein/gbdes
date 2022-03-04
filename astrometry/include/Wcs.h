// Class that extends PixelMap by including a map to true celestial coordinates.
#ifndef WCS_H
#define WCS_H

#include "PixelMap.h"
#include "Astrometry.h"
#include "AstronomicalConstants.h"

namespace astrometry {
  class Wcs: public PixelMap {
  public:
    // Constructor specifies the PixelMap to be used, and the coordinate system
    // in which the PixelMap's world coords will be interpreted.  
    // The coordinate system is defined by a SphericalCoords instance that will
    // be duplicated and saved with the Wcs.
    // Pixel coords are mapped to "world" by the PixelMap, then scaled by wScale
    // (defaults to DEGREE) to give the (lon,lat) or (xi,eta) coords
    // of the SphericalCoords object.  SphericalCoords now gives location on the sky.
    // For shareMap=false, a duplicate of the input map is owned by this class and deleted
    // in destructor.  For shareMap=true, uses the input pointer but does not delete it.
    Wcs(PixelMap* pm_, const SphericalCoords& nativeCoords_, string name="", 
	double wScale_=DEGREE, bool shareMap_=false);
    virtual ~Wcs();
    virtual Wcs* duplicate() const;  // Note that duplicate has same sharing behavior as this.

    // Extend the PixelMap interface to include maps to/from true celestial coordinates
    SphericalICRS toSky(double xpix, double ypix,
			double color=astrometry::NODATA) const;
    void fromSky(const SphericalCoords& sky, 
		 double& xpix, double& ypix,
		 double color=astrometry::NODATA) const;

    // Allow WCS to serve as PixelMap into a different projection from its native one:
    // A duplicate of the targetCoords is created and stored with the Wcs.
    void reprojectTo(const SphericalCoords& targetCoords_);
    // Return to use of native projection as the world coords:
    void useNativeProjection();

    // Now we implement the PixelMap interface:
    static string type() {return "WCS";}
    virtual string getType() const {return type();}

    // Serialization is here for convenience, since PixelMapCollection will
    // handle serialization most of the time.  This does NOT serialize
    // the wrapped PixelMap.  No create() here; this will be handled via
    // PixelMapCollection in current code.
    
#ifdef USE_YAML
    virtual void write(YAML::Emitter& os) const; 
#endif
    
    virtual void toWorld(double xpix, double ypix,
			 double& xworld, double& yworld,
			 double color=astrometry::NODATA) const;
    virtual void toPix( double xworld, double yworld,
			double &xpix, double &ypix,
			double color=astrometry::NODATA) const;
    virtual Matrix22 dPixdWorld(double xworld, double yworld,
				double color=astrometry::NODATA) const;    
    virtual Matrix22 dWorlddPix(double xpix, double ypix,
				double color=astrometry::NODATA) const;    
    virtual void toPixDerivs( double xworld, double yworld,
			      double &xpix, double &ypix,
			      DMatrix& derivs,
			      double color=astrometry::NODATA) const;
    virtual void toWorldDerivs(double xpix, double ypix,
			       double& xworld, double& yworld,
			       DMatrix& derivs,
			       double color=astrometry::NODATA) const;

    virtual void setParams(const DVector& p) {pm->setParams(p);}
    virtual DVector getParams() const {return pm->getParams();}
    virtual int nParams() const {return pm->nParams();}
    virtual double getPixelStep() const {return pm->getPixelStep();}
    virtual void setPixelStep(double p) {pm->setPixelStep(p);}
    virtual bool needsColor() const {return pm->needsColor();}

    const PixelMap* getMap() const {return pm;}
    const SphericalCoords* getNativeCoords() const {return nativeCoords;}
    const SphericalCoords* getTargetCoords() const {return targetCoords;}
    double getScale() const {return wScale;}

  private:
    PixelMap* pm;
    double wScale;
    bool shareMap;	// True if the pm is owned by someone else
    SphericalCoords* nativeCoords;
    SphericalCoords* targetCoords;

    // Hide these to avoid inadvertent ownership confusion:
    Wcs(const Wcs& rhs);
    void operator=(const Wcs& rhs);
  };
} // namespace astrometry
#endif // ifndef WCS_H
