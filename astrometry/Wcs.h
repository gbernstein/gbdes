// Class that extends PixelMap by including a map to true celestial coordinates.
#ifndef WCS_H
#define WCS_H

#include "PixelMap.h"
#include "Astrometry.h"
#include "AstronomicalConstants.h"

namespace astrometry {
  class Wcs: public PixelMap {
  public:
    // Two constructors specify the PixelMap to be used, and the coordinate system
    // in which the PixelMap's world coords will be interpreted.  
    // The coordinate system is defined either by an Orientation to be used for a 
    // TangentPlane projection (gnomonic) or can give a SphericalCoords instance
    // that will be used.
    // wScale gives the units assumed for PixelMap world coords (SphericalCoords are always
    // in degrees)
    // ownPM is true if this class owns the PixelMap and should destroy it.
    // The nativeCoords will be duplicated.  If it's a TangentPlane,
    // note that it only contains a pointer to Orientation, which is not destroyed).
    Wcs(PixelMap* pm_, const Orientation& nativeOrient_, string name="", 
	double wScale_=DEGREE, bool ownPM=false);
    Wcs(PixelMap* pm_, const SphericalCoords& nativeCoords_, string name="", 
	double wScale_=DEGREE, bool ownPM=false);
    virtual ~Wcs();

    // Extend the PixelMap interface to include maps to/from true celestial coordinates
    SphericalICRS toSky(double xpix, double ypix) const;
    void fromSky(const SphericalCoords& sky, double& xpix, double& ypix) const;

    // Allow WCS to serve as PixelMap into a different projection from its native one:
    // If this is done via SphericalCoords, this class will assume ownership of it.
    void reprojectTo(const Orientation& targetOrient_);
    void reprojectTo(const SphericalCoords& targetCoords_);
    // Return to use of native projection as the world coords:
    void useNativeProjection();

    // Now we implement the PixelMap interface:
    static string mapType() {return "WCS";}
    static PixelMap* create(std::istream& is, string name="");
    virtual void write(std::ostream& os) const;

    virtual void toWorld(double xpix, double ypix,
			 double& xworld, double& yworld) const;
    virtual void toPix( double xworld, double yworld,
			double &xpix, double &ypix) const;
    virtual Matrix22 dPixdWorld(double xworld, double yworld) const;    
    virtual Matrix22 dWorlddPix(double xpix, double ypix) const;    
    virtual void toPixDerivs( double xworld, double yworld,
			      double &xpix, double &ypix,
			      DMatrix& derivs) const;
    virtual void toWorldDerivs(double xpix, double ypix,
			       double& xworld, double& yworld,
			       DMatrix& derivs) const;

    virtual void setParams(const DVector& p) {pm->setParams(p);}
    virtual DVector getParams() const {return pm->getParams();}
    virtual int nParams() const {return pm->nParams();}
    virtual double getPixelStep() const {return pm->getPixelStep();}
    virtual void setPixelStep(double p) {pm->setPixelStep(p);}

  private:
    PixelMap* pm;
    double wScale;
    bool ownPM;
    Orientation nativeOrient;
    mutable SphericalCoords* nativeCoords;
    Orientation targetOrient;
    mutable SphericalCoords* targetCoords;
  };
} // namespace astrometry
#endif // ifndef WCS_H
