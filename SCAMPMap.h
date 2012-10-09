// $Id: SCAMPMap.h,v 1.4 2011/10/29 22:10:23 garyb Exp $
// Read SCAMP-style coordinate maps from header files.
// Map from pix to world is:
//   CRPIX-style linear map into degree units
//   PVx_y polynomial map into a gnomonic system, in degrees
//   optionally project into new system.
// *** Documentation of the transforms at http://fits.gsfc.nasa.gov/registry/tpvwcs/tpv.html **

#ifndef SCAMPMAP_H
#define SCAMPMAP_H

#include "Astrometry.h"
#include "PolyMap.h"
#include "Header.h"
#include "Bounds.h"

namespace astrometry {
  class SCAMPMap: public CompoundPixelMap {
  public:
    // Map is into tangent system at SCAMP's orientation,
    // unless a pointer is given to a preferred output projection.
    // Note that this class copies the Orientation given
    //   for the reprojection, as well as owning the component PixelMaps.
    // Degree units are adopted for world systems here.
    SCAMPMap(const img::Header& h,
	     const Orientation* reproject=0);
    ~SCAMPMap();
    SphericalICRS toICRS(double xpix, double ypix) const {
      double lon, lat;
      toWorld(xpix, ypix, lon, lat);
      tp->setLonLat(lon*DEGREE, lat*DEGREE);
      return SphericalICRS(*tp);
    }
    const LinearMap* linear() const {return lm;}
    // Note that there may be no PolyMap (pm=0) if
    // SCAMPMap was read from a header that did not give PV coeffs at all
    // and hence takes identity in place of polynomial step.
    const PolyMap* poly() const {return pm;}
    const Orientation& orientFITS() const {return orientIn;}
    const Orientation& projection() const {return orientOut;}
  private:
    LinearMap* lm;
    PolyMap* pm;
    ReprojectionMap* rm;
    TangentPlane* tp;
    Orientation orientIn;
    Orientation orientOut;
  };

  img::Header FitSCAMP(Bounds<double> b,
		       const PixelMap& pm,
		       const Orientation& pmOrient,
		       const SphericalCoords& pole,
		       double tolerance=0.0001*ARCSEC/DEGREE);
} // namespace astrometry
#endif
