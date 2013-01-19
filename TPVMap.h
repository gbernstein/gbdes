// Read/write TPV (= TAN) format FITS WCS headers into our Wcs class.
// Map from pix to world is:
//   CRPIX-style linear map into degree units
//   PVx_y polynomial map into a (xi, eta), in degrees
//   Then gnomonic (de-)projection onto the sky.
// *** Documentation of the transforms at http://fits.gsfc.nasa.gov/registry/tpvwcs/tpv.html **

#ifndef TPVMAP_H
#define TPVMAP_H

#include "Wcs.h"
#include "PolyMap.h"
#include "Header.h"
#include "Bounds.h"

namespace astrometry {
  Wcs* readTPV(const img::Header& h);
  img::Header writeTPV(const Wcs& w);
  Wcs* fitTPV(Bounds<double> b,
	      const Wcs& wcsIn,
	      const SphericalCoords& tpvPole,
	      double tolerance=0.0001*ARCSEC/DEGREE);
} // namespace astrometry
#endif
