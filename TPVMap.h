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
#include "PixelMapCollection.h"

namespace astrometry {
  // Convert a Wcs expressible as FITS "TPV" WCS standard to/from a list of header key/values
  // The Wcs will have a Gnomonic sky system, and a PixelMap that is linear followed by
  // optional polynomial.
  Wcs* readTPV(const img::Header& h, string name="");
  img::Header writeTPV(const Wcs& w);   // Will throw exception if Wcs is wrong form
  Wcs* fitTPV(Bounds<double> b,
	      const Wcs& wcsIn,
	      const Orientation& tpvOrient,
	      string name="",
	      double tolerance=0.0001*ARCSEC/DEGREE);
} // namespace astrometry
#endif
