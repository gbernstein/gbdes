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
unique_ptr<Wcs> readTPV(const img::Header &h, string name = "");
shared_ptr<Wcs> readTPVFromSIP(const map<string, double> &header, string name = "");
img::Header writeTPV(const Wcs &w);  // Will throw exception if Wcs is wrong form
// Fit a TPV model to wcsIn over the range specified by b.  The TPV system will
// have its CRVAL[01] and the tpvPole.
// The transformation will use the given color for any color terms.
// Order of the fitting polynomial can be specified (stick to <=5).  Entry
// of -1 means it will start at 3 and continue until RMS is <tolerance or 5 is
// done.
unique_ptr<Wcs> fitTPV(Bounds<double> b, const Wcs &wcsIn, const SphericalCoords &tpvPole, string name = "",
            double color = 0., double tolerance = 0.0001 * ARCSEC / DEGREE, double order = -1);
// Get N-choose-K for binomial expansion:
double NChooseK(int n, int k);
// Use binomial expansion to get exponential of binomial:
DMatrix binom(DVector xyTerm, int power);
}  // namespace astrometry
#endif
