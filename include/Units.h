// Declarations of standard units for use in catalog inputs, internal calcs, and outputs
#ifndef UNITS_H
#define UNITS_H

#include "AstronomicalConstants.h"

namespace astrometry {
  // Unit to be used input celestial coordinates and units of angular
  // distance for other world coord systems (e.g. in gnomonic projection coords).
  // Also expected for all input and internal position errors.
  // Note that Astrometry routines will be working in radians though for RA,Dec.
  // Probably a bad idea to change this one since it has to agree with FITS
  // standards that are making input WCS's.
  const double WCS_UNIT = DEGREE;
  //
  // Unit used when inputting/outputing astrometric residuals or RMS or 2d covariance
  const double RESIDUAL_UNIT = MILLIARCSEC;
  //
  // Units of parallax an proper motion (input and output)
  // Standard is that PM is always oriented in local ICRS (E,N) coordinates.
  // Internally we will use WCS_UNIT for parallax and PM (per TDB_UNIT).
  const double PARALLAX_UNIT = MILLIARCSEC;
  const double PM_UNIT = MILLIARCSEC/YEAR;
  //
  // Units of time (TDB's)
  const double TDB_UNIT = YEAR; 
  // Standard order of 5d params
  enum PM {X0, Y0, VX, VY, PAR};
  //
  // Observatory positions assumed to be in ICRS AU.

} // end namespace astrometry
#endif
