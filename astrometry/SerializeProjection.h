// Take a SphericalCoord instance and turn into a string that tell us what the projection
// from (lon,lat) coordinates to sky is.  And then rebuild an appropriate SphericalCoords from 
// such a string:
#ifndef SERIALIZE_PROJECTION_H
#define SERIALIZE_PROJECTION_H
#include "Astrometry.h"

namespace astrometry {
  string serializeProjection(const SphericalCoords* projection);
  SphericalCoords*  deserializeProjection(string s);
} // end namespace astrometry

#endif  // SERIALIZE_PROJECTION_H
