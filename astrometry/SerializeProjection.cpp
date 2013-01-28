#include <sstream>
#include "SerializeProjection.h"
#include "StringStuff.h"

using namespace astrometry;

string
astrometry::serializeProjection(const SphericalCoords* projection) {
  const Gnomonic* gn = dynamic_cast<const Gnomonic*> (projection);
  if (gn) {
    std::ostringstream oss;
    oss << "Gnomonic " << *gn->getOrient();
    return oss.str();
  }
  const SphericalICRS* icrs = dynamic_cast<const SphericalICRS*> (projection);
  if (icrs) return "ICRS";
  const SphericalEcliptic* ecl = dynamic_cast<const SphericalEcliptic*> (projection);
  if (ecl) return "Ecliptic";

  throw AstrometryError("SerializeProjection() has unknown type of SphericalCoords");
}

SphericalCoords*
astrometry::deserializeProjection(string buffer) {
  string coordType;
  bool fail = false;
  istringstream iss(buffer);
  if ( !(iss >> coordType)) fail = true;
  if (!fail && stringstuff::nocaseEqual(coordType, "Gnomonic")) {
    Orientation orient;
    if (!(iss >> orient)) fail = true;
    return new Gnomonic(orient);
  } else if (!fail && stringstuff::nocaseEqual(coordType, "ICRS")) {
    return new SphericalICRS;
  } else if (!fail && stringstuff::nocaseEqual(coordType, "Ecliptic")) {
    return new SphericalEcliptic;
  }
  // Make it here if we have a read failure on the string
  throw AstrometryError("Error deserializing projection: <" + buffer +">");
}
