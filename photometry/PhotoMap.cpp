#include <cctype>
#include "PhotoMap.h"
#include "StringStuff.h"

using namespace photometry;

// Counter used to assign names to PhotoMaps created without specified name
int 
PhotoMap::anonymousCounter=0;

PhotoMap::PhotoMap(string name_): name(name_) {
  if (name.empty()) {
    std::ostringstream oss;
    oss << "map_" << anonymousCounter++;
    name = oss.str();
  }
}

// Default implement of inverse will assume we have an additive forward
// transform that is independent of the input mag.
double
PhotoMap::inverse(double magOut, const PhotoArguments& args) const {
  double dMag = forward(magOut, args) - magOut;
  return magOut - dMag;
}

void
PhotoMap::NewtonInverse(double magOut,
			double& magIn,
			const PhotoArguments& args,
			double magTolerance) const {
  // Iterate Newton's method to invert magnitude map
  // Use the input magIn as the initial guess
  int nIter=0;
  const int maxIterations=10;
  while (nIter < maxIterations) {
    double outTrial = forward(magIn, args);
    outTrial -= magOut;
    if ( abs(outTrial)<magTolerance) return;
    magIn -= outTrial / derivative(magIn,args);
    nIter++;
  }
  throw PhotometryError("NewtonInverse did not converge");
}

double
PolyMap::forward(double magIn, const PhotoArguments& args) const {
  double dMag;
  if (useExposureCoords)
    dMag += poly(args.xExposure, args.yExposure);
  else
    dMag += poly(args.xDevice, args.yDevice);
  return magIn + dMag;
}

double
PolyMap::forwardDerivs(double magIn, 
		       const PhotoArguments& args,
		       DVector& derivs) const {
  if (useExposureCoords)
    derivs = poly.derivC(args.xExposure,args.yExposure);
  else
    derivs = poly.derivC(args.xDevice, args.yDevice);
  return forward(magIn, args);
}

void
PolyMap::setParams(const DVector& p) {
  Assert(p.size()==nParams());
  poly.setC(p);
}

PolyMap::PolyMap(const poly2d::Poly2d& p,
		 ArgumentType argType,
		 string name): PhotoMap(name), poly(p) {
  setArgumentType(argType);
}

PolyMap::PolyMap(int orderx, int ordery,
		 ArgumentType argType,
		 string name): PhotoMap(name), poly(orderx,ordery) {
  setArgumentType(argType);
}

PolyMap::PolyMap(int order, 
		 ArgumentType argType,
		 string name): PhotoMap(name), poly(order) {
  setArgumentType(argType);
}

void
PolyMap::setArgumentType(ArgumentType argType) {
  switch (argType) {
  case Device:
    useExposureCoords = false;
    break;
  case Exposure:
    useExposureCoords = true;
    break;
  default:
    throw PhotometryError("Unknown ArgumentType");
  }
}

PhotoMap*
PolyMap::create(std::istream& is, string name) {
  string buffer;
  if (!getlineNoComment(is, buffer)) 
    throw PhotometryError("PolyMap::create() is missing inputs for name " + name);

  istringstream iss(buffer);
  // use device coordinates or exposure coordinates?
  string coordType;
  iss >> coordType;
  ArgumentType argType;
  if (stringstuff::nocaseEqual(coordType,"Device"))
    argType = Device;
  else if (stringstuff::nocaseEqual(coordType,"Exposure"))
    argType = Exposure;
  else
    throw PhotometryError("PolyMap::create() does not know argument type <" + coordType + ">");

  poly2d::Poly2d* p = poly2d::Poly2d::create(is);
  PolyMap* pm =  new PolyMap(*p, argType, name);
  delete p;
  return pm;
}

void
PolyMap::write(std::ostream& os, int precision) const {
  os << (useExposureCoords ? "Exposure" : "Device") << endl;
  poly.write(os, precision);
}

PhotoMap*
ConstantMap::create(std::istream& is, string name) {
  // Serialization of ConstantMap is just a single number (the constant)
  string buffer;
  if (!getlineNoComment(is, buffer)) 
    throw PhotometryError("ConstantMap::create() is missing inputs for name " + name);

  istringstream iss(buffer);
  double c;
  if (!(iss >> c))
    throw PhotometryError("Format error for ConstantMap::create() name " + name
			  + " on input line: " + buffer);
  ConstantMap* cm =  new ConstantMap(c, name);
  return cm;
}

void
ConstantMap::write(std::ostream& os, int precision) const {
  // Save precision and format of stream before changing it:
  StreamSaver ss(os);
  os.precision(precision);
  os.setf( ios_base::showpos | ios_base::scientific);
  os << c << endl;
}
