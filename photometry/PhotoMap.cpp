#include <cctype>
#include "PhotoMap.h"
#include "StringStuff.h"

using namespace photometry;

int 
PhotoMap::anonymousCounter=0;

PhotoMap::PhotoMap(string name_): name(name_) {
  if (name.empty()) {
    std::ostringstream oss;
    oss << "map_" << anonymousCounter++;
    name = oss.str();
  }
}

double
PolyMap::forward(double magIn, const PhotoArguments& args) const {
  double dMag;
  if (useExposureCoords)
    dMag += poly(args.xExposure, args.yExposure);
  else
    dMag += poly(args.xDevice, args.yDevice);

  if (useColor) 
    dMag *= args.color;
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
  if (useColor)
    derivs *= args.color;
  return forward(magIn, args);
}

void
PolyMap::setParams(const DVector& p) {
  Assert(p.size()==nParams());
  poly.setC(p);
}

PolyMap::PolyMap(const poly2d::Poly2d& p,
		 ArgumentType argType,
		 string name,
		 bool isColorTerm): PhotoMap(name), poly(p), useColor(isColorTerm) {
  setArgumentType(argType);
}

PolyMap::PolyMap(int orderx, int ordery,
		 ArgumentType argType,
		 string name,
		 bool isColorTerm): PhotoMap(name), poly(orderx,ordery), useColor(isColorTerm) {
  setArgumentType(argType);
}

PolyMap::PolyMap(int order, 
		 ArgumentType argType,
		 string name,
		 bool isColorTerm): PhotoMap(name), poly(order), useColor(isColorTerm) {
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

  // Use color term?
  string s;
  bool ifColor=false;
  if ( iss >> s ) {
    if (stringstuff::nocaseEqual(s, "color"))
      ifColor = true;
    else 
      throw PhotometryError("PolyMap::create() expected to be told color term but got <" 
			    + s + ">");
  }
  poly2d::Poly2d* p = poly2d::Poly2d::create(is);
  PolyMap* pm =  new PolyMap(*p, argType, name, ifColor);
  delete p;
  return pm;
}

void
PolyMap::write(std::ostream& os) const {
  os << (useExposureCoords ? "Exposure" : "Device");
  if (useColor)
    os << " Color";
  os << endl;
  poly.write(os);
}
