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

void
PolyMap::setDomain(Bounds<double> domain) {
  // Will be mapping the bounds into (-1,1) intervals
  if (domain.area()==0.)
    FormatAndThrow<PhotometryError>() << "Degenerate domain for PolyMap: " << domain;
  xshift = 0.5*(domain.getXMax() + domain.getXMin());
  xscale = 2. /(domain.getXMax() - domain.getXMin());
  yshift = 0.5*(domain.getYMax() + domain.getYMin());
  yscale = 2. /(domain.getYMax() - domain.getYMin());
}

Bounds<double>
PolyMap::getDomain() const {
  return Bounds<double>( xshift - 1./xscale,
			 xshift + 1./xscale,
			 yshift - 1./yscale,
			 yshift + 1./yscale);
}

double
PolyMap::forward(double magIn, const PhotoArguments& args) const {
  double xin, yin, dMag;
  if (usePixelCoords) {
    xin = args.xDevice;
    yin = args.yDevice;
  } else {
    xin = args.xExposure;
    yin = args.yExposure;
  }
  rescale(xin, yin);
  dMag = poly(xin, yin);
  return magIn + dMag;
}

double
PolyMap::forwardDerivs(double magIn, 
		       const PhotoArguments& args,
		       DVector& derivs) const {
  double xin, yin;
  if (usePixelCoords) {
    xin = args.xDevice;
    yin = args.yDevice;
  } else {
    xin = args.xExposure;
    yin = args.yExposure;
  }
  rescale(xin, yin);
  derivs = poly.derivC(xin, yin);
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
		 Bounds<double> domain): PhotoMap(name), poly(p) {
  setDomain(domain);
  setArgumentType(argType);
}

PolyMap::PolyMap(int orderx, int ordery,
		 ArgumentType argType,
		 string name,
		 Bounds<double> domain): PhotoMap(name), poly(orderx,ordery) {
  setDomain(domain);
  setArgumentType(argType);
}

PolyMap::PolyMap(int order, 
		 ArgumentType argType,
		 string name,
		 Bounds<double> domain): PhotoMap(name), poly(order) {
  setDomain(domain);
  setArgumentType(argType);
}

void
PolyMap::setArgumentType(ArgumentType argType) {
  switch (argType) {
  case Device:
    usePixelCoords = true;
    break;
  case Exposure:
    usePixelCoords = false;
    break;
  default:
    throw PhotometryError("Unknown ArgumentType");
  }
}

////////////////////////////////////////////////////////////////////////
// YAML (de-)serializations
////////////////////////////////////////////////////////////////////////
void
PolyMap::write(YAML::Emitter& os) const {
  Bounds<double> b = getDomain();
  os << YAML::BeginMap 
     << YAML::Key << "Type" << YAML::Value << type()
     << YAML::Key << "UsePixelCoords" << YAML::Value << usePixelCoords
     << YAML::Key << "XMin" << YAML::Value << b.getXMin()
     << YAML::Key << "XMax" << YAML::Value << b.getXMax()
     << YAML::Key << "YMin" << YAML::Value << b.getYMin()
     << YAML::Key << "YMax" << YAML::Value << b.getYMax()
     << YAML::Key << "Poly";
  poly.write(os);
  os << YAML::EndMap;
}

PhotoMap*
PolyMap::create(const YAML::Node& node, string name) {
  if (!node.IsMap() || 
      !node["Type"] || node["Type"].as<string>() != type() ||
      !node["UsePixelCoords"] || !node["Poly"]) {
    cerr << (node["Type"].as<string>() != type())
	 << " " << !node["UsePixelCoords"]
	 << " " << !node["Poly"]
	 << endl;
    throw PhotometryError("PolyMap::create() is missing YAML keys for name " + name);
  }
  ArgumentType argType = node["UsePixelCoords"].as<bool>() ? Device : Exposure;

  // Set up domain, starting with default
  Bounds<double> b(-1.,1.,-1.,1.);
  if (node["XMin"]) b.setXMin(node["XMin"].as<double>());
  if (node["XMax"]) b.setXMax(node["XMax"].as<double>());
  if (node["YMin"]) b.setYMin(node["YMin"].as<double>());
  if (node["YMax"]) b.setYMax(node["YMax"].as<double>());

  bool defaulted;
  poly2d::Poly2d* p = poly2d::Poly2d::create(node["Poly"], defaulted);
  PolyMap* pm =  new PolyMap(*p, argType, name, b);
  delete p;
  return pm;
}

void
ConstantMap::write(YAML::Emitter& os) const {
  // Save precision and format of stream before changing it:
  os << YAML::BeginMap 
     << YAML::Key << "Type" << YAML::Value << type()
     << YAML::Key << "Parameter"
     << YAML::Value << c
     << YAML::EndMap;
}

PhotoMap*
ConstantMap::create(const YAML::Node& node, string name) {
  if (!node.IsMap() || 
      !node["Type"] || node["Type"].as<string>() != type())
    throw PhotometryError("ConstantMap::create() is missing YAML keys for name " + name);
  double c = 0.;  // Constant defaults to zero if not given
  if (node["Parameter"])
    c = node["Parameter"].as<double>();
  ConstantMap* cm =  new ConstantMap(c, name);
  return cm;
}

void
ColorTerm::write(YAML::Emitter& os) const {
  os << YAML::BeginMap
     << YAML::Key << "Type" << YAML::Value << ColorTerm::type()
     << YAML::Key << "Reference" << YAML::Value << refColor()
     << YAML::Key << "Function" << YAML::Value;
  pm->write(os);
  os << YAML::EndMap;
}
  
