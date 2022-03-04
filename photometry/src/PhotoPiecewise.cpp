// Implement the templated PhotoMap adjustment.
#include "PhotoPiecewise.h"
#include <sstream>
#include "StringStuff.h"

using namespace photometry;
using stringstuff::getlineNoComment;

// This constructor initializes array with all zeros
PhotoPiecewise::PhotoPiecewise(string name_,
			       double argStart_,
			       double argEnd,
			       double argStep_,
			       PhotoPiecewise::Coordinate axis_,
			       double xCenter_,
			       double yCenter_): PhotoMap(name_), 
						 argStart(argStart_),
						 argStep(argStep_),
						 axis(axis_),
						 xCenter(xCenter_),
						 yCenter(yCenter_) {
  // Table must include ArgEnd, though last point is zero.
  // The 0.0001 is to allow for some roundoff error:
  int argCount = static_cast<int> (std::ceil( (argEnd-argStart)/argStep-0.0001));
  if (argCount <= 2) 
    throw PhotometryError("PhotoPiecewise constructed with argCount<=2 at " + name);
  v.resize(argCount);
  v.setZero();
}

// This constructor is given initial array values.
PhotoPiecewise::PhotoPiecewise(string name_,
			       double argStart_,
			       double argStep_,
			       const DVector& values,
			       PhotoPiecewise::Coordinate axis_,
			       double xCenter_,
			       double yCenter_): PhotoMap(name_), 
						 argStart(argStart_),
						 argStep(argStep_),
						 axis(axis_),
						 xCenter(xCenter_),
						 yCenter(yCenter_),
						 v(values) {
  if (v.size() <= 2) 
    throw PhotometryError("PhotoPiecewise constructed with argCount<=2 at " + name);
  // Force endpoints to zero
  v[0] = 0.;
  v[v.size()-1] = 0.;
}


PhotoPiecewise*
PhotoPiecewise::duplicate() const {
  // Default copy of all members works just fine, since cache holds all tables
  return new PhotoPiecewise(*this);
}

void 
PhotoPiecewise::setParams(const DVector& p) {
  Assert(p.size()==v.size()-2);
  v.subVector(1,v.size()-1)  = p;
}

DVector
PhotoPiecewise::getParams() const {
  return v.subVector(1,v.size()-1);
}

/////////////////////////////////////////////////////////////////////////////////
// The actual calculations
/////////////////////////////////////////////////////////////////////////////////

// Extract the appropriate x, y, or r value.
double
PhotoPiecewise::getarg(const PhotoArguments& args) const {
  double x,y;
  switch (axis) {
  case X:
    return args.xDevice;
  case Y:
    return args.yDevice;
  case R:
    x = args.xDevice - xCenter;
    y = args.yDevice - yCenter;
    return sqrt(x*x+y*y);
  default:
    throw PhotometryError("Invalid axis in lookup");
  }
}

double
PhotoPiecewise::lookup(double arg) const {
  double nodeCount = (arg-argStart) / argStep;
  int index = static_cast<int> (std::floor(nodeCount));
  if (index < 0) {
    return 0.;	// Forcing to zero out of bounds
  } else if (index >= v.size()-1) {
    return 0.;	// Forcing to zero out of bounds
  } else {
    double frac = nodeCount - index;
    return v[index]*(1-frac) + v[index+1]*frac;
    // Note that I did not try anything fancy when arg is at a node.
  }
}

double
PhotoPiecewise::forwardDerivs(double magIn,
			      const PhotoArguments& args,
			      DVector& derivs) const {
  derivs.resize(v.size()-2);
  derivs.setZero();

  double nodeCount = (getarg(args)-argStart) / argStep;
  int index = static_cast<int> (std::floor(nodeCount));
  if (index < 0) {
    return magIn;
  } else if (index >= v.size()-1) {
    return magIn;
  } else {
    double frac = nodeCount - index;
    if (index>0)
      derivs[index-1] = 1-frac;  // Remember that v[0] is not a free parameter
    if (index<v.size()-2)
      derivs[index] = frac;
    return magIn + v[index]*(1-frac) + v[index+1]*frac;
  }
}

////////////////////////////////////////////////////////////////////////////////////
// YAML (de-)serialization
////////////////////////////////////////////////////////////////////////////////////

void
PhotoPiecewise::write(YAML::Emitter& os) const {
  // Write contents of map to YAML.
  os << YAML::BeginMap
     << YAML::Key << "Type" << YAML::Value << type()
     << YAML::Key << "Axis" << YAML::Value;
  switch (axis) {
  case X:
    os << "X";
    break;
  case Y:
    os << "Y";
    break;
  case R:
    os << "R"
       << YAML::Key << "XCenter" << YAML::Value << xCenter
       << YAML::Key << "YCenter" << YAML::Value << yCenter;
    break;
  default:
    throw PhotometryError("Unknown PhotoPiecewise axis");
  }
  vector<double> vv(v.size());
  for (int i=0; i<vv.size(); i++) vv[i] = v[i];
  os << YAML::Key << "ArgStart" << YAML::Value << argStart
     << YAML::Key << "ArgStep" << YAML::Value << argStep
     << YAML::Key << "Parameters" << YAML::Flow << YAML::Value << vv
     << YAML::EndMap;
}

PhotoMap*
PhotoPiecewise::create(const YAML::Node& node, string name) {
  // Construct map from a YAML entry.
  if (!node.IsMap() ||
      !node["Type"] || node["Type"].as<string>() != type() ||
      !node["Type"] ||
      !node["Axis"] ||
      !node["ArgStep"] || !node["ArgStart"])
    throw PhotometryError("PhotoPiecewise::create is missing YAML keys at: " + name);

  string ax = node["Axis"].as<string>();
  double xCenter, yCenter;
  double argStart, argStep;
  Coordinate axis;
  if (stringstuff::nocaseEqual(ax,"X")) {
    axis = X;
  } else if (stringstuff::nocaseEqual(ax,"Y")) {
    axis = Y;
  } else if (stringstuff::nocaseEqual(ax,"R")) {
    axis = R;
    if (!node["XCenter"] || !node["YCenter"])
      throw PhotometryError("Missing YAML [XY]Center key for radial PhotoPiecewise at " + name);
    xCenter = node["XCenter"].as<double>();
    yCenter = node["YCenter"].as<double>();
  }
  argStart = node["ArgStart"].as<double>();
  argStep = node["ArgStep"].as<double>();

  if (node["Parameters"]) {
    // Build structure from given node values.
    vector<double> vv;
    vv = node["Parameters"].as<vector<double> >();
    DVector v(vv.size());
    for (int i=0; i<v.size(); i++) v[i] = vv[i];

    return new PhotoPiecewise(name,
			      argStart, argStep,
			      v,
			      axis,
			      xCenter, yCenter);
  } else if (node["ArgEnd"]) {
    // Build a default (all zeros) LUT of chosen length.
    double argEnd = node["ArgEnd"].as<double>();
    return new PhotoPiecewise(name,
			      argStart, argEnd,
			      argStep,
			      axis,
			      xCenter, yCenter);
  } else {
    throw PhotometryError("PhotoPiecewise needs either Parameters or ArgEnd in YAML input");
  }
}

