// Implement the templated PixelMap adjustment.
#include "PiecewiseMap.h"
#include <sstream>
#include "StringStuff.h"

using namespace astrometry;
using stringstuff::getlineNoComment;

// This constructor initializes array with all zeros
PiecewiseMap::PiecewiseMap(string name_,
			   double argStart_,
			   double argEnd,
			   double argStep_,
			   PiecewiseMap::Coordinate axis_,
			   double xCenter_,
			   double yCenter_): PixelMap(name_), 
					     argStart(argStart_),
					     argStep(argStep_),
					     axis(axis_),
					     xCenter(xCenter_),
					     yCenter(yCenter_) {
  // Table must include ArgEnd, though last point is zero.
  // The 0.0001 is to allow for some roundoff error:
  int argCount = static_cast<int> (std::ceil( (argEnd-argStart)/argStep-0.0001));
  if (argCount <= 2) 
    throw AstrometryError("PiecewiseMap constructed with argCount<=2 at " + name);
  v.resize(argCount);
  v.setZero();
  dvda.resize(argCount);
  dvda.setZero();
}

// This constructor is given initial array values.
PiecewiseMap::PiecewiseMap(string name_,
		   double argStart_,
		   double argStep_,
		   const DVector& values,
		   PiecewiseMap::Coordinate axis_,
		   double xCenter_,
		   double yCenter_): PixelMap(name_), 
				     argStart(argStart_),
				     argStep(argStep_),
				     axis(axis_),
				     xCenter(xCenter_),
				     yCenter(yCenter_),
				     v(values),
				     dvda(values.size(),0.) {
  setup();  // Need to calculate the derivatives and make sure first/last are zero.
}


PiecewiseMap*
PiecewiseMap::duplicate() const {
  // Default copy of all members works just fine, since cache holds all tables
  return new PiecewiseMap(*this);
}

void 
PiecewiseMap::setParams(const DVector& p) {
  Assert(p.size()==v.size()-2);
  v.subVector(1,v.size()-1)  = p;
  setup();
}

DVector
PiecewiseMap::getParams() const {
  return v.subVector(1,v.size()-1);
}

/////////////////////////////////////////////////////////////////////////////////
// The actual calculations
/////////////////////////////////////////////////////////////////////////////////


// Given an argument, get the index of its lower-bound node and the
// weight to give to the upper-bound node
void 
PiecewiseMap::getIndex(double arg, int& index, double& frac) const {
  double nodeCount = (arg-argStart) / argStep;
  index = static_cast<int> (std::floor(nodeCount));
  if (index < 0) {
    index = -1; // Could just leave anything <0 as a sign
    frac = 0.;
  } else if (index >= v.size()-1) {
    index = v.size()-1;
    frac = 0.;
  } else {
    frac = nodeCount - index;
  }
}

void
PiecewiseMap::setup() {
  v[0] = 0.;
  v[v.size()-1] = 0.;
  for (int i=0; i<dvda.size()-1; i++)
    dvda[i] = (v[i+1] - v[i]) / argStep;
  dvda[v.size()-1] = 0.;
}

double
PiecewiseMap::lookup(double arg) const {
  double frac;
  int index;
  getIndex(arg,index, frac);
  if (index<0)
    return 0.;
  else
    return v[index] + frac * dvda[index] * argStep;
}

void
PiecewiseMap::valueSlope(double arg, double& value, double& slope) const {
  // return d(lookup) / d(arg)
  double frac;
  int index;
  getIndex(arg, index, frac);
  if (index < 0) {
    value = 0.;
    slope = 0.;
  } else {
    slope = dvda[index];
    value = v[index] + frac * slope * argStep;
  }
  return;
  // Note that I did not try anything fancy when arg is at a node.
}

double 
PiecewiseMap::inverse(double val) const {
  // Try to invert the map arg -> arg+lookup
  // First the cases that are beyond the (monotonically increasing) function's range:
  if ((val <= argStart) || (val >= argStart + argStep*(v.size()-1)))
    return val;

  // Need to localize the interval in which results bound val.
  int index = static_cast<int> (std::floor( (val - lookup(val) - argStart)/argStep));
  const int MAX_ITERATIONS = 5;
  for (int i = 0; i < MAX_ITERATIONS; i++) {
    if (index < 0 || index > v.size()-2)
      throw AstrometryError("Logic problem in PiecewiseMap::reverse_lookup()");
    double lowerVal = argStart + argStep*index + v[index];
    double deltaVal = argStep + (v[index+1]-v[index]);
    // Compute fraction of the way between ends of interval:
    double f = (val - lowerVal) / deltaVal; 
    if ( f < 0. ) {
      // move index down...
      index += static_cast<int> (std::floor( (val - lowerVal)/argStep ));
    } else if ( f > 1. ) {
      // move index up...
      index += 1 + static_cast<int> (std::floor( (val - (lowerVal+deltaVal))/argStep));
    } else {
      // interpolate and return
      return argStart + argStep*(index+f);
    }
  }
  throw AstrometryError("Too many iterations in PiecewiseMap::inverse()");
}

void
PiecewiseMap::toWorld(double xpix, double ypix,
		      double &xworld, double &yworld,
		      double color) const {
  switch (axis) {
  case X:
    xworld = xpix + lookup(xpix);
    yworld = ypix;
    break;
  case Y:
    xworld = xpix;
    yworld = ypix + lookup(ypix);
    break;
  case R:
    double xr = xpix - xCenter;
    double yr = ypix - yCenter;
    double rin = sqrt(xr*xr+yr*yr);
    if (rin == 0.) {
      xworld = xpix;
      yworld = ypix;
    } else {
      double dr = lookup(rin);
      xworld = xpix + dr * xr / rin;
      yworld = ypix + dr * yr / rin;
    }
    break;
  }
}

void
PiecewiseMap::toPix(double xworld, double yworld,
		    double &xpix, double &ypix,
		    double color) const {
  switch (axis) {
  case X:
    xpix = inverse(xworld);
    ypix = yworld;
    break;
  case Y:
    xpix = xworld;
    ypix = inverse(yworld);
    break;
  case R:
    double xr = xworld - xCenter;
    double yr = yworld - yCenter;
    double rout = sqrt(xr*xr+yr*yr);
    if (rout == 0.) {
      xpix = xworld;
      ypix = yworld;
    } else {
      double rin = inverse(rout);
      xpix = xCenter + xr * rin / rout;
      ypix = yCenter + yr * rin / rout;
    }
    break;
  }
}

Matrix22
PiecewiseMap::dWorlddPix(double xpix, double ypix, 
			 double color) const {
  Matrix22 derivs;
  derivs(0,0) = 1.;
  derivs(1,1) = 1.;
  derivs(0,1) = 0.;
  derivs(1,0) = 0.;

  double value, slope;
  switch (axis) {
  case X:
    valueSlope(xpix, value, slope);
    derivs(0,0) += slope;
    break;
  case Y:
    valueSlope(ypix, value, slope);
    derivs(1,1) += slope;
    break;
  case R:
    double xr = xpix - xCenter;
    double yr = ypix - yCenter;
    double rin = sqrt(xr*xr+yr*yr);
    if (rin > 0.) {
      valueSlope(rin, value, slope);
      xr /= rin;
      yr /= rin;
      value / rin;
      double t = slope - value;
      derivs(0,0) += value + xr * xr * t;
      derivs(0,1) += xr * yr * t;
      derivs(1,0) += xr * yr * t;
      derivs(1,1) += value + yr * yr *t;
    }
    break;
  }
  return derivs;
}

void 
PiecewiseMap::toWorldDerivs(double xpix, double ypix,
			    double& xworld, double& yworld,
			    DMatrix& derivs,
			    double color) const {
  Assert(derivs.cols()==v.size()-2 && derivs.rows()==2);
  derivs.setZero();

  // Keep in mind here that first and last v elements must be zero and
  // the parameter vector will skip them, so indexing is off by 1 from v.
  double frac;
  int index;
  yworld = ypix;
  xworld = xpix;
  switch (axis) {
  case X:
    getIndex(xpix, index, frac);
    if (index < 0 || index >= v.size()-1)
      return;
    if (index >0) {
      xworld += v[index] * (1-frac);
      derivs(0,index-1) = 1-frac;
    }
    if (index < v.size()-2){
      xworld += frac * v[index+1];
      derivs(0,index) = frac;
    }
    break;
  case Y:
    getIndex(ypix, index, frac);
    if (index < 0 || index >= v.size()-1)
      return;
    if (index >0) {
      yworld += v[index] * (1-frac);
      derivs(1,index-1) = 1-frac;
    }
    if (index < v.size()-2){
      yworld += frac * v[index+1];
      derivs(1,index) = frac;
    }
    break;
  case R:
    double xr = xpix - xCenter;
    double yr = ypix - yCenter;
    double rin = sqrt(xr*xr+yr*yr);
    getIndex(ypix, index, frac);
    if (index < 0 || index >= v.size()-1)
      return;
    double dr = v[index] * (1-frac) + frac * v[index+1];
    xworld += dr * xr / rin;
    yworld += dr * yr / rin;
    if (index > 0) {
      derivs(0,index-1) = (1-frac) * xr / rin;
      derivs(1,index-1) = (1-frac) * yr / rin;
    }
    if (index < v.size()-2) {
      derivs(0,index) = frac * xr / rin;
      derivs(1,index) = frac * xr / rin;
    }
    break;
  }
}

void 
PiecewiseMap::toPixDerivs(double xworld, double yworld,
			  double &xpix, double &ypix,
			  DMatrix& derivs,
			  double color) const {
  double xdum, ydum;
  // Solve inverse, then just get derivs forward
  toPix(xworld,yworld,xpix, ypix);
  toWorldDerivs(xpix, ypix, xdum, ydum, derivs);
}

////////////////////////////////////////////////////////////////////////////////////
// YAML (de-)serialization
////////////////////////////////////////////////////////////////////////////////////
#ifdef USE_YAML

void
PiecewiseMap::write(YAML::Emitter& os) const {
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
    throw AstrometryError("Unknown PiecewiseMap axis");
  }
  vector<double> vv(v.size());
  for (int i=0; i<vv.size(); i++) vv[i] = v[i];
  os << YAML::Key << "ArgStart" << YAML::Value << argStart
     << YAML::Key << "ArgStep" << YAML::Value << argStep
     << YAML::Key << "Parameters" << YAML::Flow << YAML::Value << vv
     << YAML::EndMap;
}

PixelMap*
PiecewiseMap::create(const YAML::Node& node,
		     bool& defaulted,
		     string name) {
  // Construct map from a YAML entry.
  if (!node.IsMap() ||
      !node["Type"] || node["Type"].as<string>() != type() ||
      !node["Type"] ||
      !node["Axis"] ||
      !node["ArgStep"] || !node["ArgStart"])
    throw AstrometryError("PiecewiseMap::create is missing YAML keys at: " + name);

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
      throw AstrometryError("Missing YAML [XY]Center key for radial PiecewiseMap at " + name);
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

    defaulted = false;
    return new PiecewiseMap(name,
			    argStart, argStep,
			    v,
			    axis,
			    xCenter, yCenter);
  } else if (node["ArgEnd"]) {
    // Build a default (all zeros) LUT of chosen length.
    double argEnd = node["ArgEnd"].as<double>();
    defaulted = true;
    return new PiecewiseMap(name,
			    argStart, argEnd,
			    argStep,
			    axis,
			    xCenter, yCenter);
  } else {
    throw AstrometryError("PiecewiseMap needs either Parameters or ArgEnd in YAML input");
  }
}

#endif
