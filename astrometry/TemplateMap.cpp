// Implement the templated PixelMap adjustment.
#include "TemplateMap.h"
#include <sstream>
#include "StringStuff.h"

using namespace astrometry;
using stringstuff::getlineNoComment;

const double DefaultPixelStep = 1.;
const double DefaultScaling = 1.;
TemplateMap1d::TemplateMap1d(Coordinate c, 
			     double nodeStart_, double nodeStep_, 
			     const vector<double>& deviations_,
			     string name):  PixelMap(name),
					    nodeStart(nodeStart_), nodeStep(nodeStep_), 
					    deviations(deviations_),
					    pixelStep(DefaultPixelStep),
					    scaling(DefaultScaling)
{
  switch (c) {
  case X:
    applyToX = true;
    break;
  case Y:
    applyToX = false;
    break;
  default:
    throw AstrometryError("Unknown TemplateMap1d::Coordinate value");
  }
}

TemplateMap1d*
TemplateMap1d::duplicate() const {
  return new TemplateMap1d(*this);
}

void 
TemplateMap1d::setParams(const DVector& p) {
  Assert(p.size()==1);
  scaling = p[0];
}

void
TemplateMap1d::toWorld(double xpix, double ypix,
		       double& xworld, double& yworld) const {
  xworld = xpix;
  yworld = ypix;
  if (applyToX)
    xworld += scaling * lookup(xpix);
  else
    yworld += scaling * lookup(ypix);
}

void
TemplateMap1d::toPix(double xworld, double yworld,
		     double& xpix, double& ypix) const {
  if (applyToX) {
    xpix = reverse_lookup(xworld);
    ypix = yworld; 
  } else {
    xpix = xworld; 
    ypix = reverse_lookup(yworld);
  }  
}


void 
TemplateMap1d::toWorldDerivs(double xpix, double ypix,
			     double& xworld, double& yworld,
			     DMatrix& derivs) const {
  Assert(derivs.ncols()==1 && derivs.nrows()==2);
  toWorld(xpix, ypix, xworld, yworld);
  double dev = lookup( applyToX? xpix : ypix);
  if (applyToX) {
    derivs(0,0) = dev;
    derivs(1,0) = 0.;
  } else {
    derivs(0,0) = 0.;
    derivs(1,0) = dev;
  }
}

void 
TemplateMap1d::toPixDerivs( double xworld, double yworld,
			    double &xpix, double &ypix,
			    DMatrix& derivs) const {
  Assert(derivs.ncols()==1 && derivs.nrows()==2);
  toPix(xworld, yworld, xpix, ypix);
  double dev = lookup( applyToX? xpix : ypix);
  if (applyToX) {
    derivs(0,0) = dev;
    derivs(1,0) = 0.;
  } else {
    derivs(0,0) = 0.;
    derivs(1,0) = dev;
  }
}

Matrix22
TemplateMap1d::dWorlddPix(double xpix, double ypix) const {
  Matrix22 m;
  m.setToIdentity();
  if (applyToX) {
    m(0,0) = 1. + scaling * slope(xpix);
  } else {
    m(1,1) = 1. + scaling * slope(ypix);
  }
  return m;
}

Matrix22
TemplateMap1d::dPixdWorld(double xworld, double yworld) const {
  Matrix22 m;
  m.setToIdentity();
  double xpix, ypix;
  toPix(xworld, yworld, xpix, ypix);
  if (applyToX) {
    m(0,0) = 1./(1. + scaling * slope(xpix));
  } else {
    m(1,1) = 1./(1. + scaling * slope(ypix));
  }
  return m;
}

double
TemplateMap1d::lookup(double arg) const {
  double nodeCount = (arg-nodeStart) / nodeStep;
  int index = static_cast<int> (std::floor(nodeCount));
  if (index < 0) {
    return deviations.front();
  } else if (index >= deviations.size()-1) {
    return deviations.back();
  } else {
    double f = nodeCount - index;
    return deviations[index] * (1.-f) + deviations[index+1] * f;
  }
}

double
TemplateMap1d::slope(double arg) const {
  double nodeCount = (arg-nodeStart) / nodeStep;
  int index = static_cast<int> (std::floor(nodeCount));
  // return d(lookup) / d(arg)
  if (index < 0) {
    return 0.;
  } else if (index >= deviations.size()-1) {
    return 0.;
  } else {
    return (deviations[index+1] - deviations[index]) / nodeStep;
  }
  // Note that I did not try anything fancy when arg is at a node.
}

double 
TemplateMap1d::reverse_lookup(double val) const {
  // Try to invert the map arg -> arg+lookup*scaling
  // First the cases that are beyond the (monotonically increasing) function's range:
  if (val <= nodeStart + scaling * deviations.front())
    return val - scaling * deviations.front();
  if (val >= nodeStart + nodeStep*(deviations.size()-1)+ scaling * deviations.back())
    return val - scaling * deviations.back();

  // Need to localize the nodes bounding
  int index = static_cast<int> (std::floor( (val - scaling*lookup(val) - nodeStart)/nodeStep));
  const int MAX_ITERATIONS = 5;
  for (int i = 0; i < MAX_ITERATIONS; i++) {
    if (index < 0 || index > deviations.size()-2)
      throw AstrometryError("Logic problem in TemplateMap1d::reverse_lookup()");
    if ( val - nodeStart > nodeStep*(index+1) + scaling*deviations[index+1]) {
      // move index up...
      index = static_cast<int> (std::floor( (val - scaling*deviations[index+1] - nodeStart)/nodeStep));
    } else if ( val - nodeStart < nodeStep*(index) + scaling*deviations[index]) {
      // move index down...
      index = static_cast<int> (std::floor( (val - scaling*deviations[index] - nodeStart)/nodeStep));
    } else {
      // interpolate and return
      double f = (val - nodeStart - scaling*deviations[index]) /
	(scaling*(deviations[index+1]-deviations[index])+nodeStep);
      return nodeStart + nodeStep*(index+f);
    }
  }
  throw AstrometryError("Too many iterations in TemplateMap1d::reverse_lookup()");
}

void
TemplateMap1d::write(std::ostream& os) const {
  os << (applyToX ? "X" : "Y")
     << " " << deviations.size()
     << " " << nodeStart
     << " " << nodeStep
     << endl;
  os << scaling << endl;
  const int NodesPerLine = 5;
  int nodesThisLine = 0;
  for (int i=0; i<deviations.size(); i++) {
    if (nodesThisLine >= NodesPerLine) {
      os << endl;
      nodesThisLine = 0;
    }
    os << deviations[i] << " ";
    nodesThisLine++;
  }
  os << endl;
}

PixelMap*
TemplateMap1d::create(std::istream& is, string name) {
  string axis;
  int nNodes;
  double start;
  double step;
  double scaling;
  vector<double> devs;
  Coordinate c;

  string buffer;
  if (!getlineNoComment(is, buffer))
    throw AstrometryError("Missing input for Template1d::create for name <" + name + ">");
  {
    istringstream iss(buffer);
    if (! (iss >> axis >> nNodes >> start >> step))
    throw AstrometryError("Bad spec inputs for Template1d::create at line: <" + buffer + ">");
    if (axis=="X" || axis=="x")
      c = X;
    else if (axis=="Y" || axis=="y")
      c = X;
    else
      throw AstrometryError("Input to Template1d::create must be X or Y axis, have <" + axis + "<");
  }
  if (!getlineNoComment(is, buffer))
    throw AstrometryError("Missing scaling input for Template1d::create for name <" + name + ">");
  {
    istringstream iss(buffer);
    if (! (iss >> scaling))
    throw AstrometryError("Bad scaling input for Template1d::create at line: <" + buffer + ">");
  }

  devs.reserve(nNodes);
  while (devs.size() < nNodes) {
    /**/cerr << devs.size() << endl;
    if (!(stringstuff::getlineNoComment(is, buffer)))
      throw AstrometryError("Out of input reading nodes for TemplateMap1d: " + buffer);
    /**/cerr << "read input line <" + buffer + ">" << endl;
    double c;
    istringstream iss(buffer);
    do {
      iss >> c;
      /**/cerr << "read " << c << endl;
      if (iss.fail())
	throw AstrometryError("Error reading nodes for TemplateMap1d: " + buffer);
      devs.push_back(c);
      /**/cerr << "pushed" << endl;
      if (iss.eof()) 
	break;
    } while (devs.size() < nNodes);
  }
  
  TemplateMap1d* out = new TemplateMap1d(c, start, step, devs, name);
  out->scaling = scaling;
  return out;
}

