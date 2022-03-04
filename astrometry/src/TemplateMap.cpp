// Implement the templated PixelMap adjustment.
#include "TemplateMap.h"
#include <sstream>
#include "StringStuff.h"

using namespace astrometry;
using stringstuff::getlineNoComment;
using namespace lookup;

const double DefaultScaling = 1.;

const string cachename="astro";

TemplateMap::TemplateMap(string tableName,
			 string filename_,
			 string name_): PixelMap(name_), hasSplit(false), xSplit(0.),
					tableLow(0), tableHigh(0),
					tableLowName(tableName), tableHighName(""),
					filename(filename_),
					scaling(DefaultScaling) {
  // Load the file if not already:
  Lookup1d::ingestFile(filename, cachename);
  tableLow  = Lookup1d::find(tableLowName, cachename);
}

// This constructor is to use 2 LUTs astride a split in X values
TemplateMap::TemplateMap(double xSplit_,
			 string lowName, string highName,
			 string filename_,
			 string name_): PixelMap(name_), hasSplit(true), xSplit(xSplit_),
					tableLow(0), tableHigh(0),
					tableLowName(lowName), tableHighName(highName),
					filename(filename_),
					scaling(DefaultScaling) {
  // Load the file if not already:
  Lookup1d::ingestFile(filename, cachename);
  tableLow  = Lookup1d::find(tableLowName, cachename);
  tableHigh = Lookup1d::find(tableHighName, cachename);
}


TemplateMap*
TemplateMap::duplicate() const {
  // Default copy of all members works just fine, since cache holds all tables
  return new TemplateMap(*this);
}

void 
TemplateMap::setParams(const DVector& p) {
  Assert(p.size()==1);
  scaling = p[0];
}

/////////////////////////////////////////////////////////////////////////////////
// The actual calculations
/////////////////////////////////////////////////////////////////////////////////


void
TemplateMap::toWorld(double xpix, double ypix,
		     double &xworld, double &yworld,
		     double color) const {
  const lookup::Lookup1d* lu= (hasSplit && xpix>xSplit) ? tableHigh : tableLow;
  double delta = scaling*lu->value(xpix,ypix);

  xworld = xpix;
  yworld = ypix;
  switch (lu->getAxis()) {
  case Lookup1d::X:
    xworld += delta;
    return;
  case Lookup1d::Y:
    yworld += delta;
    return;
  case Lookup1d::R:
    double xc, yc;
    lu->getCenter(xc,yc);
    double xr = xpix - xc;
    double yr = ypix - yc;
    double rin = sqrt(xr*xr+yr*yr);
    if (rin > 0.) {
      xworld += delta * xr / rin;
      yworld += delta * yr / rin;
    }
    return;
  }
  throw AstrometryError("Bad Lookup1d axis");
}

void
TemplateMap::toPix(double xworld, double yworld,
		   double &xpix, double &ypix,
		   double color) const {
  // Note that inverse is not handle properly at xSplit???
  const lookup::Lookup1d* lu= (hasSplit && xpix>xSplit) ? tableHigh : tableLow;
  double arg = lu->inverse(lu->getArg(xworld,yworld),scaling);

  switch (lu->getAxis()) {
  case Lookup1d::X:
    xpix = arg;
    ypix = yworld;
    return;
  case Lookup1d::Y:
    xpix = xworld;
    ypix = arg;
    return;
  case Lookup1d::R:
    double xc, yc;
    lu->getCenter(xc,yc);
    double xr = xworld - xc;
    double yr = yworld - yc;
    double rout = sqrt(xr*xr+yr*yr);
    if (rout == 0.) {
      xpix = xworld;
      ypix = yworld;
    } else {
      xpix = xc + arg * xr / rout;
      ypix = yc + arg * yr / rout;
    }
    return;
  }
  throw AstrometryError("Bad Lookup1d axis");
}

void 
TemplateMap::toWorldDerivs(double xpix, double ypix,
			   double& xworld, double& yworld,
			   DMatrix& derivs,
			   double color) const {
  Assert(derivs.cols()==1 && derivs.rows()==2);

  const lookup::Lookup1d* lu= (hasSplit && xpix>xSplit) ? tableHigh : tableLow;
  double delta = lu->value(xpix,ypix);

  xworld = xpix;
  yworld = ypix;
  switch (lu->getAxis()) {
  case Lookup1d::X:
    xworld += scaling*delta;
    derivs(0,0) = delta;
    derivs(1,0) = 0.;
    return;
  case Lookup1d::Y:
    yworld += scaling*delta;
    derivs(0,0) = 0.;
    derivs(1,0) = delta;
    return;
  case Lookup1d::R:
    double xc, yc;
    lu->getCenter(xc,yc);
    double xr = xpix - xc;
    double yr = ypix - yc;
    double rin = sqrt(xr*xr+yr*yr);
    if (rin == 0.) {
      derivs(0,0) = 0.;
      derivs(1,0) = 0.;
    } else {
      xworld += scaling * delta * xr / rin;
      yworld += scaling * delta * yr / rin;
      derivs(0,0) = delta * xr / rin;
      derivs(1,0) = delta * yr / rin;
    }
    return;
  }
  throw AstrometryError("Bad Lookup1d axis");
}

void 
TemplateMap::toPixDerivs(double xworld, double yworld,
			 double &xpix, double &ypix,
			 DMatrix& derivs,
			 double color) const {
  double xdum, ydum;
  // Do this lazily and somewhat wastefully by
  // getting inverse first
  toPix(xworld,yworld,xpix,ypix);
  // and then get derivs from fwd transform
  toWorldDerivs(xpix, ypix, xdum, ydum, derivs);
}

Matrix22
TemplateMap::dWorlddPix(double xpix, double ypix, 
			double color) const {
  Matrix22 derivs;
  derivs(0,0) = 1.;
  derivs(1,1) = 1.;
  derivs(0,1) = 0.;
  derivs(1,0) = 0.;

  const lookup::Lookup1d* lu= (hasSplit && xpix>xSplit) ? tableHigh : tableLow;
  double value, slope;
  lu->valueSlope(xpix,ypix,value,slope);
  
  switch (lu->getAxis()) {
  case Lookup1d::X:
    derivs(0,0) += scaling*slope;
    return derivs;
  case Lookup1d::Y:
    derivs(1,1) += scaling*slope;
    return derivs;
  case Lookup1d::R:
    double xc, yc;
    lu->getCenter(xc,yc);
    double xr = xpix - xc;
    double yr = ypix - yc;
    double rin = sqrt(xr*xr+yr*yr);
    if (rin > 0.) {
      xr /= rin;
      yr /= rin;
      value / rin;
      double t = slope - value;
      derivs(0,0) += scaling*(value + xr * xr * t);
      derivs(0,1) += scaling*(xr * yr * t);
      derivs(1,0) += scaling*(xr * yr * t);
      derivs(1,1) += scaling*(value + yr * yr *t);
    }
    return derivs;
  }
  throw AstrometryError("Bad Lookup1d axis");
}


////////////////////////////////////////////////////////////////////////////////////
// YAML (de-)serialization
////////////////////////////////////////////////////////////////////////////////////
#ifdef USE_YAML

void
TemplateMap::write(YAML::Emitter& os) const {
  os << YAML::BeginMap
     << YAML::Key << "Type" << YAML::Value << type()
     << YAML::Key << "HasSplit" << YAML::Value << hasSplit
     << YAML::Key << "Filename" << YAML::Value << filename
     << YAML::Key << "LowTable" << YAML::Value << tableLowName;
  if (hasSplit) {
    os << YAML::Key << "XSplit" << YAML::Value << xSplit
       << YAML::Key << "HighTable" << YAML::Value << tableHighName;
  }
  os << YAML::Key << "Parameter" << YAML::Value << scaling
     << YAML::EndMap;
}

PixelMap*
TemplateMap::create(const YAML::Node& node,
		    bool& defaulted,
		    string name) {
  if (!node.IsMap() ||
      !node["Type"] || node["Type"].as<string>() != type() ||
      !node["HasSplit"] ||
      !node["Filename"] ||
      !node["LowTable"])
    throw AstrometryError("TemplateMap::create is missing YAML keys at: " + name);
  bool hasSplit = node["HasSplit"].as<bool>();
  string filename = node["Filename"].as<string>();
  double scaling = 1.;
  if (node["Parameter"]) {
    defaulted = false;
    scaling = node["Parameter"].as<double>();
  } else {
    defaulted = true;
  }
  string lowName = node["LowTable"].as<string>();
  TemplateMap* out=0;
  if (hasSplit) {
    if (!node["HighTable"] || !node["XSplit"]) 
      throw AstrometryError("Missing YAML keys creating a split TemplateMap at " + name);
    string highName = node["HighTable"].as<string>();
    double xSplit = node["XSplit"].as<double>();
    out = new TemplateMap(xSplit, lowName, highName, filename, name);
  } else {
    out = new TemplateMap(lowName, filename, name);
  }
  out->scaling = scaling;
  return out;
}

#endif
