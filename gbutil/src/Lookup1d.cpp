#include "Lookup1d.h"
#include "StringStuff.h"

using namespace lookup;

// Static variables used for cache

string
Lookup1d::CAL_PATH="CAL_PATH";

std::map<string, Lookup1d::Cache>
Lookup1d::caches;
std::map<string,std::set<string> >
Lookup1d::cachedFilenames;

/////////////////////////////////////////
// (De-)serialization routines
/////////////////////////////////////////
#ifdef USE_YAML

void
Lookup1d::write(YAML::Emitter& os) const {
  // Write contents of LUT to YAML.
  os << YAML::BeginMap
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
    throw LookupError("Unknown Lookup1d.axis");
  }
  os << YAML::Key << "ArgStart" << YAML::Value << argStart
     << YAML::Key << "ArgStep" << YAML::Value << argStep
     << YAML::Key << "Values" << YAML::Flow << YAML::Value << v
     << YAML::EndMap;
}

Lookup1d::Lookup1d(const YAML::Node& node) {
  // Construct LUT from a YAML entry.
  if (!node.IsMap() ||
      !node["Axis"] ||
      !node["ArgStep"] || !node["ArgStart"] ||
      !node["Values"])
    throw LookupError("Missing YAML keys for Lookup1d");

  string ax = node["Axis"].as<string>();
  if (stringstuff::nocaseEqual(ax,"X")) {
    axis = X;
  } else if (stringstuff::nocaseEqual(ax,"Y")) {
    axis = Y;
  } else if (stringstuff::nocaseEqual(ax,"R")) {
    axis = R;
    if (!node["XCenter"] || !node["YCenter"])
      throw LookupError("Missing YAML [XY]Center key for radial Lookup1d");
    xCenter = node["XCenter"].as<double>();
    yCenter = node["YCenter"].as<double>();
  }
  argStart = node["ArgStart"].as<double>();
  argStep = node["ArgStep"].as<double>();
  v = node["Values"].as<vector<double> >();
  // Calculate and keep the slopes:
  dvda.resize(v.size());
  for (int i=0; i<dvda.size()-1; i++)
    dvda[i] = (v[i+1]-v[i])/argStep;

}
#endif 

void
Lookup1d::ingestFile(const string& filename, const string& cachename) {
#ifdef USE_YAML
  // Create new cache if needed
  if (caches.count(cachename)==0) {
    caches[cachename] = Cache();
    cachedFilenames[cachename] = std::set<string>();
  }

  if (cachedFilenames[cachename].count(filename)==0) {
    // Need to ingest this file.  First, find and open it.
    string filepath = findFileOnPath(filename, CAL_PATH);

    if (filepath.empty())
      throw LookupError("Could not open Lookup1d file " + filename);

    // Read file into YAML
    ifstream ifs(filepath.c_str());
    YAML::Node root = YAML::Load(ifs);
    if (!root.IsMap())
      throw LookupError("Lookup1d file " + filepath + " is not a YAML map");

    // Decode each map entry as a Lookup1d and place into cache.
    // Exception for duplicate table names within cache
    for (YAML::const_iterator i = root.begin();
	 i != root.end();
	 ++i) {
      string tablename = i->first.as<string>();
      if (caches[cachename].count(tablename)>0) 
	throw LookupError("Duplicate entry for table " + tablename +
			  " found in file " + filepath);
      caches[cachename][tablename] = new Lookup1d(i->second);
    }
    // Add this filename to set of ingested ones
    cachedFilenames[cachename].insert(filename);
  }
#else
  // No YAML: this can't be done
  cerr << "ERROR: Lookup1d cannot be built and used without YAML libraries" << endl;
  exit(1);
#endif
}

const Lookup1d*
Lookup1d::find(const string& tablename, const string& cachename) {
#ifdef USE_YAML
  if (caches.count(cachename)==0)
    throw LookupError("No such Lookup1d cache name " + cachename);
  auto i = caches[cachename].find(tablename);
  if (i==caches[cachename].end())
    throw LookupError("No Lookup1d with name " + tablename +
		      " in cache " + cachename);
  return i->second;
#else
  // No YAML: this can't be done
  cerr << "ERROR: Lookup1d cannot be built and used without YAML libraries" << endl;
  exit(1);
#endif
}

/////////////////////////////////////////
// Interpolation routines
/////////////////////////////////////////

double
Lookup1d::getArg(double x, double y) const {
  switch (axis) {
  case X:
    return x;
  case Y:
    return y;
  case R:
    x = x - xCenter;
    y = y - yCenter;
    return sqrt(x*x+y*y);
  default:
    throw LookupError("Invalid axis in Lookup1d");
  }
}

void
Lookup1d::valueSlope(double arg, double& val, double& slope) const {
  double nodeCount = (arg-argStart) / argStep;
  int index = static_cast<int> (std::floor(nodeCount));
  // return d(lookup) / d(arg)
  if (index < 0) {
    val = v.front();
    slope = 0.;
  } else if (index >= v.size()-1) {
    val = v.back();
    slope = 0.;
  } else {
    double f = nodeCount - index;
    slope = dvda[index];
    val = v[index] + f * slope * argStep;
  }
  return;
  // Note that I did not try anything fancy for slope when arg is at a node.
}

// Find location where arg + scale*f(arg) = val:
double
Lookup1d::inverse(double val, double scale) const {
  // Try to invert the map arg -> arg+lookup*scale
  // First the cases that are beyond the (monotonically increasing) function's range:
  if (val <= argStart + scale * v.front())
    return val - scale * v.front();
  if (val >= argStart + argStep*(v.size()-1)+ scale * v.back())
    return val - scale * v.back();

  // Need to localize the interval in which results bound val.
  int index = static_cast<int> (std::floor( (val - scale*value(val) - argStart)/argStep));
  const int MAX_ITERATIONS = 5;
  for (int i = 0; i < MAX_ITERATIONS; i++) {
    if (index < 0 || index > v.size()-2)
      throw LookupError("Logic problem in Lookup1d::inverse()");
    double lowerVal = argStart + argStep*index + scale*v[index];
    double deltaVal = argStep + scale*(v[index+1]-v[index]);
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
  throw LookupError("Too many iterations in Lookup1d::inverse()");
}
