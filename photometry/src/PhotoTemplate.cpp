// Radial and linear lookup-table photometric adjustments
#include <fstream>
#include <sstream>

#include "PhotoTemplate.h"
#include "StringStuff.h"

using namespace photometry;
using namespace stringstuff;
using namespace lookup;

const double DefaultScaling = 1.;
const string cachename="photo";

PhotoTemplate::PhotoTemplate(string tableName,
			     string filename_,
			     string name_): PhotoMap(name_), hasSplit(false), xSplit(0.),
					    tableLow(nullptr), tableHigh(nullptr),
					    tableLowName(tableName), tableHighName(""),
					    filename(filename_),
					    scaling(DefaultScaling) {
  // Load the file if not already:
  Lookup1d::ingestFile(filename,cachename);
  tableLow = Lookup1d::find(tableLowName, cachename);
}

// This constructor is to use 2 LUTs astride a split in X values
PhotoTemplate::PhotoTemplate(double xSplit_,
			     string lowName, string highName,
			     string filename_,
			     string name_): PhotoMap(name_), hasSplit(true), xSplit(xSplit_),
					    tableLow(nullptr), tableHigh(nullptr),
					    tableLowName(lowName), tableHighName(highName),
					    filename(filename_),
					    scaling(DefaultScaling) {
  // Load the file if not already:
  Lookup1d::ingestFile(filename,cachename);
  tableLow  = Lookup1d::find(tableLowName, cachename);
  tableHigh = Lookup1d::find(tableHighName, cachename);
}

////////////////////////////////////////////////////////////////////////////////////
// YAML (de-)serialization
////////////////////////////////////////////////////////////////////////////////////

void
PhotoTemplate::write(YAML::Emitter& os) const {
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

PhotoMap*
PhotoTemplate::create(const YAML::Node& node, string name) {
  if (!node.IsMap() ||
      !node["Type"] || node["Type"].as<string>() != type() ||
      !node["HasSplit"] ||
      !node["Filename"] ||
      !node["LowTable"])
    throw PhotometryError("PhotoTemplate::create is missing YAML keys at: " + name);
  bool hasSplit = node["HasSplit"].as<bool>();
  string filename = node["Filename"].as<string>();
  string lowName = node["LowTable"].as<string>();
  PhotoTemplate* out=0;
  if (hasSplit) {
    if (!node["HighTable"] || !node["XSplit"]) 
      throw PhotometryError("Missing YAML keys creating a split PhotoTemplate at " + name);
    string highName = node["HighTable"].as<string>();
    double xSplit = node["XSplit"].as<double>();
    out = new PhotoTemplate(xSplit, lowName, highName, filename, name);
  } else {
    out = new PhotoTemplate(lowName, filename, name);
  }
  if (node["Parameter"]) {
    double scaling = 
    out->scaling = node["Parameter"].as<double>();
  }
  return out;
}


