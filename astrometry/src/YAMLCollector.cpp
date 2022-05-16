#include "YAMLCollector.h"
#include "StringStuff.h"

#ifdef USE_YAML  // None of this is built without YAML

using namespace astrometry;
using namespace stringstuff;

string
YAMLCollector::CAL_PATH = "CAL_PATH";

YAMLCollector::YAMLCollector(string specs, string magic_): magic(magic_) {
  // Split specs at commas
  auto speclist = split(specs, ',');
  for (auto spec : speclist) {
    // Allow empty inputs:
    if (spec.size() == 0) {
      continue;
    }
    // Split each spec at '@' sign
    auto regname = split(spec,'@');
    string expr = "";  // Default will be to match anything
    string filename;
    if (regname.size()==1) {
      filename = regname.front();
    } else if (regname.size()==2) {
      expr = regname.front();
      filename = regname.back();
    } else {
      throw AstrometryError("Bad YAMLCollector regex/name spec: " + spec);
    }
    // Find filename using CAL_PATH environment variable
    string filepath = stringstuff::findFileOnPath(filename, CAL_PATH);
    if (filepath.empty())
      throw AstrometryError("YAMLCollector did not find input YAML file "
			    + filename);
    std::ifstream ifs(filepath.c_str());
    /**/cerr << "Ready to read YAML input " << filepath << endl;
    addInput(ifs, expr, false);
    /**/cerr << "...done with " << filepath << endl;
    ifs.close();
  }
  // Add magic word to root node
  out[magic] = "Magic word";
};

void
YAMLCollector::addInput(istream& is, string filter, bool prepend) {
  if (prepend) 
    inputs.push_front(InFile(is,filter));
  else
    inputs.push_back(InFile(is,filter));
  return;
}

YAMLCollector::InFile::InFile(istream& is, string _filter) {
  stripWhite(_filter);
  // Default filter is to match anything
  if (_filter.empty())
    _filter = ".*";
  filter = std::regex(_filter, std::regex::icase); // Will do case-insensitive filtering

  // Read the input YAML - save the name and serialized string of each node.
  auto root = YAML::Load(is);
  for (auto node : root) {
    if (node.first.as<string>() =="WCS") continue; // Not using WCS's here
    std::ostringstream oss;
    oss << node.second;
    (*this)[node.first.as<string>()] = oss.str();
  }
}

string
YAMLCollector::InFile::findName(string name, const Replacer& reps, int& minSubCount) const {
  // If the name does not pass the input filter, we are done:
  string matchedName;

  if (!std::regex_match(name, filter)) return matchedName;
  
  // First just check for name with no substitutions, as it's v fast
  if (this->count(name)) {
    minSubCount = 0;
    return name;
  }

  if (minSubCount > 1) {
    // we have a chance of finding a match with fewer substitutions.  Try all.
    for (auto pr : *this) {
      //    Do substitutions on name (counting them)
      string testname = pr.first;
      int substitutions = reps.process(testname, minSubCount);
      // Is this a match with the fewest subs so far?
      if ( (substitutions < minSubCount) &&
            stringstuff::nocaseEqual(name, testname)) {
        minSubCount = substitutions;
        matchedName = pr.first;
        if (minSubCount<=1) {
          // We have found a match, we're not going to beat that
          break;
        }
      }
    }
  }
  return matchedName;
}

YAMLCollector::Replacer::Replacer(const Dictionary& dict_) {
  for (auto i : dict_)
    dict.push_back( std::pair<std::regex,string>(std::regex(i.first, std::regex::icase),
						 i.second) );
}

int
YAMLCollector::Replacer::process(string& src, int maxSubs) const {
  int substitutions = 0;
  for (auto& entry : dict) {
    if (std::regex_search(src, entry.first)) {
      // Replacing from the dictionary
      ++substitutions;
      string temp = std::regex_replace(src, entry.first, entry.second);
      src = temp;
    }
    if (substitutions >= maxSubs) break;
  }
  return substitutions;
}

bool
YAMLCollector::addMap(string name, const Dictionary& dict,
		      double* criticalTime) {
  // If we already have an output node with this name, just return success

  // Single thread only since count() may fail if out() is altered simultaneously?
  int outcount;
#ifdef _OPENMP
#pragma omp critical(yamlout)
#endif
  outcount = out.count(name);
  
  if (outcount) return true;

  string serialized; // Where we'll put the serialized version of map
  Replacer reps(dict); // The regexes used for dictionary substitution

  // First do a quick check for exact matches
  for (auto& i : inputs) {
    if (i.count(name)) {
      serialized = i[name];
      reps.process(serialized);
      break;
    }
  }

  // If there's no exact match, search all files for
  // first & most economical match w/substitutions
  if (serialized.empty()) {
    int minSubCount = 99; // Comfortably above any real sub count
    
    for (auto& i : inputs) {
      // Loop through input files, looking for matching mapname
      string matchedName = i.findName(name, reps, minSubCount);
      if (!matchedName.empty()) {
        // Do variable substitution on the serialized specs
        serialized = i[matchedName];
        reps.process(serialized);
      }
      if (minSubCount<=1) break;  // Not going to do better than 1 sub
    }
  }

  // If we have not found any match at all, signal on output:
  if (serialized.empty())
    return false;

  // Create a node from the serialized output.
  YAML::Node mapNode = YAML::Load(serialized);
  if (!mapNode.IsMap() || !mapNode["Type"])
    throw AstrometryError("YAMLCollector found ill-formed map at " + name);
  
  // Save this node for output - but only one thread should write a node at a time:
#ifdef _OPENMP
#pragma omp critical(yamlout)
#endif
  {
    // If multithreading, it's possible the node will have been created since
    // we checked that it's absent.  That should be harmless as it will be replaced
    // here with a fresh copy.
    Stopwatch timer;
    timer.start();
    out[name] = mapNode;
    timer.stop();
    if (criticalTime)
      *criticalTime += timer;
  }
  if (mapNode["Type"].as<string>() == "Composite") {
    // If this node is a compound map, call recursively for all constituent elements
    if (!mapNode["Elements"] || !mapNode["Elements"].IsSequence())
      throw AstrometryError("YAMLCollector is missing Elements node for map " + name);
    for (auto element : mapNode["Elements"]) {
      if (!addMap(element.as<string>(), dict)) 
	throw AstrometryError("Did not find map for composite element "
			      + element.as<string>());
    }
  }

  return true;
};

bool
YAMLCollector::removeMap(string name) {
  return out.erase(name)>0; // Return true if anything was erased
}

void
YAMLCollector::clearMaps() {
  out.clear();
  return;
}

string
YAMLCollector::dump() const {
  YAML::Emitter os;
  os << YAML::BeginMap
     << YAML:: Key << magic
     << YAML:: Value << "Magic word";
  for (auto m : out)
    os << YAML::Key << m.first << YAML::Value << m.second;
  os << YAML::EndMap;
  return os.c_str();
}

#endif
