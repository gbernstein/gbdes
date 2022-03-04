// Class for building new PixelMap (or PhotoMap) by selecting maps from existing
// serialized collections. 
#ifndef YAMLCOLLECTOR_H
#define YAMLCOLLECTOR_H
#ifdef USE_YAML // Only have this class with YAML

#include "Std.h"
#include "Astrometry.h"
#include "yaml-cpp/yaml.h"
#include <regex>
#include <map>
#include <list>
#include "Stopwatch.h" /**/
namespace astrometry {

  class YAMLCollector {
  public:
    typedef std::map<string,string> Dictionary;

    // Create the class with a string specifying the source files from which
    // serialized maps will be taken.  This string has the format
    // [regex@]file1 [, [regex@]file2], ...
    // where fileN gives the name of a YAML file, and if a regex is given is
    // specifies the map names that can be taken from that file.
    // Regexes and filenames cannot have '@' or ',' characters in them.
    // Leading/trainling whitespace will be stripped.
    // The CAL_PATH environment variable, if present, will specify search path for
    // files.
    // The output YAML map will be initialized with a node with key=magic
    YAMLCollector(string specs, string magic);

    // Add input YAML, either at front or back of priority list.
    // YAML will be decoded from stream.  Only input names matching
    // the regex filter will read from this YAML (default is anything matches).
    // The prepend argument specifies whether to add to front or back of priority
    void addInput(istream& is, string filter="", bool prepend=false);

    // Search the input YAML files' root maps for a Node matching the given name.
    // If found, return true, and add the selected Node to internal YAML structure.
    // If the named Node is of type "Composite" then we will recursively add
    // each map named in the "Elements" sequence.
    // All input files are searched for a match.  The one with fewest substitutions
    // to match the requested name gets priority.  At same priority,
    // the one found in the earliest input file gets priority.  Priority within
    // a file of names with equal substitution count is indeterminate.
    // Substitution and name matching are case-insensitive of
    // of the keys in the Dictionary.
    // Returns false if no match is found.
    // Adds time in non-multi-threaded parts of loop to criticalTime if it's given
    bool addMap(string name, const Dictionary& dict=Dictionary(),
		double* criticalTime=nullptr);

    // Remove node of this name:
    bool removeMap(string name);

    // Remove all maps
    void clearMaps();
    
    // Return serialized version of all the maps that have been added.
    string dump() const;

    // If an environment variable matching this path is defined, it will
    // be used as colon-separated list of paths in which to search for YAML files.
    static string CAL_PATH;

  private:
    // Note that the yaml-cpp code seems very slow for large maps, the indexing
    // time seems to grow linearly with number of elements.
    // So I am storing the nodes in a C++ std::map indexed by string (whereas YAML
    // has to accomodate any type as a key).  On output, put the nodes into one
    // big YAML map.
    std::map<string, YAML::Node> out;
    class Replacer {
      // Class that will replace strings found in input string
      // Initialized with a Dictionary (map) where key/value are the
      // target/replacement strings
    public:
      Replacer(const Dictionary& dict_);
      // Do all replacements in-place on src string.  Return number
      // of distinct words replaced. 
      // Quit subs when number reaches maxSubs.
      int process(string& src, int maxSubs=99) const;
    private:
      std::list<std::pair<std::regex,string>> dict;
    };
    class InFile : public map<string,string> {
      // Class holding input from one file.  Will be a pair of strings: name,
      // and serialized content of the node.  Augment with a regex that 
    public:
      // Build the map from an input YAML file
      InFile(istream& is, string _filter="");
      // Return the name that matches input name with <minSubCount replacements.
      // Returns empty string if none.
      string findName(string name, const Replacer& reps, int& minSubCount) const;
    private:
      std::regex filter;
    };
    std::list<InFile> inputs;
    string magic;  // Magic key that tells type of collection
  };
    
} // end namespace astrometry


#endif  // USE_YAML
#endif  // YAMLCOLLECTOR_H
