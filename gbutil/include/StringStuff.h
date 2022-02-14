// Common convenience functions with strings
#ifndef STRINGSTUFF_H
#define STRINGSTUFF_H

#include <iostream>
#include <sstream>
#include <string>
#include <set>
#include <list>
#include <ios>

// For regular-expression matching:
#include <regex>

using std::string;

namespace stringstuff {
  extern bool isComment(const string& instr);

  extern std::istream& getlineNoComment(std::istream& is, string& s);

  bool nocaseEqual(char c1, char c2);
  bool nocaseEqual(const string& s1, const string& s2);

  // Comparison function suitable for standard containers
  class LessNoCase {
  public:
    bool operator()(const string& s1, const string& s2) const;
  };
    
  // Remove trailing whitespace
  void stripTrailingBlanks(string& s);
  // Strip leading and trailing white space from a string:
  void stripWhite(string& s);

  // Remove everything after & including the last "." period
  void stripExtension(string& s);

  // Split string at designated character;
  // Default is to split at whitespace.
  // Trailing whitespace is dropped but trailing
  // non-white will add null string at end.
  std::list<string> split(const string& s, char c=0);

  // Make a string that holds the current time and the command line
  string taggedCommandLine(int argc, char *argv[]);

  // Basic regular expression match (wrapping Boost).  Must match entire input string.
  // Using POSIX extended regular expressions
  bool regexMatch(const std::string& regex_,
		  const string& s,
		  bool caseSensitive=false);

  // Match string to any of several regular expressions
  template <class Container>
  bool regexMatchAny(const Container& regexes,
		     const string& s,
		     bool caseSensitive=false) {
    for (auto& r : regexes) 
      if (regexMatch(r, s)) return true;
    return false;
  }

  // Basic regular expression search and replace (wrapping Boost).  
  // Using POSIX extended regular expression and sed replacement string syntax.
  string regexReplace(const std::string& regex,
		      const std::string& replacement,
		      const string& s);

  // Return a set of all strings in container that match a regular expression.
  // Use POSIX extended definition of regex.  Use Boost Regex package.
  template <class Container>
  std::set<std::string> findMatches(const std::string& regex_,
				    const Container& c,
				    bool caseSensitive=false) {
    std::set<std::string> result;
    auto flags = std::regex::basic;
    if (!caseSensitive) flags |= std::regex::icase;
    std::regex e(regex_, flags);
    for (auto& s : c)
      if (std::regex_match(s, e)) result.insert(s);
    return result;
  }

  // Return a set of all strings in container that match any of multiple regular expressions.
  template <class Container, class Container2>
  std::set<std::string> findMatches(const Container2& regexes,
				    const Container& c,
				    bool caseSensitive=false) {
    std::set<std::string> result;
    for (auto& s : regexes) {
      auto r = findMatches(s, c, caseSensitive);
      result.insert(r.begin(), r.end());
    }
    return result;
  }

  // Class containing a list of regular expressions and replacement strings.
  class RegexReplacements {
  public:
    // Replace s according to rule of any matching regex.  Returns true if any match was found.
    bool operator()(string& s) const; 
    void addRegex(string regex, string replacement);
    bool empty() {return regexes.empty();}
  private:
    std::list<string> regexes;
    std::list<string> replacements;
  };

  // This class saves and clears formatting state of a stream on creation, restores on destruction
  class StreamSaver {
  public:
    StreamSaver(std::ostream& os_): os(os_) {
      oldFormat = os.flags(std::ios_base::fmtflags());
      oldPrecision = os.precision();
    }
    ~StreamSaver() {
      os.flags(oldFormat);
      os.precision(oldPrecision);
    }
  private:
    std::ostream& os;
    std::ios_base::fmtflags oldFormat;
    std::streamsize oldPrecision;
  };

  // This wraps the POSIX C glob function to provide matching filenames
  // to an input pattern.  More than one pattern may be in the string,
  // separated by ",".  Leading/trailing white space are stripped from
  // patterns.
  extern std::list<string>
  fileGlob(const string patterns);

  // This finds first file that is open-able following path given by
  // some environment variable.  Path is just current directory if
  // not given in environment.
  // Returns path of the file, null if none is found.
  // Path is ignored if the filename is absolute.
  // Leading & trailing white space are stripped from filename.
  extern string
  findFileOnPath(const string& filename, const string& environmentVar);
  
} // namespace stringstuff

using namespace stringstuff;

#endif
