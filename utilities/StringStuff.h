// $Id: StringStuff.h,v 1.1.1.1 2009/10/30 21:20:52 garyb Exp $
// Common convenience functions with strings
#ifndef STRINGSTUFF_H
#define STRINGSTUFF_H

#include <iostream>
#include <sstream>
#include <string>
#include <set>

// For regular-expression matching:
#include <boost/regex.hpp>

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
  // Remove everything after & including the last "." period
  void stripExtension(string& s);

  // Make a string that holds the current time and the command line
  string taggedCommandLine(int argc, char *argv[]);

  // Return a set of all strings in container that match a regular expression.
  // Use POSIX definition of regex.  Use Boost Regex package.
  template <class Container>
  std::set<std::string> findMatches(const std::string& regex_,
				    const Container& c,
				    bool caseSensitive=false) {
    std::set<std::string> result;
    boost::regex e(regex_, boost::regex::basic | (caseSensitive ? 0 : boost::regex::icase));
    for (typename Container::const_iterator i = c.begin();
	   i != c.end();
	   ++i)
      if (boost::regex_match(*i, e)) result.insert(*i);
    return result;
  }

  // Return a set of all strings in container that match any of multiple regular expressions.
  template <class Container, class Container2>
  std::set<std::string> findMatches(const Container2& regexes,
				    const Container& c,
				    bool caseSensitive=false) {
    std::set<std::string> result;
    for (typename Container2::const_iterator i = regexes.begin();
	 i != regexes.end();
	 ++i) {
      std::set<std::string> r = findMatches(*i, c);
      result.insert(r.begin(), r.end());
    }
    return result;
  }
} // namespace stringstuff

using namespace stringstuff;

#endif
