// $Id: StringStuff.h,v 1.1.1.1 2009/10/30 21:20:52 garyb Exp $
// Common convenience functions with strings
#ifndef STRINGSTUFF_H
#define STRINGSTUFF_H

#include <iostream>
#include <sstream>
#include <string>
#include <set>
#include <list>

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
  // Strip leading and trailing white space from a string:
  void stripWhite(string& s);

  // Remove everything after & including the last "." period
  void stripExtension(string& s);

  // Split string at designated character
  std::list<string> split(string s, char c);

  // Make a string that holds the current time and the command line
  string taggedCommandLine(int argc, char *argv[]);

  // Basic regular expression match (wrapping Boost)
  bool regexMatch(const std::string& regex_,
		  const string& s,
		  bool caseSensitive=false) {
    boost::regex e(regex_, boost::regex::basic | (caseSensitive ? 0 : boost::regex::icase));
    return boost::regex_match(s, e);
  }

  // Match string to any of several regular expressions
  template <class Container>
  bool regexMatchAny(const Container& regexes,
		     const string& s,
		     bool caseSensitive=false) {
    for (typename Container::const_iterator i = regexes.begin();
	   i != regexes.end();
	   ++i)
      if (regexMatch(*i, s)) return true;
    return false;
  }

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
      std::set<std::string> r = findMatches(*i, c, caseSensitive);
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

} // namespace stringstuff

using namespace stringstuff;

#endif
