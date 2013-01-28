// $Id: StringStuff.cpp,v 1.1.1.1 2009/10/30 21:20:52 garyb Exp $
#include "StringStuff.h"
#include <cctype>
#include <ctime>

namespace stringstuff {
  bool isComment(const string& instr) {
    std::istringstream is(instr);
    string word1;
    is >> word1;
    return (!is || word1.empty() || word1[0]=='#');
  }

  std::istream& getlineNoComment(std::istream& is, string& s) {
    do {
      if (!getline(is,s)) return is;
    } while (isComment(s));
    return is;
  }

  bool nocaseEqual(char c1, char c2) {
    return std::toupper(c1)==std::toupper(c2); 
  }

  bool nocaseEqual(const string& s1, const string& s2) {
    if (s1.size() != s2.size()) return false;
    string::const_iterator p1=s1.begin();
    string::const_iterator p2=s2.begin();
    for ( ; p1!=s1.end(); ++p1, ++p2)
      if (!nocaseEqual(*p1, *p2)) return false;
    return true;
  }

  bool 
  LessNoCase::operator()(const string& s1, const string& s2) const {
    size_t nMin = std::min(s1.size(), s2.size());
    string::const_iterator p1=s1.begin();
    string::const_iterator p2=s2.begin();
    for (size_t i=0; i<nMin; ++i, ++p1, ++p2) {
      if (std::toupper(*p1) < std::toupper(*p2)) return true;
      if (std::toupper(*p1) > std::toupper(*p2)) return false;
    }
    // All characters tied, shorter string is lesser:
    return s1.size() < s2.size();
  }

  void stripTrailingBlanks(string& s) {
    string::iterator tail;
    while (!s.empty()) {
      tail=s.end()-1;
      if (!std::isspace(*tail)) return;
      s.erase(tail);
    }
  }

  // Strip leading and trailing white space from a string:
  void stripWhite(string& s) {
    while (!s.empty() && std::isspace(s[0]))
      s.erase(0,1);
    while (!s.empty() && std::isspace(s[s.size()-1]))
      s.erase(s.size()-1);
  }

  void stripExtension(string& s) {
    size_t dot=s.find_last_of(".");
    if (dot==string::npos) return;	// No extension
    s.erase(dot);
  }

  std::list<string>
  split(const string& s, char c) {
    std::list<string> out;
    size_t subStart=0;
    size_t subEnd=0;
    while (subEnd < s.size()) {
      if (s[subEnd]==c) {
	out.push_back(s.substr(subStart, subEnd-subStart));
	subEnd++;
	subStart = subEnd;
      } else {
	subEnd++;
      }
    }
    if (subStart < s.size()) 
      out.push_back(s.substr(subStart));
    else
      out.push_back("");
    return out;
  }

  std::string taggedCommandLine(int argc, char *argv[]) {
    time_t now;
    time(&now);
    string output = ctime(&now);
    // get rid of last character, which is linefeed:
    output.erase(output.size()-1);  // ????
    output += ": ";
    for (int i=0; i<argc; i++) {
      output += " "; output += argv[i];
    }
    return output;
  }

  bool 
  regexMatch(const std::string& regex_,
	     const string& s,
	     bool caseSensitive) {
    boost::regex e(regex_, boost::regex::basic | (caseSensitive ? 0 : boost::regex::icase));
    return boost::regex_match(s, e);
  }

  void
  RegexReplacements::addRegex(string regex, string replacement) {
    regexes.push_back(regex);
    replacements.push_back(replacement);
  }

  bool
  RegexReplacements::operator()(string& s) const {
    // This is not the fastest way to do this:
    std::list<string>::const_iterator ireg = regexes.begin();
    std::list<string>::const_iterator irep = replacements.begin();
    for ( ; ireg != regexes.end(); ++ireg, ++irep) {
      boost::regex e(*ireg, boost::regex::extended);
      if (boost::regex_match(s, e)) {
	s = regex_replace(s, e, *irep, boost::match_default | boost::format_sed);
	return true;
      }
    }
    return false;  // string is unchanged if there are no matches.
  }
}
