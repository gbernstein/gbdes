#include "StringStuff.h"
#include <cctype>
#include <ctime>
#include <glob.h>
#include <list>
#include <cstdlib>
#include <fstream>

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
      if (s[subEnd]==c || (c==0 && std::isspace(s[subEnd]))) {
	// Do not save a word on consecutive whitespace
	if (c!=0 || subEnd>subStart) 
	  out.push_back(s.substr(subStart, subEnd-subStart));
	subEnd++;
	subStart = subEnd;
      } else {
	subEnd++;
      }
    }
    if (c==0 && subStart==s.size()) {
      // don't add trailing whitespace
      return out;
    } else if (subStart < s.size()) {
      out.push_back(s.substr(subStart));
    } else {
      // Return empty string after non-white
      // trailing separator
      out.push_back("");
    }
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
    auto flags = std::regex::extended;
    if (!caseSensitive) flags |= std::regex::icase;
    std::regex e(regex_, flags); 
    return std::regex_match(s, e);
  }

  string
  regexReplace(const string& regex_,
	       const string& replacement,
	       const string& s) {
  std::regex e(regex_, std::regex::extended);
  return std::regex_replace(s, e, replacement,
			    std::regex_constants::match_default | std::regex_constants::format_sed);
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
      std::regex e(*ireg, std::regex::extended);
      if (std::regex_match(s, e)) {
	s = regex_replace(s, e, *irep,
			  std::regex_constants::match_default | std::regex_constants::format_sed);
	return true;
      }
    }
    return false;  // string is unchanged if there are no matches.
  }

  std::list<string>
  fileGlob(const string patterns) {
    std::list<string> out;
    std::list<string> plist = split(patterns,',');
    for (std::list<string>::iterator i=plist.begin();
	 i != plist.end();
	 ++i) {
      stripWhite(*i);
      glob_t gt;
      glob(i->c_str(), 0, NULL, &gt);
      for (int j=0; j<gt.gl_pathc; j++)
	out.push_back( gt.gl_pathv[j] );
      globfree(&gt);
    }
    return out;
  }


  string
  findFileOnPath(const string& filename, const string& environmentVar) {
    string filepath = filename;
    stripWhite(filepath);
    // Return input if it's empty or already an absolute path
    if (filepath.empty() || filepath[0]=='/') return filepath;
      
    std::list<string> paths;
    // Is there an environment variable giving paths?
    auto pathspec = std::getenv(environmentVar.c_str());
    if (pathspec==0) {
      // Search only current directory
      paths.push_back(".");
    } else {
      paths = stringstuff::split(pathspec, ':');
    }

    for (auto i : paths) {
      // Try each path in turn
      filepath = i + "/" + filename;
      std::ifstream ifs(filepath.c_str());
      if (ifs.good())
	break;
      else
	filepath.clear();
    }
    return filepath;
  }
}
