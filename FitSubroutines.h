// Constants and subroutines shared by the astrometric and photometric matching/fitting codes
#ifndef FITSUBROUTINES_H
#define FITSUBROUTINES_H

#include <list>
#include <set>
#include "StringStuff.h"
#include "Pset.h"

using namespace std;

const int REF_INSTRUMENT=-1;	// Instrument for reference objects (no fitting)
const int TAG_INSTRUMENT=-2;	// Exposure number for tag objects (no fitting nor contrib to stats)
const int NO_INSTRUMENT=-3;

const string stellarAffinity="STELLAR";

const double NO_MAG_DATA = -100.;	// Value entered when there is no valid mag or color
const string colorColumnName = "COLOR";
const string colorErrorColumnName = "COLOR_ERR";

// Here is the default character at which to split lists given in parameters strings
const char DefaultListSeperator=',';

// A helper function that strips white space from front/back of a string and replaces
// internal white space with underscores:
void
spaceReplace(string& s);

// Another helper function to split up a string into a list of whitespace-trimmed strings.
// Get rid of any null ones.
list<string> 
splitArgument(string input, const char listSeperator = DefaultListSeperator);

// Read parameters from files on command line and from std input.
// First nRequiredArgs must be present and are not names of parameter files.
// Print the usage string (to cerr) and the default parameters (to cout) if
// there are not enough parameters.
// Then dumps final parameter values to cout.
void  processParameters(Pset& parameters, 
			const string& usage,
			int nRequiredArgs,
			int argc,
			char *argv[]);

// Take a string with this format:
// <regex>=<replace> [, <regex>=<replace>] ...
// and parse it into a stringstuff::RegexReplacement object that will execute
// the regular expression replacements on an input string.
stringstuff::RegexReplacements  parseTranslator(string specString, string errorDescription);

// This class reads a set of exposure/object number pairs from a file on construction.
// Then it can be used as a function that returns true if a given pair was on the list.
// First it needs the class holding the pairs (user does not need this)

class EOPair {
public:
  EOPair(int exp, long obj): extensionNumber(exp), objectNumber(obj) {}
  bool operator<(const EOPair& rhs) const {
    return extensionNumber < rhs.extensionNumber
      || (extensionNumber==rhs.extensionNumber && objectNumber < rhs.objectNumber);
  }
private:
  int extensionNumber;
  long objectNumber;
};

class ExtensionObjectSet {
public:
  ExtensionObjectSet(string filename);
  bool operator()(int extensionNumber, long objectNumber) const;
  void clear() {pairs.clear();}
private:
  set<EOPair> pairs;
};

#endif
