// Constants and subroutines shared by the astrometric and photometric matching/fitting codes
#include "FitSubroutines.h"
#include "StringStuff.h"
#include <list>
#include "Pset.h"
#include "PhotoTemplate.h"
#include "TemplateMap.h"
#include "PixelMapCollection.h"
#include "PhotoMapCollection.h"

// A helper function that strips white space from front/back of a string and replaces
// internal white space with underscores:
void
spaceReplace(string& s) {
  stripWhite(s);
  s = regexReplace("[[:space:]]+","_",s);
}

// Another helper function to split up a string into a list of whitespace-trimmed strings.
// Get rid of any null ones.
list<string> 
splitArgument(string input, const char listSeperator) {
  list<string> out = split(input, listSeperator);
  for (list<string>::iterator i = out.begin();
       i != out.end(); )  {
    stripWhite(*i);
    if (i->empty()) {
      i = out.erase(i);
    } else {
      ++i;
    }
  }
  return out;
}

// Read parameters from files on command line and from std input.
// First nRequiredArgs must be present and are not names of parameter files.
// Print the usage string (to cerr) and the default parameters (to cout) if
// there are not enough parameters.
// Then dumps final parameter values to cout.
void
processParameters(Pset& parameters, 
		  const string& usage,
		  int nRequiredArgs,
		  int argc,
		  char *argv[]) {

  // Read parameters
  if (argc<nRequiredArgs+1) {
    cerr << usage << endl;
    cout << "#--------- Default Parameters: ---------" << endl;
    parameters.setDefault();
    parameters.dump(cout);
    exit(1);
  }

  parameters.setDefault();
  for (int i=nRequiredArgs+1; i<argc; i++) {
    // Open & read all specified input files
    ifstream ifs(argv[i]);
    if (!ifs) {
      cerr << "Can't open parameter file " << argv[i] << endl;
      exit(1);
    }
    try {
      parameters.setStream(ifs);
    } catch (std::runtime_error &m) {
      cerr << "In file " << argv[i] << ":" << endl;
      quit(m,1);
    }
  }
  // and stdin:
  try {
    parameters.setStream(cin);
  } catch (std::runtime_error &m) {
    cerr << "In stdin:" << endl;
    quit(m,1);
  }

  // List parameters in use
  parameters.dump(cout);
}

// Take a string with this format:
// <regex>=<replace> [, <regex>=<replace>] ...
// and parse it into a stringstuff::RegexReplacement object that will execute
// the regular expression replacements on an input string.
RegexReplacements
parseTranslator(string specString, string errorDescription) {
  list<string> ls = split(specString, DefaultListSeperator);
  RegexReplacements translator;
  for (list<string>::const_iterator i = ls.begin();
       i != ls.end();
       ++i) {
    if (i->empty()) continue;
    list<string> ls2 = split(*i, '=');
    if (ls2.size() != 2) {
      cerr << errorDescription << " has bad translation spec: " << *i << endl;
      exit(1);
    }
    string regex = ls2.front();
    string replacement = ls2.back();
    stripWhite(regex);
    stripWhite(replacement);
    translator.addRegex(regex, replacement);
  }
  return translator;
}

ExtensionObjectSet::ExtensionObjectSet(string filename) {
  if (filename.empty()) return;
  ifstream ifs(filename.c_str());
  if (!ifs) {
    cerr << "Could not open Extension/Object pair file " << filename << endl;
    exit(1);
  }
  string buffer;
  while (stringstuff::getlineNoComment(ifs, buffer)) {
    int extensionNumber;
    long objectNumber;
    istringstream iss(buffer);
    if (!(iss >> extensionNumber >> objectNumber)) {
      cerr << "Bad Extension/Object pair in " << filename
	   << ": <" << buffer << ">" << endl;
      exit(1);
    }
    pairs.insert(EOPair(extensionNumber, objectNumber));
  }
}

bool
ExtensionObjectSet::operator()(int extensionNumber, long objectNumber) const {
  return pairs.count( EOPair(extensionNumber, objectNumber)) >0;
}

void
loadPixelMapParser() {
  // put all new kinds of PixelMap atoms here:
  astrometry::PixelMapCollection::registerMapType<astrometry::TemplateMap1d>();
}

void
  loadPhotoMapParser() {
  // put all new kinds of PhotoMap atoms here:
  photometry::PhotoMapCollection::registerMapType<photometry::PhotoTemplate1d>();
  photometry::PhotoMapCollection::registerMapType<photometry::PhotoRings>();
}

int 
elementNumber(string& key) {
  int out = -1;
  list<string> l = stringstuff::split( stringstuff::regexReplace("^(.*)\\[([0-9]*)\\]$",
								 "\\1;\\2",
								 key),
				       ';');
  if (l.size()==2) {
    key = l.front();
    stripWhite(key);
    out = atoi(l.back().c_str());
  }
  return out;
}

bool 
isDouble(img::FTable f, string key, int elementNumber) {
  bool out = true;
  try {
    if (elementNumber < 0) {
      double d;
      f.readCell(d, key, 0);
    } else {
      vector<double> d;
      f.readCell(d, key, 0);
    }
  } catch (img::FTableError& e) {
	out = false;
  }
  return out;
}

double 
getTableDouble(img::FTable f, string key, int elementNumber, bool isDouble, long irow) {
  double out;
  if (isDouble) {
    if (elementNumber < 0) {
      f.readCell(out, key, irow);
    } else {
      vector<double> v;
      f.readCell(v, key, irow);
      if (elementNumber >= v.size()) {
	cerr << "Requested element " << elementNumber
	     << " of array at key " << key
	     << " but size of array is " << v.size();
	exit(1);
      }
      out = v[elementNumber];
    }
  } else {
    if (elementNumber < 0) {
      float fl;
      f.readCell(fl, key, irow);
      out = fl;
    } else {
      vector<float> v;
      f.readCell(v, key, irow);
      if (elementNumber >= v.size()) {
	cerr << "Requested element " << elementNumber
	     << " of array at key " << key
	     << " but size of array is " << v.size();
	exit(1);
      }
      out = v[elementNumber];
    }
  }
  return out;
}