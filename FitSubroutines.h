// Constants and subroutines shared by the astrometric and photometric matching/fitting codes
// WCSFit and PhotoFit

#ifndef FITSUBROUTINES_H
#define FITSUBROUTINES_H

// Number of threads to use for reading catalogs
#define CATALOG_THREADS 4

#include <list>
#include <set>
#include "StringStuff.h"

// Load the kinds of maps we'll want
#include "PixelMapCollection.h"
#include "PhotoMapCollection.h"
#include "PhotoTemplate.h"
#include "TemplateMap.h"

#include "Pset.h"
#include "FTable.h"
#include "FitsTable.h"
#include "Instrument.h"
#include "YAMLCollector.h"

#include "Match.h"
#include "PhotoMatch.h"

using namespace std;

const int REF_INSTRUMENT=-1;	// Instrument for reference objects (no fitting)
const int TAG_INSTRUMENT=-2;	// Exposure number for tag objects (no fitting nor contrib to stats)
const int NO_INSTRUMENT=-3;

const string stellarAffinity="STELLAR";

const double NO_MAG_DATA = -100.;	// Value entered when there is no valid mag or color
const double NODATA = photometry::PhotoArguments::NODATA;  // Signal for unknown color

const string colorColumnName = "COLOR";
const string colorErrorColumnName = "COLOR_ERR";

// Here is the default character at which to split lists given in parameters strings
const char DefaultListSeperator=',';

// These structures are used so that we can write templated code that will work
// for either astrometric or photometric fitting.
struct Astro {
  typedef astrometry::Detection Detection;
  typedef astrometry::Match Match;
  typedef astrometry::SubMap SubMap;
  typedef ExtensionBase<SubMap, Detection> Extension;
  typedef ColorExtensionBase<Match> ColorExtension;
  typedef astrometry::PixelMapCollection Collection;
  static void fillDetection(Detection* d,
			    img::FTable& table, long irow,
			    double weight,
			    string xKey, string yKey, string errKey,
			    string magKey, string magErrKey,
			    int magKeyElement, int magErrKeyElement,
			    bool errorColumnIsDouble,
			    bool magColumnIsDouble, bool magErrColumnIsDouble,
			    double magshift,
			    const astrometry::PixelMap* startWcs,
			    double sysErrorSq,
			    bool isTag);
  static void setColor(Detection* d, double color) {
    //** ?? Not implemented yet: d->color = color;
  }
  static double getColor(const Detection* d) {
    //** ?? Not implemented yet: return d->color;
    return 0.;
  }
  static const int isAstro = 1;
};
struct Photo {
  typedef photometry::Detection Detection;
  typedef photometry::Match Match;
  typedef photometry::SubMap SubMap;
  typedef ExtensionBase<SubMap, Detection> Extension;
  typedef ColorExtensionBase<Match> ColorExtension;
  typedef photometry::PhotoMapCollection Collection;
  static void fillDetection(Detection* d,
			    img::FTable& table, long irow,
			    double weight,
			    string xKey, string yKey, string errKey,
			    string magKey, string magErrKey,
			    int magKeyElement, int magErrKeyElement,
			    bool errorColumnIsDouble,
			    bool magColumnIsDouble, bool magErrColumnIsDouble,
			    double magshift,
			    const astrometry::PixelMap* startWcs,
			    double sysErrorSq,
			    bool isTag);
  static void setColor(Detection* d, double color) {
    d->args.color = color;
  }
  static double getColor(const Detection* d) {
    return d->args.color;
  }
  static const int isAstro = 0;
};

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

// These functions fill static tables with all the types of Photo/PixelMaps that we
// may be deserializing:
void
loadPixelMapParser(); 
void
loadPhotoMapParser(); 

// Routines to handle table entries that may be float or double, and might be array elements
// First one: parse a key to see if it has form COLUMN[#].  The integer # is returned and the
// key is return as COLUMN.  If no [#], key is unchanged and elementNumber=-1 is returned.
int elementNumber(string& key);
// Next one: see if the table column is double (true) else assume float (false)
bool isDouble(img::FTable f, string key, int elementNumber);
// Retrieve a double-valued number from either float or double column, element of array or scalar cell
double getTableDouble(img::FTable f, string key, int elementNumber, bool isDouble, long irow);

// Class to produce ordering for vector of points to objects
// that have a public "name" member.
template <class T>
class NameSorter {
public:
  NameSorter(const vector<T*>& exposures): ve(exposures) {}
  bool operator()(int i1, int i2) const {
    LessNoCase lnc;
    return lnc(ve[i1]->name, ve[i2]->name);
  }
private:
  const vector<T*>& ve;
};

// This function is used to find degeneracies between exposures and device maps.
// Start with list of free & fixed devices as initial degen/ok, same for exposures.
// Will consider as "ok" any device used in an "ok" exposure and vice-versa.
// The last argument says which exposure/device pairs are used together.
void
findDegeneracies(set<int>& degenerateDevices,
		 set<int>& okDevices,
		 set<int>& degenerateExposures,
		 set<int>& okExposures,
		 const vector<set<int>>& exposuresUsingDevice);

// Figure out which extensions of the FITS file inputTables
// are Instrument or MatchCatalog extensions.
void
inventoryFitsTables(string inputTables,
		    vector<int>& instrumentHDUs,
		    vector<int>& catalogHDUs);

// Read in all the instrument extensions and their device info from input
// FITS file, save useful ones and write to output FITS file.
// The useInstrumentList entries are regexes, empty means use all.
// The final bool argument is set true if we have already created
// the outCatalog FITS file.
vector<Instrument*> readInstruments(vector<int>& instrumentHDUs,
				    list<string>& useInstrumentList,
				    string inputTables,
				    string outCatalog,
				    bool& outputCatalogAlreadyOpen);

// Read the Exposure table into an array.
// Uses the instruments table.
// Fills the exposureColorPriorities table using the useColorPriorities
// list of regexes to establish priority order.
// useReference exposures is set if we want to use exposures from REFERENCE instrument.
// last Boolean is as above.
vector<Exposure*>
readExposures(const vector<Instrument*>& instruments,
	      vector<int>& exposureColorPriorities,
	      const list<string>&  useColorList,
	      string inputTables,
	      string outCatalog,
	      bool useReferenceExposures,
	      bool& outputCatalogAlreadyOpen);

// Read extensions from the table.
// colorExtensions will get filled with ColorExtension objects for color data
// inputYAML is set to produce YAML for all extensions being fit.
template <class S>
vector<typename S::Extension*>
readExtensions(img::FTable& extensionTable,
	       const vector<Instrument*>& instruments,
	       const vector<Exposure*>& exposures,
	       const vector<int>& exposureColorPriorities,
	       vector<typename S::ColorExtension*>& colorExtensions,
	       astrometry::YAMLCollector& inputYAML);

// fix all maps in a photo/pixelMapCollection whose names match
// any regex in the fixMapList.  Also any instrument whose name
// matches gets all device maps fixed.
template <class S>
void
fixMapComponents(typename S::Collection& pmc,
		 const list<string>& fixMapList,
		 const vector<Instrument*>& instruments);

template <class S>
int
findCanonical(Instrument& instr,
	      int iInst,
	      vector<Exposure*>& exposures,
	      vector<typename S::Extension*>& extensions,
	      typename S::Collection& pmc);

template <class S>
void
createMapCollection(const vector<Instrument*>& instruments,
		    const vector<Exposure*>& exposures,
		    const vector<typename S::Extension*> extensions,
		    astrometry::YAMLCollector& inputYAML,
		    typename S::Collection& pmc);

// Set the bool in each Extension indicating whether its map uses color
template <class S>
void
whoNeedsColor(vector<typename S::Extension*> extensions);
	      
// Read a MatchCatalog Extension, recording in each extension the
// objects from it that need to be read from catalog.
// skipSet can give objects to ignore;
// Matches with less than minMatches useful inputs are discarded too.
// Discards Detection if requires color but doesn't have it.
template <class S>
void
readMatches(img::FTable& table,
	    list<typename S::Match*>& matches,
	    vector<typename S::Extension*>& extensions,
	    vector<typename S::ColorExtension*>& colorExtensions,
	    const ExtensionObjectSet& skipSet,
	    int minMatches);

// Read each Extension's objects' data from it FITS catalog
// and place into Detection structures.
template <class S>
void readObjects(const img::FTable& extensionTable,
		 const vector<Exposure*>& exposures,
		 vector<typename S::Extension*>& extensions,
		 double sysError,
		 double referenceSysError);

// Read color information from files marked as holding such, insert into
// relevant Matches.
template <class S>
void
readColors(img::FTable extensionTable,
	   vector<typename S::ColorExtension*> colorExtensions);

// Find all matched Detections that exceed allowable error, then
// delete them from their Match and delete the Detection.
// Does not apply to reference/tag detections.
template <class S>
void
purgeNoisyDetections(double maxError,
		     list<typename S::Match*>& matches,
		     const vector<Exposure*>& exposures,
		     const vector<typename S::Extension*>& extensions);

// Get rid of Matches with too few Detections being fit: delete
// the Match and all of its Detections.
template <class S>
void
purgeSparseMatches(int minMatches,
		   list<typename S::Match*>& matches);

// Get rid of Matches with color outside of specified range.
// Color is always ok if it has NODATA value.
// Note that this even kills Detections that do not need a color for their maps.
template <class S>
void
purgeBadColor(double minColor, double maxColor,
	      list<typename S::Match*>& matches);

template <class S>
void
reserveMatches(list<typename S::Match*>& matches,
	       double reserveFraction,
	       int randomNumberSeed);

template <class S>
map<string, long>
findUnderpopulatedExposures(long minFitExposure,
			    const list<typename S::Match*> matches,
			    const vector<Exposure*> exposures,
			    const vector<typename S::Extension*> extensions,
			    const typename S::Collection& pmc);

// Fix the parameters of a map, and mark all Detections making
// use of it as clipped so they will not be used in fitting
template <class S>
void
freezeMap(string mapName,
	  const list<typename S::Match*> matches,
	  const vector<typename S::Extension*> extensions,
	  typename S::Collection& pmc);

// Report number of unclipped matches and their chisq
template <class S>
void
matchCensus(const list<typename S::Match*>& matches);
#endif
