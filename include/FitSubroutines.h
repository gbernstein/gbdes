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
#include "Accum.h"

using namespace std;

const int REF_INSTRUMENT=-1;	// Instrument for reference objects (no fitting)
const int PM_INSTRUMENT =-2;    // Reference catalog with 5d proper motion / parallax solutions
const int TAG_INSTRUMENT=-3;	// Exposure number for tag objects (no fitting nor contrib to stats)
const int NO_INSTRUMENT=-4;

const string stellarAffinity="STELLAR";

const double NO_MAG_DATA = -100.;	// Value entered when there is no valid mag or color

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
  typedef astrometry::CoordAlign Align;
  typedef astrometry::MCat MCat;
  static void fillDetection(Detection* d, const Exposure* e,
			    astrometry::SphericalCoords& fieldProjection,
			    img::FTable& table, long irow,
			    string xKey, string yKey,
			    vector<string>&  xyErrKeys,
			    string magKey, string magErrKey,
			    int magKeyElement, int magErrKeyElement,
			    bool xColumnIsDouble, bool yColumnIsDouble,
			    bool errorColumnIsDouble,
			    bool magColumnIsDouble, bool magErrColumnIsDouble,
			    double magshift,
			    const astrometry::PixelMap* startWcs,
			    bool isTag);
  static void setColor(Detection* d, double color) {
    d->color = color;
  }
  static double getColor(const Detection* d) {
    return d->color;
  }

  // saving results wants an array giving the field projection to use for each catalog.
  static void saveResults(const astrometry::MCat& matches,
			  string outCatalog,
			  string starCatalog,
			  vector<astrometry::SphericalCoords*> catalogProjections);


  static void reportStatistics(const list<typename Astro::Match*>& matches,
			       const vector<Exposure*>& exposures,
			       const vector<typename Astro::Extension*>& extensions,
			       ostream& os);

  // Specialized function for reading detection from a catalog
  // holding proper motion solution.
  // It will incorporate the basic info in (non-PM) Detection d
  // into a new PMDetection object and return it.
  static astrometry::PMDetection*
  makePMDetection(astrometry::Detection* d, const Exposure* e,
		img::FTable& table, long irow,
		string xKey, string yKey,
		string pmRaKey, string pmDecKey, string parallaxKey, string pmCovKey,
		bool xColumnIsDouble, bool yColumnIsDouble,
		bool errorColumnIsDouble,
		const astrometry::PixelMap* startWcs);

  static void handlePMDetection(astrometry::PMDetection* pmd, Detection* d);
  
  static Match*
  makeNewMatch(Detection* d, bool usePM);

  static const int isAstro = 1;
};
struct Photo {
  typedef photometry::Detection Detection;
  typedef photometry::Match Match;
  typedef list<Match*> MCat;
  typedef photometry::SubMap SubMap;
  typedef ExtensionBase<SubMap, Detection> Extension;
  typedef ColorExtensionBase<Match> ColorExtension;
  typedef photometry::PhotoMapCollection Collection;
  typedef photometry::PhotoAlign Align;
  static void fillDetection(Detection* d, const Exposure* e,
			    astrometry::SphericalCoords& fieldProjection,
			    img::FTable& table, long irow,
			    string xKey, string yKey,
			    vector<string>&  xyErrKeys,
			    string magKey, string magErrKey,
			    int magKeyElement, int magErrKeyElement,
			    bool xColumnIsDouble, bool yColumnIsDouble,
			    bool errorColumnIsDouble,
			    bool magColumnIsDouble, bool magErrColumnIsDouble,
			    double magshift,
			    const astrometry::PixelMap* startWcs,
			    bool isTag);
  static void setColor(Detection* d, double color) {
    d->args.color = color;
  }
  static double getColor(const Detection* d) {
    return d->args.color;
  }
  static void saveResults(const list<Match*>& matches,
			  string outCatalog);

  static void reportStatistics(const list<typename Photo::Match*>& matches,
			       const vector<Exposure*>& exposures,
			       const vector<typename Photo::Extension*>& extensions,
			       ostream& os);

  // This is a no-op for photometry:
  static astrometry::PMDetection*
  makePMDetection(photometry::Detection* d, const Exposure* e,
		  img::FTable& table, long irow,
		  string xKey, string yKey,
		  string pmRaKey, string pmDecKey, string parallaxKey, string pmCovKey,
		  bool xColumnIsDouble, bool yColumnIsDouble,
		  bool errorColumnIsDouble,
		  const astrometry::PixelMap* startWcs) {return nullptr;}

  static void handlePMDetection(astrometry::PMDetection* pmd,
				Detection* d) {};   // This is a no-op for photo

  static Match*
  makeNewMatch(Detection* d, bool usePM) {
    return new Match(d);
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

// Read the Fields table from input, copy to output, extract needed info
// If a usable default is given, unassigned field epochs are set to default.
void
readFields(string inputTables,
	   string outCatalog,
	   NameIndex& fieldNames,
	   vector<astrometry::SphericalCoords*>& fieldProjections,
	   vector<double>& fieldEpochs,
	   double defaultEpoch=0.);

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
	      const vector<double>& fieldEpochs,
	      vector<int>& exposureColorPriorities,
	      const list<string>&  useColorList,
	      string inputTables,
	      string outCatalog,
	      const list<string>& skipExposureList,
	      bool useReferenceExposures,
	      bool& outputCatalogAlreadyOpen);

// Read extensions from the table.
// colorExtensions will get filled with ColorExtension objects for color data
// inputYAML is set to produce YAML for all extensions being fit.
// Systematic error for Detections in this Extension will be assigned from
// the column in the extensionTable named sysErrorColumn, if it is non-Null,
// otherwise taken from the value of sysError or referenceSysError as appropriate
template <class S>
vector<typename S::Extension*>
readExtensions(img::FTable& extensionTable,
	       const vector<Instrument*>& instruments,
	       const vector<Exposure*>& exposures,
	       const vector<int>& exposureColorPriorities,
	       vector<typename S::ColorExtension*>& colorExtensions,
	       astrometry::YAMLCollector& inputYAML,
	       bool logging=true);

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

// Add every extension's map to the YAMLCollector and then emit the
// YAML and read into the MapCollection.
// The names of all maps are already in the extension list.
// Returns time spent in critical regions of the addMap() routine.
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
	    typename S::MCat& matches,
	    vector<typename S::Extension*>& extensions,
	    vector<typename S::ColorExtension*>& colorExtensions,
	    const ExtensionObjectSet& skipSet,
	    int minMatches,
	    bool usePM=false);        // If true, create PMMatches

// Read each Extension's objects' data from it FITS catalog
// and place into Detection structures.
template <class S>
void readObjects(const img::FTable& extensionTable,
		 const vector<Exposure*>& exposures,
		 vector<typename S::Extension*>& extensions,
		 vector<astrometry::SphericalCoords*> fieldProjections,
		 bool logging=true); //Give progress updates?

// Read color information from files marked as holding such, insert into
// relevant Matches.
template <class S>
void
readColors(img::FTable extensionTable,
	   vector<typename S::ColorExtension*> colorExtensions,
	   bool logging=true); // Progress reports

// Find all matched Detections that exceed allowable error, then
// delete them from their Match and delete the Detection.
// Does not apply to reference/tag detections.
template <class S>
void
purgeNoisyDetections(double maxError,
		     typename S::MCat& matches,
		     const vector<Exposure*>& exposures,
		     const vector<typename S::Extension*>& extensions);

// Get rid of Matches with too few Detections being fit: delete
// the Match and all of its Detections.
template <class S>
void
purgeSparseMatches(int minMatches,
		   typename S::MCat& matches);

// Get rid of Matches with color outside of specified range.
// Color is always ok if it has NODATA value.
// Note that this even kills Detections that do not need a color for their maps.
template <class S>
void
purgeBadColor(double minColor, double maxColor,
	      typename S::MCat& matches);

template <class S>
void
reserveMatches(typename S::MCat& matches,
	       double reserveFraction,
	       int randomNumberSeed);

template <class S>
map<string, long>
findUnderpopulatedExposures(long minFitExposure,
			    const typename S::MCat& matches,
			    const vector<Exposure*> exposures,
			    const vector<typename S::Extension*> extensions,
			    const typename S::Collection& pmc);

// Fix the parameters of a map, and mark all Detections making
// use of it as clipped so they will not be used in fitting
template <class S>
void
freezeMap(string mapName,
	  const typename S::MCat matches,
	  vector<typename S::Extension*> extensions,
	  typename S::Collection& pmc);

// Report number of unclipped matches and their chisq
template <class S>
void
matchCensus(const typename S::MCat& matches, ostream& os);

// Map and clip reserved matches
template <class S>
void
clipReserved(typename S::Align& ca,
	     double clipThresh,
	     double minimumImprovement,
	     bool clipEntireMatch,
	     bool reportToCerr);

// Function to produce a list of PhotoPriors from an input file
// (for PhotoFit only - has its own source file)
list<photometry::PhotoPrior*>
readPriors(string filename, 
	   const vector<Instrument*>& instruments, 
	   const vector<Exposure*>& exposures, 
	   const vector<Photo::Extension*>& extensions);

#endif
