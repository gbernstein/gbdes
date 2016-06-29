// Constants and subroutines shared by the astrometric and photometric matching/fitting codes
#include "FitSubroutines.h"
#include "StringStuff.h"
#include "Bounds.h"
#include <list>
#include "Pset.h"
#include "PhotoTemplate.h"
#include "PhotoPiecewise.h"
#include "TemplateMap.h"
#include "PiecewiseMap.h"
#include "PixelMapCollection.h"
#include "PhotoMapCollection.h"
#include "FitsTable.h"
#include "Match.h"
#include "PhotoMatch.h"
#include "Random.h"

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

  int positionalArguments;
  try {
    // First read the command-line arguments so we know how many positional
    // arguments precede them.
    positionalArguments = parameters.setFromArguments(argc, argv);
  } catch (std::runtime_error &m) {
    // An error here might indicate someone entered "-help" or something
    cerr << usage << endl;
    cout << "#---- Parameter defaults: ----" << endl;
    parameters.dump(cerr);
    quit(m,1);
  }
  parameters.setDefault();
  if (positionalArguments < nRequiredArgs+1) {
    cerr << usage << endl;
    cout << "#---- Parameter defaults: ----" << endl;
    parameters.dump(cerr);
    exit(1);
  }

  for (int i=nRequiredArgs+1; i<positionalArguments; i++) {
    // Open & read all specified input files
    ifstream ifs(argv[i]);
    if (!ifs) {
      cerr << "Can't open parameter file " << argv[i] << endl;
      cerr << usage << endl;
      exit(1);
    }
    try {
      parameters.setStream(ifs);
    } catch (std::runtime_error &m) {
      cerr << "In file " << argv[i] << ":" << endl;
      quit(m,1);
    }
  }
  // And now re-read the command-line arguments so they take precedence
  parameters.setFromArguments(argc, argv);

  cout << "#" << stringstuff::taggedCommandLine(argc,argv) << endl;
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
  astrometry::PixelMapCollection::registerMapType<astrometry::TemplateMap>();
  astrometry::PixelMapCollection::registerMapType<astrometry::PiecewiseMap>();
}

void
  loadPhotoMapParser() {
  // put all new kinds of PhotoMap atoms here:
  photometry::PhotoMapCollection::registerMapType<photometry::PhotoTemplate>();
  photometry::PhotoMapCollection::registerMapType<photometry::PhotoPiecewise>();
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

// This function is used to find degeneracies between exposures and device maps.
// Start with list of free & fixed devices as initial degen/ok, same for exposures.
// Will consider as "ok" any device used in an "ok" exposure and vice-versa.
// The last argument says which exposure/device pairs are used together.
void
findDegeneracies(set<int>& degenerateDevices,
		 set<int>& okDevices,
		 set<int>& degenerateExposures,
		 set<int>& okExposures,
		 const vector<set<int>>& exposuresUsingDevice) {
  
  // Device/exposures with degeneracy broken this round
  auto nowOkDevices = okDevices;
  auto nowOkExposures = okExposures;
  while (!(nowOkDevices.empty() && nowOkExposures.empty())) {
    // Mark as ok any exposures using an ok device
    for (auto iDev : nowOkDevices) {
      for (auto iExpo: exposuresUsingDevice[iDev]) {
	if (degenerateExposures.erase(iExpo)>0) {
	  // Newly non-degenerate exposure.  Change category
	  nowOkExposures.insert(iExpo);
	  okExposures.insert(iExpo);
	}
      }
    }
    nowOkDevices.clear();

    // Now mark as ok any device that is used by a non-degenerate exposure
    for (auto iDev : degenerateDevices) {
      for (auto iExpo : exposuresUsingDevice[iDev]) {
	if (nowOkExposures.count(iExpo)>0) {
	  // Mark device as not degenerate after all
	  nowOkDevices.insert(iDev);
	  break;
	}
      }
    }
    nowOkExposures.clear();
    for (auto iDev: nowOkDevices) {
      degenerateDevices.erase(iDev);
      okDevices.insert(iDev);
    }
  }
  return;
}

// Figure out which extensions of the FITS file inputTables
// are Instrument or MatchCatalog extensions.
void
inventoryFitsTables(string inputTables,
		    vector<int>& instrumentHDUs,
		    vector<int>& catalogHDUs)  {
  FITS::FitsFile ff(inputTables);
  for (int i=1; i<ff.HDUCount(); i++) {
    FITS::Hdu h(inputTables, FITS::HDUAny, i);
    if (stringstuff::nocaseEqual(h.getName(), "Instrument"))
      instrumentHDUs.push_back(i);
    else if (stringstuff::nocaseEqual(h.getName(), "MatchCatalog"))
      catalogHDUs.push_back(i);
  }
}

// Read in all the instrument extensions and their device info from input
// FITS file, save useful ones and write to output FITS file.
// The useInstrumentList entries are regexes, empty means use all.
// The final bool argument is set true if we have already created
// the outCatalog FITS file.
vector<Instrument*> readInstruments(vector<int>& instrumentHDUs,
				    list<string>& useInstrumentList,
				    string inputTables,
				    string outCatalog,
				    bool& outputCatalogAlreadyOpen) {
  vector<Instrument*> instruments;
  for (int iextn : instrumentHDUs) {
    FITS::FitsTable ft(inputTables, FITS::ReadOnly, iextn);
    Assert( stringstuff::nocaseEqual(ft.getName(), "Instrument"));

    string instrumentName;
    int instrumentNumber;
    if (!ft.header()->getValue("Name", instrumentName)
	|| !ft.header()->getValue("Number", instrumentNumber)) {
      cerr << "Could not read name and/or number of instrument at extension "
	   << iextn << endl;
    }
    spaceReplace(instrumentName);

    if (regexMatchAny(useInstrumentList, instrumentName))  {
      // This is an instrument we will use.  Save to the output file first.
      if (instrumentNumber >= instruments.size())
	instruments.resize(instrumentNumber+1);

      // Copy extension to output file:
      FITS::Flags outFlags = FITS::ReadWrite+FITS::Create;
      if (!outputCatalogAlreadyOpen) {
	outFlags = outFlags + FITS::OverwriteFile;
	outputCatalogAlreadyOpen = true;
      }
      FITS::FitsTable out(outCatalog, outFlags, -1);
      out.setName("Instrument");
      img::FTable ff=ft.extract();
      out.adopt(ff);

      Instrument* instptr = new Instrument(instrumentName);
      instruments[instrumentNumber] = instptr;
      string band;
      if (!ff.header()->getValue("Band",band)) {
	instptr->band = instptr->name;  // Use instrument name for BAND if not found
      } else {
	spaceReplace(band);
	instptr->band = band;
      }
      
      vector<string> devnames;
      vector<double> vxmin;
      vector<double> vxmax;
      vector<double> vymin;
      vector<double> vymax;
      ff.readCells(devnames, "Name");
      ff.readCells(vxmin, "XMin");
      ff.readCells(vxmax, "XMax");
      ff.readCells(vymin, "YMin");
      ff.readCells(vymax, "YMax");
      for (int j=0; j<devnames.size(); j++) {
	spaceReplace(devnames[j]);
	instptr->addDevice( devnames[j],
			    Bounds<double>( vxmin[j], vxmax[j], vymin[j], vymax[j]));
      }
      instruments[instrumentNumber] = instptr;
    } // Done reading an instrument.
  }

  return instruments;
}

// Read the Exposure table into an array.
vector<Exposure*>
readExposures(const vector<Instrument*>& instruments,
	      vector<int>& exposureColorPriorities,
	      const list<string>&  useColorList,
	      string inputTables,
	      string outCatalog,
	      const list<string>& skipExposureList,
	      bool useReferenceExposures,
	      bool& outputCatalogAlreadyOpen) {

  // Read the exposure table
  FITS::FitsTable ft(inputTables, FITS::ReadOnly, "Exposures");
  FITS::Flags outFlags = FITS::ReadWrite+FITS::Create;
  if (!outputCatalogAlreadyOpen) {
    outFlags = outFlags + FITS::OverwriteFile;
    outputCatalogAlreadyOpen = true;
  }
  FITS::FitsTable out(outCatalog, outFlags, "Exposures");
  img::FTable ff = ft.extract();
  out.adopt(ff);
  vector<string> names;
  vector<double> ra;
  vector<double> dec;
  vector<int> fieldNumber;
  vector<int> instrumentNumber;
  vector<double> airmass;
  vector<double> exptime;
  // ??? Add MJD??
  ff.readCells(names, "Name");
  ff.readCells(ra, "RA");
  ff.readCells(dec, "Dec");
  ff.readCells(fieldNumber, "fieldNumber");
  ff.readCells(instrumentNumber, "InstrumentNumber");
  ff.readCells(airmass, "Airmass");
  ff.readCells(exptime, "Exptime");

  // Initialize our output arrays to not-in-use values
  vector<Exposure*> exposures(names.size(), 0);
  exposureColorPriorities = vector<int>(names.size(), -1);
  for (int i=0; i<names.size(); i++) {
    spaceReplace(names[i]);

    if (regexMatchAny(skipExposureList, names[i]))
      continue;  // do not use exposures on the list.
    
    // See if this exposure name matches any of the color exposures
    int priority = 0;
    for (auto j = useColorList.begin();
	 j != useColorList.end();
	 ++j, ++priority) {
      if (stringstuff::regexMatch(*j, names[i])) {
	// Yes: give this exposure a priority according to order of the exposure
	// it first matches.
	exposureColorPriorities[i] = priority;
	break;  // First priority is given to lower numbers.
      }
    }
    // Are we going to want use data this exposure?
    // Yes, if it's a reference exposure and we are using references
    bool useThisExposure = instrumentNumber[i]==REF_INSTRUMENT 
      && useReferenceExposures;
    // or if it's a normal exposure using an instrument that has been included
    useThisExposure = useThisExposure || (instrumentNumber[i]>=0 && 
					  instruments[instrumentNumber[i]]);

    if (useThisExposure) {
      // The projection we will use for this exposure:
      astrometry::Gnomonic gn(astrometry::Orientation(astrometry::SphericalICRS(ra[i]*DEGREE,
										dec[i]*DEGREE)));
      auto expo = new Exposure(names[i],gn);
      expo->field = fieldNumber[i];
      expo->instrument = instrumentNumber[i];
      expo->airmass = airmass[i];
      expo->exptime = exptime[i];
      // ??? MJD
      exposures[i] = expo;
    }
  } // End exposure loop
  return exposures;
}

template <class S>
void
fixMapComponents(typename S::Collection& pmc,
		 const list<string>& fixMapList,
		 const vector<Instrument*>& instruments) {
  set<string> fixTheseMaps;
  for (auto iName : pmc.allMapNames()) {
    if (stringstuff::regexMatchAny(fixMapList, iName))
      fixTheseMaps.insert(iName);
  }
  // Add names of all devices of instruments on the fixMapList
  for (auto instptr : instruments) {
    if (!instptr) continue; // Not in use
    if (stringstuff::regexMatchAny(fixMapList, instptr->name)) {
      // Loop through all devices
      for (int i=0; i<instptr->nDevices; i++) {
	string devMap = instptr->name + "/" + instptr->deviceNames.nameOf(i);
	// Freeze the device's map if it's in use
	if (pmc.mapExists(devMap))
	  fixTheseMaps.insert(devMap);
      }
    }
  }
  pmc.setFixed(fixTheseMaps);
}

template <class S>
vector<typename S::Extension*>
readExtensions(img::FTable& extensionTable,
	       const vector<Instrument*>& instruments,
	       const vector<Exposure*>& exposures,
	       const vector<int>& exposureColorPriorities,
	       vector<typename S::ColorExtension*>& colorExtensions,
	       astrometry::YAMLCollector& inputYAML) {
      
  vector<typename S::Extension*> extensions(extensionTable.nrows(), 0);
  colorExtensions = vector<typename S::ColorExtension*>(extensionTable.nrows(), 0);
  for (int i=0; i<extensionTable.nrows(); i++) {
    int iExposure;
    extensionTable.readCell(iExposure, "Exposure", i);

    if (iExposure < 0 || iExposure >= exposures.size()) {
      cerr << "Extension " << i << " has invalid exposure number " << iExposure << endl;
      exit(1);
    }

    // Determine whether this extension might be used to provide colors
    int colorPriority = exposureColorPriorities[iExposure];
    if (colorPriority >= 0) {
      auto* ce = new typename S::ColorExtension;
      ce->priority = colorPriority;
      colorExtensions[i] = ce;
    }
	
    if (!exposures[iExposure]) continue;
    
    typename S::Extension* extn = new typename S::Extension;
    extn->exposure = iExposure;
    const Exposure& expo = *exposures[iExposure];

    int iDevice;
    extensionTable.readCell(iDevice, "Device", i);
    extn->device = iDevice;
    extn->airmass = expo.airmass;
    extn->magshift = +2.5*log10(expo.exptime);
      
    // Create the starting WCS for the exposure
    string s;
    extensionTable.readCell(s, "WCSIn", i);
    if (stringstuff::nocaseEqual(s, "_ICRS")) {
      // Create a Wcs that just takes input as RA and Dec in degrees;
      astrometry::IdentityMap identity;
      astrometry::SphericalICRS icrs;
      extn->startWcs = new astrometry::Wcs(&identity, icrs, "ICRS_degrees", DEGREE);
    } else {
      istringstream iss(s);
      astrometry::PixelMapCollection pmcTemp;
      if (!pmcTemp.read(iss)) {
	cerr << "Could not deserialize starting WCS for extension #" << i << endl;
	exit(1);
      }
      string wcsName = pmcTemp.allWcsNames().front();
      extn->startWcs = pmcTemp.cloneWcs(wcsName);
    }

    // destination projection for startWCS is the Exposure projection,
    // so that any exposure-level magnitude corrections use this coord system
    extn->startWcs->reprojectTo(*expo.projection);


    // Extract the map specifications for this extension from the input
    // YAML files.
    if ( expo.instrument < 0) {
      // Tag & reference exposures have no instruments and no fitting
      // being done.  Coordinates are fixed to xpix = xw.
      // WCS will be the same for everything in this field
      // ???      extn->wcsName = fieldNames.nameOf(ifield);
      extn->mapName = astrometry::IdentityMap().getName();
    } else {
      // Real instrument, make a map combining its exposure with its Device map:
      astrometry::YAMLCollector::Dictionary d;
      d["INSTRUMENT"] = instruments[expo.instrument]->name;
      d["DEVICE"] = instruments[expo.instrument]->deviceNames.nameOf(extn->device);
      d["EXPOSURE"] = expo.name;
      d["BAND"] = instruments[expo.instrument]->band;
      // This will be the name of the extension's final WCS and photometry map:
      extn->wcsName = d["EXPOSURE"] + "/" + d["DEVICE"];
      // And this is the name of the PixelMap underlying the WCS, or
      // the name of the photometric map (for which we do not want the "base" part).
      if (S::isAstro)
	extn->mapName = extn->wcsName + "/base";
      else
	extn->mapName = extn->wcsName;
	
	  
      if (!inputYAML.addMap(extn->mapName,d)) {
	cerr << "Input YAML files do not have complete information for map "
	     << extn->mapName
	     << endl;
	exit(1);
      }
    }
    extensions[i] = extn;
  }  // End extension loop
  return extensions;
}


// Routine to find an exposure that should have exposure map set to identity
// to resolve exposure/device degeneracies in the model.  Returns <0 if
// no canonical exposure is necessary for this instrument.  Takes
// as input references to an Instrument (and its id number), the map collection
// being studied, and the exposure/extension tables.

template <class S>
int
findCanonical(Instrument& instr,
	      int iInst,
	      vector<Exposure*>& exposures,
	      vector<typename S::Extension*>& extensions,
	      typename S::Collection& pmc)
{
  // Classify the device maps for this instrument
  set<int> fixedDevices;
  set<int> freeDevices;
  set<int> unusedDevices;

  // And the exposure maps as well:
  set<int> fixedExposures;
  set<int> freeExposures;
  set<int> unusedExposures;
  set<int> itsExposures;  // All exposure numbers using this instrument
      
  for (int iDev=0; iDev<instr.nDevices; iDev++) {
    string mapName = instr.name + "/" + instr.deviceNames.nameOf(iDev);
    if (!pmc.mapExists(mapName)) {
      unusedDevices.insert(iDev);
      continue;
    }
    instr.mapNames[iDev] = mapName;
    if (pmc.getFixed(mapName)) {
      fixedDevices.insert(iDev);
    } else {
      freeDevices.insert(iDev);
    }
  }

  for (int iExpo=0; iExpo < exposures.size(); iExpo++) {
    if (!exposures[iExpo]) continue; // Skip unused
    if (exposures[iExpo]->instrument==iInst) {
      itsExposures.insert(iExpo);
      string mapName = exposures[iExpo]->name;
      if (!pmc.mapExists(mapName)) {
	unusedExposures.insert(iExpo);
      } else if (pmc.getFixed(mapName)) {
	fixedExposures.insert(iExpo);
      } else {
	freeExposures.insert(iExpo);
      }
    }
  }

  //**/cerr << "Done collecting fixed/free" << endl;

  // Now take an inventory of all extensions to see which device
  // solutions are used in coordination with which exposure solutions
  vector<set<int>> exposuresUsingDevice(instr.nDevices);
  for (auto extnptr : extensions) {
    if (!extnptr) continue; // Extension not in use.
    int iExpo = extnptr->exposure;
    if (itsExposures.count(iExpo)>0) {
      // Extension is from one of the instrument's exposures
      int iDev = extnptr->device;
      if (pmc.dependsOn(extnptr->mapName, exposures[iExpo]->name)
	  && pmc.dependsOn(extnptr->mapName, instr.mapNames[iDev])) {
	// If this extension's map uses both the exposure and device
	// maps, then enter this dependence into our sets
	exposuresUsingDevice[iDev].insert(iExpo);
	if (unusedExposures.count(iExpo)>0 ||
	    unusedDevices.count(iDev)>0) {
	  cerr << "ERROR: Logic problem: extension map "
	       << extnptr->mapName
	       << " is using allegedly unused exposure or device map"
	       << endl;
	  exit(1);
	}
      }
    }
  }

  //**/cerr << "Done building exposure/device graph" << endl;

  // We have a degeneracy if there is a set of exposures and devices such
  // that
  // * the device maps and exposure maps are free
  // * and the device maps are used only in exposures from this set
  // * and the exposures contain only devices from this set.
  // See here if such a subset exists, in which case we need
  // to fix one of the exposure maps as a "canonical" exposure with
  // Identity map.

  // Split in-use devices into those potentially degenerate and those
  // either fixed or tied to a fixed solution
  auto degenerateDevices = freeDevices;
  auto okDevices = fixedDevices;
  auto degenerateExposures = freeExposures;
  auto okExposures = fixedExposures;

  // propagate "ok-ness" from devices to exposures and back until no more.
  findDegeneracies(degenerateDevices,
		   okDevices,
		   degenerateExposures,
		   okExposures,
		   exposuresUsingDevice);

  // If there are no degenerate device maps, we are done!
  // Return a value indicating no canonical needed at all.
  if (degenerateDevices.empty())
    return -1;

  if (degenerateExposures.empty()) {
    cerr << "Logic problem: Instrument " << instr.name
	 << " came up with degenerate devices but not exposures:"
	 << endl;
    exit(1);
  }

  int canonicalExposure = -1;
  // Find the exposure using the most degenerate devices
  {
    int maxDevices = 0;
    for (auto iExpo : degenerateExposures) {
      int nDevices = 0;
      for (auto iDev : degenerateDevices) {
	if (exposuresUsingDevice[iDev].count(iExpo)>0)
	  nDevices++;
      }
      if (nDevices > maxDevices) {
	// Save this exposure as best to use, no need
	// to continue if it uses all devices.
	maxDevices = nDevices;
	canonicalExposure = iExpo;
	if (maxDevices==degenerateDevices.size())
	  break;
      }
    }
  }

  if (canonicalExposure < 0) {
    cerr << "Failed to locate a canonical exposure for " << instr.name
	 << endl;
    exit(1);
  }

  // Check that fixing this exposure map will resolve degeneracies.
  degenerateExposures.erase(canonicalExposure);
  okExposures.insert(canonicalExposure);
  findDegeneracies(degenerateDevices,
		   okDevices,
		   degenerateExposures,
		   okExposures,
		   exposuresUsingDevice);
  if (!degenerateDevices.empty()) {
    cerr << "But canonical did not resolve exposure/device degeneracy."
	 << endl;
    exit(1);
  }
  return canonicalExposure;
}

// Add every extension's map to the YAMLCollector and then emit the
// YAML and read into the MapCollection.
// The names of all maps are already in the extension list.
template <class S>
void
createMapCollection(const vector<Instrument*>& instruments,
		    const vector<Exposure*>& exposures,
		    const vector<typename S::Extension*> extensions,
		    astrometry::YAMLCollector& inputYAML,
		    typename S::Collection& pmc) {
  for (auto extnptr : extensions) {
    if (!extnptr) continue;  // Not in use.
    auto& expo = *exposures[extnptr->exposure];
    // Extract the map specifications for this extension from the input
    // YAML files.
    astrometry::YAMLCollector::Dictionary d;
    if ( expo.instrument >= 0) {
      // Real instrument, need the translation dictionary
      d["INSTRUMENT"] = instruments[expo.instrument]->name;
      d["DEVICE"] = instruments[expo.instrument]->deviceNames.nameOf(extnptr->device);
      d["EXPOSURE"] = expo.name;
      d["BAND"] = instruments[expo.instrument]->band;
      }
      if (!inputYAML.addMap(extnptr->mapName,d)) {
	cerr << "Input YAML files do not have complete information for map "
	     << extnptr->mapName
	     << endl;
	exit(1);
      }
  }  // End extension loop
      
  {
    istringstream iss(inputYAML.dump());
    if (!pmc.read(iss)) {
      cerr << "Failure parsing the final YAML map specs" << endl;
      exit(1);
    }
  }
}

template <class S>
void
whoNeedsColor(vector<typename S::Extension*> extensions) {
  for (auto extnptr : extensions) {
    if (extnptr) {
      extnptr->needsColor = extnptr->map->needsColor();
    }
  }
}
	      
// Read a MatchCatalog Extension, recording in each extension the
// objects from it that need to be read from catalog.
template <class S>
void
readMatches(img::FTable& table,
	    list<typename S::Match*>& matches,
	    vector<typename S::Extension*>& extensions,
	    vector<typename S::ColorExtension*>& colorExtensions,
	    const ExtensionObjectSet& skipSet,
	    int minMatches) {

  vector<int> seq;
  vector<LONGLONG> extn;
  vector<LONGLONG> obj;
  table.readCells(seq, "SequenceNumber");
  table.readCells(extn, "Extension");
  table.readCells(obj, "Object");

  // Smaller collections for each match
  vector<long> matchExtns;
  vector<long> matchObjs;
  // These variables determine what the highest-priority color information available
  // in the match is so far.
  long matchColorExtension = -1;
  long matchColorObject = 0;
  int colorPriority = -1;

  for (int i=0; i<=seq.size(); i++) {
    if (i>=seq.size() || seq[i]==0) {
      // Processing the previous (or final) match.
      int nValid = matchExtns.size();
      if (matchColorExtension < 0) {
	// There is no color information for this match. 
	// Discard any detection which requires a color to produce its map:
	nValid = 0;
	for (int j=0; j<matchExtns.size(); j++) {
	  Assert(extensions[matchExtns[j]]);
	  if (extensions[matchExtns[j]]->needsColor) {
	    // Mark this detection as useless
	    matchExtns[j] = -1;
	  } else {
	    nValid++;
	  }
	}
      }

      if (nValid >= minMatches) {
	// Make a match from the valid entries, and note need to get data for the detections and color
	typename S::Match* m=0;
	for (int j=0; j<matchExtns.size(); j++) {
	  if (matchExtns[j]<0) continue;  // Skip detections starved of their color
	  auto d = new typename S::Detection;
	  d->catalogNumber = matchExtns[j];
	  d->objectNumber = matchObjs[j];
	  extensions[matchExtns[j]]->keepers.insert(std::pair<long,typename S::Detection*>
						    (matchObjs[j], d));
	  if (m)
	    m->add(d);
	  else
	    m = new typename S::Match(d);
	}
	matches.push_back(m);

	if (matchColorExtension >=0) {
	  // Tell the color catalog that it needs to look this guy up:
	  Assert(colorExtensions[matchColorExtension]);
	  colorExtensions[matchColorExtension]->keepers.insert
	    (std::pair<long,typename S::Match*>(matchColorObject,
						matches.back()));
	}
      }
      // Clear out previous Match:
      matchExtns.clear();
      matchObjs.clear();
      matchColorExtension = -1;
      colorPriority = -1;
    } // Finished processing previous match

    // If we done reading entries, quit this loop
    if (i >= seq.size()) break;

    // continue loop if this is an object to ignore
    if ( skipSet(extn[i],obj[i]) ) continue;

    // Note extn/obj number of detections with useful mag data
    if (extensions[extn[i]]) {
      // Record a Detection in a useful extension:
      matchExtns.push_back(extn[i]);
      matchObjs.push_back(obj[i]);
    }

    // Record if we've got color information here
    if (colorExtensions[extn[i]]) {
      int newPriority = colorExtensions[extn[i]]->priority;
      if (newPriority >= 0 && (colorPriority < 0 || newPriority < colorPriority)) {
	// This detection holds color info that we want
	colorPriority = newPriority;
	matchColorExtension = extn[i];
	matchColorObject = obj[i];
      }
    }
  } // End loop of catalog entries
}

// Subroutine to get what we want from a catalog entry for WCS fitting
inline
void
Astro::fillDetection(Astro::Detection* d,
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
		     bool isTag) {
  table.readCell(d->xpix, xKey, irow);
  table.readCell(d->ypix, yKey, irow);

  double sigma = isTag ? 0. : getTableDouble(table, errKey, -1, magColumnIsDouble,irow);
  sigma = std::sqrt(sysErrorSq + sigma*sigma);
  d->sigma = sigma;

  startWcs->toWorld(d->xpix, d->ypix, d->xw, d->yw);
  Matrix22 dwdp = startWcs->dWorlddPix(d->xpix, d->ypix);
  // no clips on tags
  double wt = isTag ? 0. : pow(sigma,-2.);
  d->clipsqx = wt / (dwdp(0,0)*dwdp(0,0)+dwdp(0,1)*dwdp(0,1));
  d->clipsqy = wt / (dwdp(1,0)*dwdp(1,0)+dwdp(1,1)*dwdp(1,1));
  d->wtx = d->clipsqx * weight;
  d->wty = d->clipsqy * weight;
}

// And a routine to get photometric information too
inline
void
Photo::fillDetection(Photo::Detection* d,
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
		     bool isTag) {
  table.readCell(d->args.xDevice, xKey, irow);
  table.readCell(d->args.yDevice, yKey, irow);
  startWcs->toWorld(d->args.xDevice, d->args.yDevice,
		    d->args.xExposure, d->args.yExposure);

  // Get the mag input and its error
  d->magIn = getTableDouble(table, magKey, magKeyElement, magColumnIsDouble,irow)
    + magshift;
  double sigma = getTableDouble(table, magErrKey, magErrKeyElement, magErrColumnIsDouble,irow);

  // Map to output and estimate output error
  d->magOut = d->map->forward(d->magIn, d->args);

  sigma = std::sqrt(sysErrorSq + sigma*sigma);
  d->sigma = sigma;
  sigma *= d->map->derivative(d->magIn, d->args);

  // no clips on tags
  double wt = isTag ? 0. : pow(sigma,-2.);
  d->clipsq = wt;
  d->wt = wt * weight;
}

inline
void
Astro::matchMean(const Astro::Match& m,
		 double& mx, double& my, double& mmag,
		 double& wtot) {
  double wx,wy;
  m.centroid(mx,my,wx,wy);
  wtot = 0.5*(wx+wy); // Use mean weight
}
inline
void
Photo::matchMean(const Photo::Match& m,
		 double& mx, double& my, double& mmag,
		 double& wtot) {
  m.getMean(mmag,wtot);
}

inline
void
Astro::getOutputA(const Astro::Detection& d,
		  double mx, double my, double wtot,
		  double& xp, double& yp, double& sigma,
		  double& xerrpix, double& yerrpix,
		  double& xw, double& yw, double& sigmaw,
		  double& xerrw, double& yerrw,
		  double& wf) {
  wf = d.isClipped ? 0. : 0.5 * (d.wtx + d.wty) / wtot;  
  sigmaw = 0.5*(d.clipsqx + d.clipsqy);
  sigmaw = (sigmaw > 0.) ? 1./sqrt(sigmaw) : 0.;
  xp = d.xpix;
  yp = d.ypix;
  xw = d.xw;
  yw = d.yw;
  sigma = 0.5*d.sigma;
  
  // Calculate residuals if we have a centroid for the match:
  if (mx!=0. || my!=0.) {
    double xcpix=xp, ycpix=yp;
    Assert (d.map);
    try {
      d.map->toPix(mx, my, xcpix, ycpix);
    } catch (astrometry::AstrometryError& e) {
      cerr << "WARNING: toPix failure in map " << d.map->getName()
	   << " at world coordinates (" << mx << "," << my << ")"
	   << " approx pixel coordinates (" << xp << "," << yp << ")"
	   << endl;
      // Just set inverse map to pixel coords (i.e. zero pixel error reported)
      xcpix = xp;
      ycpix = yp;
    }
    xerrw = xw - mx;
    yerrw = yw - my;
    xerrpix = xp - xcpix;
    yerrpix = yp - ycpix;
  }
  return;
}
inline
void
Photo::getOutputP(const Photo::Detection& d,
		  double mmag, double wtot,
		  double& xDevice, double& yDevice,
		  double& xExposure, double& yExposure,
		  double& mag, double& magerr, double& sigma,
		  double& wf) {
  xDevice = d.args.xDevice;
  yDevice = d.args.yDevice;
  xExposure = d.args.xExposure;
  yExposure = d.args.yExposure;
  mag = d.magOut;
  magerr = mmag==0. ? 0. : mag - mmag;
  sigma = d.clipsq;
  sigma = sigma>0. ? 1./sqrt(sigma) : 0.;
  wf = d.isClipped ? 0. : d.wt / wtot;
}
  

// Read each Extension's objects' data from it FITS catalog
// and place into Detection structures.
template <class S>
void readObjects(const img::FTable& extensionTable,
		 const vector<Exposure*>& exposures,
		 vector<typename S::Extension*>& extensions,
		 double sysError,
		 double referenceSysError) {

  // Should be safe to multithread this loop as different threads write
  // only to distinct parts of memory.  Protect the FITS table read though.
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic,CATALOG_THREADS)
#endif
  for (int iext = 0; iext < extensions.size(); iext++) {
    if (!extensions[iext]) continue; // Skip unused 
    // Relevant structures for this extension
    typename S::Extension& extn = *extensions[iext];
    Exposure& expo = *exposures[extn.exposure];
    if (extn.keepers.empty()) continue; // or useless

    string filename;
    extensionTable.readCell(filename, "Filename", iext);
    /**/if (iext%50==0) cerr << "# Reading object catalog " << iext
			     << "/" << extensions.size()
			     << " from " << filename 
			     << " seeking " << extn.keepers.size()
			     << " objects" << endl;
    int hduNumber;
    string xKey;
    string yKey;
    string idKey;
    string magKey;
    string magErrKey;
    double weight;
    string errKey;
    extensionTable.readCell(hduNumber, "Extension", iext);
    extensionTable.readCell(xKey, "xKey", iext);
    extensionTable.readCell(yKey, "yKey", iext);
    extensionTable.readCell(idKey, "idKey", iext);

    if (S::isAstro) {
      extensionTable.readCell(errKey, "errKey", iext);
      extensionTable.readCell(weight, "Weight", iext);
    } else {
      extensionTable.readCell(magKey, "magKey", iext);
      extensionTable.readCell(magErrKey, "magErrKey", iext);
      extensionTable.readCell(weight, "magWeight", iext);
    }
    
    const typename S::SubMap* sm=extn.map;

    bool isReference = (expo.instrument == REF_INSTRUMENT);
    bool isTag = (expo.instrument == TAG_INSTRUMENT);

    double sysErrorSq = pow(isReference ? referenceSysError : sysError, 2.);
    
    astrometry::Wcs* startWcs = extn.startWcs;

    if (!startWcs) {
      cerr << "Failed to find initial Wcs for exposure " << expo.name
	   << endl;
      exit(1);
    }

    bool useRows = stringstuff::nocaseEqual(idKey, "_ROW");
    // What we need to read from the FitsTable:
    vector<string> neededColumns;
    if (!useRows)
      neededColumns.push_back(idKey);
    neededColumns.push_back(xKey);
    neededColumns.push_back(yKey);
      
    int magKeyElement;
    int magErrKeyElement;
    if (S::isAstro) {
      if (!isTag)
	neededColumns.push_back(errKey);
    } else {
      // Be willing to get an element of array-valued bintable cell
      // for mag or magerr.  Syntax would be
      // MAGAPER[4]  to get 4th (0-indexed) element of MAGAPER column
      magKeyElement = elementNumber(magKey);
      magErrKeyElement = elementNumber(magErrKey);
      neededColumns.push_back(magKey);
      neededColumns.push_back(magErrKey);
    }
    img::FTable ff;
#ifdef _OPENMP
#pragma omp critical(fitsio)
#endif
    {
      FITS::FitsTable ft(filename, FITS::ReadOnly, hduNumber);
      ff = ft.extract(0, -1, neededColumns);
    }
    vector<long> id;
    if (useRows) {
      id.resize(ff.nrows());
      for (long i=0; i<id.size(); i++)
	id[i] = i;
    } else {
      ff.readCells(id, idKey);
    }
    Assert(id.size() == ff.nrows());

    bool magColumnIsDouble;
    bool magErrColumnIsDouble;
    bool errorColumnIsDouble;
    double magshift = extn.magshift;

    if (S::isAstro) {
      errorColumnIsDouble = isDouble(ff, errKey, -1);
    } else {
      magColumnIsDouble = isDouble(ff, magKey, magKeyElement);
      magErrColumnIsDouble = isDouble(ff, magErrKey, magErrKeyElement);
    }

    for (long irow = 0; irow < ff.nrows(); irow++) {
      auto pr = extn.keepers.find(id[irow]);
      if (pr == extn.keepers.end()) continue; // Not a desired object

      // Have a desired object now.  Fill its Detection structure
      typename S::Detection* d = pr->second;
      extn.keepers.erase(pr);
      d->map = sm;

      S::fillDetection(d, ff, irow,
		       weight,
		       xKey, yKey, errKey, magKey, magErrKey,
		       magKeyElement, magErrKeyElement,
		       errorColumnIsDouble, magColumnIsDouble, magErrColumnIsDouble,
		       magshift,
		       startWcs, sysErrorSq, isTag);
    } // End loop over catalog objects

    if (!extn.keepers.empty()) {
      cerr << "Did not find all desired objects in catalog " << filename
	   << " extension " << hduNumber
	   << endl;
      exit(1);
    }
  } // end loop over catalogs to read
  
}

// Read color information from files marked as holding such, insert into
// relevant Matches.
template <class S>
void
readColors(img::FTable extensionTable,
	   vector<typename S::ColorExtension*> colorExtensions) {

  for (int iext = 0; iext < colorExtensions.size(); iext++) {
    if (!colorExtensions[iext]) continue; // Skip unused 
    auto& extn = *colorExtensions[iext];
    if (extn.keepers.empty()) continue; // Not using any colors from this catalog
    string filename;
    extensionTable.readCell(filename, "Filename", iext);
    /**/ cerr << "# Reading color catalog " << iext
	      << "/" << colorExtensions.size()
	      << " from " << filename << endl;
    int hduNumber;
    extensionTable.readCell(hduNumber, "Extension", iext);
    string idKey;
    extensionTable.readCell(idKey, "idKey", iext);
    string colorExpression;
    try {
      extensionTable.readCell(colorExpression, "colorExpression", iext);
    } catch (img::FTableError& e) {
      // If there is no colorExpression column, use a default:
      colorExpression = "COLOR";
    }
    stripWhite(colorExpression);
    if (colorExpression.empty()) {
      cerr << "No colorExpression specified for filename " << filename
	   << " HDU " << hduNumber
	   << endl;
      exit(1);
    }

    // Read the entire catalog for this extension
    FITS::FitsTable ft(filename, FITS::ReadOnly, hduNumber);
    img::FTable ff = ft.use();

    bool useRows = stringstuff::nocaseEqual(idKey, "_ROW");
    vector<long> id;
    if (useRows) {
      id.resize(ff.nrows());
      for (long i=0; i<id.size(); i++)
	id[i] = i;
    } else {
      ff.readCells(id, idKey);
    }
    Assert(id.size() == ff.nrows());

    vector<double> color(id.size(), 0.);
    ff.evaluate(color, colorExpression);

    for (long irow = 0; irow < ff.nrows(); irow++) {
      auto pr = extn.keepers.find(id[irow]);
      if (pr == extn.keepers.end()) continue; // Not a desired object

      // Have a desired object. Put the color into everything it matches
      typename S::Match* m = pr->second;
      extn.keepers.erase(pr);

      for ( auto detptr : *m)
	S::setColor(detptr,color[irow]);

    } // End loop over catalog objects

    if (!extn.keepers.empty()) {
      cerr << "Did not find all desired objects in catalog " << filename
	   << " extension " << hduNumber
	   << " " << extn.keepers.size() << " left, first ID is "
	   << extn.keepers.begin()->first
	   << endl;
      exit(1);
    }
  } // end loop over catalogs to read
}

// Find all matched Detections that exceed allowable error, then
// delete them from their Match and delete the Detection.
// Does not apply to reference/tag detections.
template <class S>
void
purgeNoisyDetections(double maxError,
		     list<typename S::Match*>& matches,
		     const vector<Exposure*>& exposures,
		     const vector<typename S::Extension*>& extensions) {
  for (auto mptr : matches) {
    auto j=mptr->begin(); 
    while (j != mptr->end()) {
      auto d = *j; // Yields pointer to each detection
      // Keep it if pixel error is small or if it's from a tag or reference "instrument"
      if ( d->sigma > maxError
	   && exposures[ extensions[d->catalogNumber]->exposure]->instrument >= 0) {
	j = mptr->erase(j, true);  // will delete the detection too.
      } else {
	++j;
      }
    }
  }
}

// Get rid of Matches with too few Detections being fit: delete
// the Match and all of its Detections.
template <class S>
void
purgeSparseMatches(int minMatches,
		   list<typename S::Match*>& matches) {
  auto im = matches.begin();
  while (im != matches.end() ) {
    (*im)->countFit();
    if ( (*im)->fitSize() < minMatches) {
      // Remove entire match if it's too small, and kill its Detections too
      (*im)->clear(true);
      im = matches.erase(im);
    } else {
      ++im;
    }
  }
}

// Get rid of Matches with color outside of specified range.
// Color is always ok if it has NODATA value.
// Note that this even kills Detections that do not need a color for their maps.
template <class S>
void
purgeBadColor(double minColor, double maxColor,
	      list<typename S::Match*>& matches) {
  auto im = matches.begin();
  while (im != matches.end() ) {
    if ((*im)->size() > 0) {
      // See if color is in range, using color from first Detection (they should
      // all have the same color).
      double color = S::getColor(*(*im)->begin());
      if (color != NODATA && (color < minColor || color > maxColor)) {
	// Remove entire match if it's too small, and kill its Detections too
	(*im)->clear(true);
	im = matches.erase(im);
	continue;
      }
    }
    ++im;
  }
}

template <class S>
void
reserveMatches(list<typename S::Match*>& matches,
	       double reserveFraction,
	       int randomNumberSeed) {
  ran::UniformDeviate u;
  if (randomNumberSeed >0) u.Seed(randomNumberSeed);

  for (auto mptr : matches)
    mptr->setReserved( u < reserveFraction );
}

// Return a map of names of non-frozen exposure maps that
// have fewer than minFitExposure fittable Detections using
// them, also gives the number of fittable Detections they have.
template <class S>
map<string, long>
findUnderpopulatedExposures(long minFitExposure,
			    const list<typename S::Match*> matches,
			    const vector<Exposure*> exposures,
			    const vector<typename S::Extension*> extensions,
			    const typename S::Collection& pmc) {
  // First count up useful Detections in each extension:
  vector<long> extnCounts(extensions.size(), 0);
  for (auto mptr : matches)
    if (!mptr->getReserved())
      for (auto dptr : *mptr)
	if (!dptr->isClipped)
	  extnCounts[dptr->catalogNumber]++;

  // Now for each exposure, add up extnCounts of all extensions
  // that depend on a map with name of the exposure.
  map<string, long> bad;
  for (int iExpo=0; iExpo<exposures.size(); iExpo++) {
    if (!exposures[iExpo]) continue; // Not in use
    string expoName = exposures[iExpo]->name;
    if (!pmc.mapExists(expoName) || pmc.getFixed(expoName)) continue; // No free map
    long expoCounts=0;
    for (int iExtn=0; iExtn<extensions.size(); iExtn++) {
      if (extnCounts[iExtn]==0) continue;
      if (pmc.dependsOn(extensions[iExtn]->mapName, expoName))
	expoCounts += extnCounts[iExtn];
    }
    if (expoCounts < minFitExposure)
      bad[expoName] = expoCounts;
  }
  return bad;
}
    

// Fix the parameters of a map, and mark all Detections making
// use of it as clipped so they will not be used in fitting
template <class S>
void
freezeMap(string mapName,
	  const list<typename S::Match*> matches,
	  const vector<typename S::Extension*> extensions,
	  typename S::Collection& pmc) {
  // Nothing to do if map is already fixed or doesn't exist
  if (!pmc.mapExists(mapName) || pmc.getFixed(mapName)) return;
  
  set<long> badExtensions;  // values of extensions using the map
  for (int iExtn=0; iExtn<extensions.size(); iExtn++) {
    if (!extensions[iExtn]) continue; // Not in use
    if (pmc.dependsOn(extensions[iExtn]->mapName, mapName))
      badExtensions.insert(iExtn);
  }

  // Now freeze the relevant map
  pmc.setFixed(mapName);

  // And now clip all of the affected Detections
  if (badExtensions.empty()) return;  // nothing to do
  for (auto mptr : matches)
    if (!mptr->getReserved()) {
      bool recount = false;
      for (auto dptr : *mptr)
	if (badExtensions.count(dptr->catalogNumber)>0) {
	  dptr->isClipped = true;
	  recount = true;
	}
      if (recount) mptr->countFit();
    }
}
  
template <class S>
void
matchCensus(const list<typename S::Match*>& matches) {
  // ??? Choice of whether to count reserved??
  long dcount = 0;
  int dof = 0;
  double chi = 0.;
  double maxdev = 0.;
  for (auto mptr : matches) {
    mptr->countFit();
    dcount += mptr->fitSize();
    chi += mptr->chisq(dof, maxdev);
  }
  cout << "Using " << matches.size() 
       << " matches with " << dcount << " total detections." << endl;
  cout << " chisq " << chi << " / " << dof << " dof maxdev " << maxdev << endl;
}

// Map and clip reserved matches
template <class S>
void
clipReserved(typename S::Align& ca,
	     double clipThresh,
	     double minimumImprovement,
	     bool clipEntireMatch,
	     bool reportToCerr) {

  double oldthresh = 0.;
  int nclip;
  do {
    // Report number of active Matches / Detections in each iteration:
    long int mcount=0;
    long int dcount=0;
    ca.count(mcount, dcount, true, 2);
    double max;
    int dof=0;
    double chisq= ca.chisqDOF(dof, max, true);
    if (reportToCerr) 
      cerr << "Clipping " << mcount << " matches with " << dcount << " detections "
	   << " chisq " << chisq << " / " << dof << " dof,  maxdev " << max 
	   << " sigma" << endl;
      
    double thresh = sqrt(chisq/dof) * clipThresh;
    if (reportToCerr) 
      cerr << "  new clip threshold: " << thresh << " sigma"
	   << endl;
    if (thresh >= max) break;
    if (oldthresh>0. && (1-thresh/oldthresh)<minimumImprovement) break;
    oldthresh = thresh;
    nclip = ca.sigmaClip(thresh, true, clipEntireMatch);
    if (reportToCerr) 
      cerr << "Clipped " << nclip
	   << " matches " << endl;
  } while (nclip>0);
}

// Save fitting results (residual) to output FITS table.
template <class S>
void
saveResults(const list<typename S::Match*>& matches,
	    string outCatalog) {
      // Open the output bintable
  string tablename = S::isAstro? "WCSOut" : "PhotoOut";
  FITS::FitsTable ft(outCatalog, FITS::ReadWrite + FITS::Create, tablename);
  img::FTable outTable = ft.use();;

  // Create vectors that will be put into output table
  // (create union of photo and astro columns)
  vector<int> sequence;
  outTable.addColumn(sequence, "SequenceNumber");
  vector<long> catalogNumber;
  outTable.addColumn(catalogNumber, "Extension");
  vector<long> objectNumber;
  outTable.addColumn(objectNumber, "Object");
  vector<bool> clip;
  outTable.addColumn(clip, "Clip");
  vector<bool> reserve;
  outTable.addColumn(reserve, "Reserve");
  vector<float> xpix;
  outTable.addColumn(xpix, "xPix");
  vector<float> ypix;
  outTable.addColumn(ypix, "yPix");
  vector<float> wtFrac;
  outTable.addColumn(wtFrac, "wtFrac");
  vector<float> color;
  outTable.addColumn(color, "Color");

  // Astro columns:
  vector<float> sigpix;
  vector<float> xrespix;
  vector<float> yrespix;
  vector<float> xworld;
  vector<float> yworld;
  vector<float> sigw;
  vector<float> xresw;
  vector<float> yresw;

  // Photo columns
  vector<float> xExposure;
  vector<float> yExposure;
  vector<float> magOut;
  vector<float> sigMag;
  vector<float> magRes;
  if (S::isAstro) {
    outTable.addColumn(xrespix, "xresPix");
    outTable.addColumn(sigpix, "sigPix");
    outTable.addColumn(yrespix, "yresPix");
    outTable.addColumn(xworld, "xW");
    outTable.addColumn(yworld, "yW");
    outTable.addColumn(sigw, "sigW");
    outTable.addColumn(xresw, "xresW");
    outTable.addColumn(yresw, "yresW");
  } else {
    outTable.addColumn(xExposure, "xExpo");
    outTable.addColumn(yExposure, "yExpo");
    outTable.addColumn(magOut, "magOut");
    outTable.addColumn(sigMag, "sigMag");
    outTable.addColumn(magRes, "magRes");
  }  

  // Cumulative counter for rows written to table:
  long pointCount = 0;
  // Write vectors to table when this many rows accumulate:
  const long WriteChunk = 100000;

  // Write all matches to output catalog, deleting them along the way
  // and accumulating statistics of each exposure.
  // 
  for (auto im = matches.begin(); true; ++im) {
    // First, write vectors to table if they've gotten big or we're at end
    if ( sequence.size() >= WriteChunk || im==matches.end()) {
      long nAdded = sequence.size();
      // Common to photo and astro:
      outTable.writeCells(sequence, "SequenceNumber", pointCount);
      sequence.clear();
      outTable.writeCells(catalogNumber, "Extension", pointCount);
      catalogNumber.clear();
      outTable.writeCells(objectNumber, "Object", pointCount);
      objectNumber.clear();
      outTable.writeCells(clip, "Clip", pointCount);
      clip.clear();
      outTable.writeCells(reserve, "Reserve", pointCount);
      reserve.clear();
      outTable.writeCells(xpix, "xPix", pointCount);
      xpix.clear();
      outTable.writeCells(ypix, "yPix", pointCount);
      ypix.clear();
      outTable.writeCells(wtFrac, "wtFrac", pointCount);
      wtFrac.clear();
      outTable.writeCells(color, "Color", pointCount);
      color.clear();

      if (S::isAstro) {
	outTable.writeCells(xrespix, "xresPix", pointCount);
	xrespix.clear();
	outTable.writeCells(yrespix, "yresPix", pointCount);
	yrespix.clear();
	outTable.writeCells(xworld, "xW", pointCount);
	xworld.clear();
	outTable.writeCells(yworld, "yW", pointCount);
	yworld.clear();
	outTable.writeCells(xresw, "xresW", pointCount);
	xresw.clear();
	outTable.writeCells(yresw, "yresW", pointCount);
	yresw.clear();
	outTable.writeCells(sigpix, "sigPix", pointCount);
	sigpix.clear();
	outTable.writeCells(sigw, "sigW", pointCount);
	sigw.clear();
      } else {
	outTable.writeCells(xExposure, "xExpo", pointCount);
	xExposure.clear();
	outTable.writeCells(yExposure, "yExpo", pointCount);
	yExposure.clear();
	outTable.writeCells(magOut, "magOut", pointCount);
	magOut.clear();
	outTable.writeCells(magRes, "magRes", pointCount);
	magRes.clear();
	outTable.writeCells(sigMag, "sigMag", pointCount);
	sigMag.clear();
      }
      pointCount += nAdded;
    }	// Done flushing the vectors to Table

    if (im==matches.end()) break;

    const typename S::Match& m = **im;
    
    // Get mean values and weights for the whole match
    // (generic call that gets both mag & position)
    double mx, my, mmag;
    double wtot;
    S::matchMean(m, mx,my,mmag, wtot);
    // Calculate number of DOF per Detection coordinate after
    // allowing for fit to centroid:
	    
    int detcount=0;
    for (auto detptr : m) {
      // Save common quantities
      sequence.push_back(detcount);
      ++detcount;
      catalogNumber.push_back(detptr->catalogNumber);
      objectNumber.push_back(detptr->objectNumber);
      clip.push_back(detptr->isClipped);
      reserve.push_back(m.getReserved());
      color.push_back(S::getColor(detptr));

      if (S::isAstro) {
	// Variables needed to save for astrometry
	double xp=0., yp=0., sigma=0., xerrpix=0., yerrpix=0.;
	double xw=0., yw=0., xerrw=0., yerrw=0., sigmaw=0.;
	double wf=0.;
	S::getOutputA(*detptr, mx, my, wtot,
		      xp, yp, sigma, xerrpix, yerrpix,
		      xw, yw, sigmaw, xerrw, yerrw,
		      wf);
	xpix.push_back(xp);
	ypix.push_back(yp);
	sigpix.push_back(sigma);
	xworld.push_back(xw);
	yworld.push_back(yw);
	xrespix.push_back(xerrpix);
	yrespix.push_back(yerrpix);
	// Put world residuals in milliarcsec
	xresw.push_back(xerrw*1000.*DEGREE/ARCSEC);
	yresw.push_back(yerrw*1000.*DEGREE/ARCSEC);
	sigw.push_back(sigmaw*1000.*DEGREE/ARCSEC);
	wtFrac.push_back(wf);
      } else {
	double xDevice=0., yDevice=0., xExpo=0., yExpo=0.;
	double mag=0., magerr=0., sigma=0.;
	double wf=0.;
	S::getOutputP(*detptr, mmag, wtot,
		      xDevice, yDevice, xExpo, yExpo,
		      mag, magerr, sigma,
		      wf);
	xpix.push_back(xDevice);
	ypix.push_back(yDevice);
	xExposure.push_back(xExpo);
	yExposure.push_back(yExpo);
	magOut.push_back(mag);
	magRes.push_back(magerr);
	sigMag.push_back(sigma);
	wtFrac.push_back(wf);
      }
    } // End detection loop
  } // end match loop
}

template<class S>
void
reportStatistics(const list<typename S::Match*>& matches,
		 const vector<Exposure*>& exposures,
		 const vector<typename S::Extension*>& extensions,
		 ostream& os) {
  // Create Accum instances for fitted and reserved Detections on every
  // exposure, plus total accumulator for all reference
  // and all non-reference objects.
  typedef Accum<S> Acc;
  vector<Acc> vaccFit(exposures.size());
  vector<Acc> vaccReserve(exposures.size());
  Acc refAccFit;
  Acc refAccReserve;
  Acc accFit;
  Acc accReserve;

  // Accumulate stats over all detections.
  for (auto mptr : matches) {
    // Get mean values and weights for the whole match
    // (generic call that gets both mag & position)
    double xc, yc, mmag;
    double wtot;
    S::matchMean(*mptr, xc, yc, mmag, wtot);
    if (!S::isAstro) xc = mmag;  // For Photometry put mean mag into xc
    
    // Calculate number of DOF per Detection coordinate after
    // allowing for fit to centroid:
    int nFit = mptr->fitSize();
    double dofPerPt = (nFit > 1) ? 1. - 1./nFit : 0.; // ???

    for (auto dptr : *mptr) {
      // Accumulate statistics for meaningful residuals
      if (dofPerPt >= 0. && !dptr->isClipped) {
	int exposureNumber = extensions[dptr->catalogNumber]->exposure;
	Exposure* expo = exposures[exposureNumber];
	if (mptr->getReserved()) {
	  if (expo->instrument==REF_INSTRUMENT) {
	    refAccReserve.add(dptr, xc, yc, dofPerPt);
	    vaccReserve[exposureNumber].add(dptr, xc, yc, dofPerPt);
	  } else if (expo->instrument==TAG_INSTRUMENT) {
	    // do nothing
	  } else {
	    accReserve.add(dptr, xc, yc, dofPerPt);
	    vaccReserve[exposureNumber].add(dptr, xc, yc, dofPerPt);
	  }
	} else {
	  // Not a reserved object:
	  if (expo->instrument==REF_INSTRUMENT) {
	    refAccFit.add(dptr, xc, yc, dofPerPt);
	    vaccFit[exposureNumber].add(dptr, xc, yc, dofPerPt);
	  } else if (expo->instrument==TAG_INSTRUMENT) {
	    // do nothing
	  } else {
	    accFit.add(dptr, xc, yc, dofPerPt);
	    vaccFit[exposureNumber].add(dptr, xc, yc, dofPerPt);
	  }
	}
      }
    } // end Detection loop
  } // end match loop

  // Sort the exposures by name lexical order
  vector<int> ii;
  for (int i=0; i<exposures.size(); i++)
    if (exposures[i])
      ii.push_back(i);
  std::sort(ii.begin(), ii.end(), NameSorter<Exposure>(exposures));
  
  // Print summary statistics for each exposure
  os << "#   Exposure           | " << Accum<S>::header() << endl;
  for (int i=0; i<ii.size(); i++) {
    int iexp = ii[i];
    os << setw(3) << iexp
       << " " << setw(10) << exposures[iexp]->name
       << " Fit     | " << vaccFit[iexp].summary()
       << endl;
      if (vaccReserve[iexp].n > 0) 
	os << setw(3) << iexp
	   << " " << setw(10) << exposures[iexp]->name
	   << " Reserve | " << vaccReserve[iexp].summary()
	   << endl;
  } // exposure summary loop

  // Output summary data for reference catalog and detections

  cout << "# ------------------------------------------------" << endl;
  if (S::isAstro) {
    // Only Astrometric fits have references
    os << "Reference fit          | " << refAccFit.summary() << endl;
    if (refAccReserve.n>0)
      os << "Reference reserve      | " << refAccReserve.summary() << endl;
  }
  os << "Detection fit          | " << accFit.summary() << endl;
  if (accReserve.n>0)
    os << "Detection reserve      | " << accReserve.summary() << endl;
}

//////////////////////////////////////////////////////////////////
// For those routines differing for Astro/Photo, here is where
// we instantiate the two cases.
//////////////////////////////////////////////////////////////////

#define INSTANTIATE(AP) \
template \
vector<AP::Extension*> \
readExtensions<AP> (img::FTable& extensionTable, \
		    const vector<Instrument*>& instruments,	\
		    const vector<Exposure*>& exposures,		   \
		    const vector<int>& exposureColorPriorities,	     \
		    vector<AP::ColorExtension*>& colorExtensions,    \
		    astrometry::YAMLCollector& inputYAML); \
template int \
findCanonical<AP>(Instrument& instr,	\
		  int iInst,			\
		  vector<Exposure*>& exposures,	  \
		  vector<AP::Extension*>& extensions,	\
		  AP::Collection& pmc); \
template void \
fixMapComponents<AP>(AP::Collection&, \
		     const list<string>&, \
		     const vector<Instrument*>&); \
template void \
createMapCollection<AP>(const vector<Instrument*>& instruments, \
			const vector<Exposure*>& exposures, \
			const vector<AP::Extension*> extensions, \
			astrometry::YAMLCollector& inputYAML,	 \
			AP::Collection& pmc); \
template void							\
whoNeedsColor<AP>(vector<AP::Extension*> extensions); \
template void \
readMatches<AP>(img::FTable& table, \
		list<AP::Match*>& matches, \
		vector<AP::Extension*>& extensions, \
		vector<AP::ColorExtension*>& colorExtensions, \
		const ExtensionObjectSet& skipSet, \
		int minMatches); \
template void \
readObjects<AP>(const img::FTable& extensionTable, \
		const vector<Exposure*>& exposures, \
		vector<AP::Extension*>& extensions, \
		double sysError, \
		double referenceSysError); \
template void \
readColors<AP>(img::FTable extensionTable, \
	       vector<AP::ColorExtension*> colorExtensions);  \
template void  \
purgeNoisyDetections<AP>(double maxError,  \
			 list<AP::Match*>& matches,  \
			 const vector<Exposure*>& exposures,  \
			 const vector<AP::Extension*>& extensions);  \
template void  \
purgeSparseMatches<AP>(int minMatches,  \
		       list<AP::Match*>& matches);  \
template void  \
purgeBadColor<AP>(double minColor, double maxColor,  \
		  list<AP::Match*>& matches);  \
  \
template void  \
reserveMatches<AP>(list<AP::Match*>& matches,  \
		   double reserveFraction,  \
		   int randomNumberSeed);  \
  \
template map<string, long>  \
findUnderpopulatedExposures<AP> (long minFitExposure,  \
				 const list<AP::Match*> matches,  \
				 const vector<Exposure*> exposures,  \
				 const vector<AP::Extension*> extensions,  \
				 const AP::Collection& pmc);  \
template void  \
freezeMap<AP>(string mapName,  \
	      const list<AP::Match*> matches,  \
	      const vector<AP::Extension*> extensions,  \
	      AP::Collection& pmc);  \
template void  \
matchCensus<AP>(const list<AP::Match*>& matches);  \
template void  \
clipReserved<AP>(AP::Align& ca,  \
		 double clipThresh,  \
		 double minimumImprovement,  \
		 bool clipEntireMatch,  \
		 bool reportToCerr); \
template void \
saveResults<AP>(const list<AP::Match*>& matches, \
		string outCatalog); \
template void \
reportStatistics<AP>(const list<AP::Match*>& matches,  \
		     const vector<Exposure*>& exposures,  \
		     const vector<AP::Extension*>& extensions,  \
		     ostream& os);

INSTANTIATE(Astro)
INSTANTIATE(Photo)
