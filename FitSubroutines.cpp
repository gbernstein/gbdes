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

template <class T>
void
fixMapComponents(T& pmc,
		 const list<string>& fixMapList,
		 const vector<Instrument*>& instruments) {
  set<string> fixTheseMaps;
  for (auto iName : pmc.allMapNames()) {
    if (stringstuff::regexMatchAny(fixMapList, iName))
      fixTheseMaps.insert(iName);
  }
  // Add names of all devices of instruments on the fixMapList
  for (auto instptr : instruments) {
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

// Instantiate our 2 cases
template
void fixMapComponents<astrometry::PixelMapCollection>(astrometry::PixelMapCollection&,
						      const list<string>&,
						      const vector<Instrument*>&);
/**
template
void fixMapComponents<photometry::PhotoMapCollection>(photometry::PhotoMapCollection&,
						      const list<string>&,
						      const vector<Instrument*>&);
**/
template <class TExtn, class TColor>
vector<TExtn*>
readExtensions(img::FTable& extensionTable,
	       const vector<Instrument*>& instruments,
	       const vector<Exposure*>& exposures,
	       const vector<int>& exposureColorPriorities,
	       vector<TColor*>& colorExtensions,
	       astrometry::YAMLCollector& inputYAML) {
      
  vector<TExtn*> extensions(extensionTable.nrows(), 0);
  colorExtensions = vector<TColor*>(extensionTable.nrows(), 0);
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
      auto* ce = new TColor;
      ce->priority = colorPriority;
      colorExtensions[i] = ce;
    }
	
    if (!exposures[iExposure]) continue;
    
    TExtn* extn = new TExtn;
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
      extn->basemapName = astrometry::IdentityMap().getName();
    } else {
      // Real instrument, make a map combining its exposure with its Device map:
      astrometry::YAMLCollector::Dictionary d;
      d["INSTRUMENT"] = instruments[expo.instrument]->name;
      d["DEVICE"] = instruments[expo.instrument]->deviceNames.nameOf(extn->device);
      d["EXPOSURE"] = expo.name;
      d["BAND"] = instruments[expo.instrument]->band;
      // This will be the name of the extension's final WCS and photometry map:
      extn->wcsName = d["EXPOSURE"] + "/" + d["DEVICE"];
      // And this is the name of the PixelMap underlying the WCS:
      extn->basemapName = extn->wcsName + "/base";
	  
      if (!inputYAML.addMap(extn->basemapName,d)) {
	cerr << "Input YAML files do not have complete information for map "
	     << extn->basemapName
	     << endl;
	exit(1);
      }
    }
  }  // End extension loop
}

#define INSTANTIATE(ns)				            \
template                                                    \
vector<ExtensionBase<ns::SubMap, ns::Detection>*>           \
 readExtensions<ExtensionBase<ns::SubMap, ns::Detection>,   \
                ColorExtensionBase<ns::Match>>              \
              (img::FTable& extensionTable,		    \
	       const vector<Instrument*>& instruments,      \
	       const vector<Exposure*>& exposures,          \
	       const vector<int>& exposureColorPriorities,  \
	       vector<ColorExtensionBase<ns::Match>*>& colorExtensions,  \
	       astrometry::YAMLCollector& inputYAML);

INSTANTIATE(astrometry)
INSTANTIATE(photometry)
