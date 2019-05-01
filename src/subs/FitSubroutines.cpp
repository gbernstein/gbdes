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
#include "Stopwatch.h"

using astrometry::PMMatch;
using astrometry::PMDetection;
using astrometry::AstrometryError;
using astrometry::WCS_UNIT;
using astrometry::RESIDUAL_UNIT;
using astrometry::PM_UNIT;

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
  for (auto& i: ls) {
    if (i.empty()) continue;
    list<string> ls2 = split(i, '=');
    if (ls2.size() != 2) {
      cerr << errorDescription << " has bad translation spec: " << i << endl;
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

// Read info from fields table, write to a new output file
void
readFields(string inputTables,
	   string outCatalog,
	   NameIndex& fieldNames,
	   vector<astrometry::SphericalCoords*>& fieldProjections,
	   vector<double>& fieldEpochs,
	   double defaultEpoch) {
  const double MINIMUM_EPOCH=1900.;
  
  FITS::FitsTable in(inputTables, FITS::ReadOnly, "Fields");
  FITS::FitsTable out(outCatalog,
		      FITS::ReadWrite + FITS::OverwriteFile + FITS::Create,
		      "Fields");
  img::FTable ft = in.extract();
  out.adopt(ft);
  vector<double> ra;
  vector<double> dec;
  vector<string> name;
  ft.readCells(name, "Name");
  ft.readCells(ra, "RA");
  ft.readCells(dec, "Dec");
  for (int i=0; i<name.size(); i++) {
    spaceReplace(name[i]);
    fieldNames.append(name[i]);
    astrometry::Orientation orient(astrometry::SphericalICRS(ra[i]*WCS_UNIT, dec[i]*WCS_UNIT));
    fieldProjections.push_back( new astrometry::Gnomonic(orient));
  }
  if (ft.hasColumn("Epoch")) {
    ft.readCells(fieldEpochs, "Epoch");
    // Set reference epochs for each field to defaultEpoch if
    // they have a signal value like 0 or -1 
    for (int i=0; i<fieldEpochs.size(); i++) {
      if (fieldEpochs[i] < MINIMUM_EPOCH && defaultEpoch>MINIMUM_EPOCH) {
	// Update to reference epoch if it's sensible
	fieldEpochs[i] = defaultEpoch;
	ft.writeCell(defaultEpoch, "Epoch", i);
      }
    }
  } else {
    // Assign default epoch to every field
    fieldEpochs.clear();
    fieldEpochs.resize(name.size(), defaultEpoch);
    // And save back to output if this is a sensible default
    if (defaultEpoch > MINIMUM_EPOCH) {
      ft.addColumn(fieldEpochs, "Epoch");
    }
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
    if (instrumentNumber >= instruments.size()) {
      int oldsize=instruments.size();
      instruments.resize(instrumentNumber+1);
      for ( ; oldsize < instruments.size(); oldsize++)
	instruments[oldsize] = 0;
    }
    spaceReplace(instrumentName);

    // Send every instrument extension to output file
    FITS::Flags outFlags = FITS::ReadWrite+FITS::Create;
    if (!outputCatalogAlreadyOpen) {
      outFlags = outFlags + FITS::OverwriteFile;
      outputCatalogAlreadyOpen = true;
    }
    FITS::FitsTable out(outCatalog, outFlags, -1);
    out.setName("Instrument");
    img::FTable ff=ft.extract();
    out.adopt(ff);

    if (regexMatchAny(useInstrumentList, instrumentName))  {
      // This is an instrument we will use.  Make an Instrument instance
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
	      const vector<double>& fieldEpochs,
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
  vector<double> mjd;
  vector<double> apcorr;
  vector<string> epoch;

  vector<vector<double>> observatory;
  vector<vector<double>> astrometricCovariance;
  vector<double> photometricVariance;
  vector<double> weight;
  vector<double> magWeight;

  ff.readCells(names, "Name");
  ff.readCells(ra, "RA");
  ff.readCells(dec, "Dec");
  ff.readCells(fieldNumber, "fieldNumber");
  ff.readCells(instrumentNumber, "InstrumentNumber");
  ff.readCells(airmass, "Airmass");
  ff.readCells(exptime, "Exptime");

  // Columns that might not always be present:
  try {
    ff.readCells(mjd, "MJD");
  } catch (img::FTableError& e) {
    mjd.clear();
  }
  try {
    ff.readCells(epoch, "Epoch");
  } catch (img::FTableError& e) {
    epoch.clear();
  }
  try {
    ff.readCells(apcorr, "apcorr");
  } catch (img::FTableError& e) {
    apcorr.clear();
  }
  try {
    ff.readCells(weight, "weight");
  } catch (img::FTableError& e) {
    weight.clear();
  }
  try {
    ff.readCells(magWeight, "magWeight");
  } catch (img::FTableError& e) {
    magWeight.clear();
  }
  try {
    ff.readCells(photometricVariance, "sysMagVar");
  } catch (img::FTableError& e) {
    photometricVariance.clear();
  }
  try {
    ff.readCells(astrometricCovariance, "sysAstCov");
  } catch (img::FTableError& e) {
    astrometricCovariance.clear();
  }
  try {
    ff.readCells(observatory, "observatory");
  } catch (img::FTableError& e) {
    observatory.clear();
  }

  // Initialize our output arrays to not-in-use values
  vector<Exposure*> exposures(names.size(), nullptr);
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
    bool useThisExposure = (instrumentNumber[i]==REF_INSTRUMENT ||
			    instrumentNumber[i]==PM_INSTRUMENT)
      && useReferenceExposures;
    // or if it's a normal exposure using an instrument that has been included
    useThisExposure = useThisExposure || (instrumentNumber[i]>=0 && 
					  instruments[instrumentNumber[i]]);

    if (useThisExposure) {
      // The projection we will use for this exposure:
      astrometry::Gnomonic gn(astrometry::Orientation(astrometry::SphericalICRS(ra[i]*WCS_UNIT,
										dec[i]*WCS_UNIT)));
      auto expo = new Exposure(names[i],gn);
      expo->field = fieldNumber[i];
      expo->instrument = instrumentNumber[i];
      expo->airmass = airmass[i];
      expo->exptime = exptime[i];
      exposures[i] = expo;
      if (!mjd.empty()) {
	expo->mjd = mjd[i];
	// Also calculate time in years after field's reference epoch
	astrometry::UT ut;
	ut.setMJD(expo->mjd);
	// Note that getTTyr returns years since J2000.
	expo->pmTDB = ut.getTTyr() - (fieldEpochs[expo->field] - 2000.);
      }
      if (!epoch.empty())
	expo->epoch = epoch[i];
      if (!apcorr.empty())
	expo->apcorr = apcorr[i];
      if (!weight.empty())
	expo->weight = weight[i];
      if (!magWeight.empty())
	expo->magWeight = magWeight[i];
      if (!observatory.empty())
	for (int j=0; j<3; j++) expo->observatory[j] = observatory[i][j];
      if (!astrometricCovariance.empty()) {
	expo->astrometricCovariance(0,0) = astrometricCovariance[i][0];
	expo->astrometricCovariance(1,1) = astrometricCovariance[i][1];
	expo->astrometricCovariance(0,1) = astrometricCovariance[i][2];
	expo->astrometricCovariance(1,0) = astrometricCovariance[i][2];
      }
      if (!photometricVariance.empty())
	expo->photometricVariance = photometricVariance[i];
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
      
  vector<typename S::Extension*> extensions(extensionTable.nrows(), nullptr);
  colorExtensions = vector<typename S::ColorExtension*>(extensionTable.nrows(), nullptr);
  int processed=0;

#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic,100)
#endif
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
    
#ifdef _OPENMP
#pragma omp critical(processed)
#endif
    ++processed; 
    typename S::Extension* extn = new typename S::Extension;
    extn->exposure = iExposure;
    const Exposure& expo = *exposures[iExposure];

    bool isReference = (expo.instrument == REF_INSTRUMENT ||
			expo.instrument == PM_INSTRUMENT);

    int iDevice;
    extensionTable.readCell(iDevice, "Device", i);
    if (processed%1000==0) {
      cerr << "...Extn " << i << "/" << extensions.size() 
	   << " " << expo.name << endl;
    }
    extn->device = iDevice;
    extn->airmass = expo.airmass;
    extn->apcorr = expo.apcorr;
    extn->magshift = +2.5*log10(expo.exptime);

    // Create the starting WCS for the exposure
    string s;
    extensionTable.readCell(s, "WCSIn", i);
    if (stringstuff::nocaseEqual(s, "_ICRS")) {
      // Create a Wcs that just takes input as RA and Dec in degrees;
      astrometry::IdentityMap identity;
      astrometry::SphericalICRS icrs;
      extn->startWcs = new astrometry::Wcs(&identity, icrs, "ICRS_degrees", WCS_UNIT);
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
      // Add Epoch to the dictionary if it is not empty
      if (!expo.epoch.empty())
	d["EPOCH"] = expo.epoch;
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
      // Add Epoch to the dictionary if it is not empty
      if (!expo.epoch.empty())
	d["EPOCH"] = expo.epoch;
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
	    int minMatches,
	    bool usePM,
	    double parallaxPrior) {

  vector<int> seq;
  vector<LONGLONG> extn;
  vector<LONGLONG> obj;
  table.readCells(seq, "sequenceNumber");
  table.readCells(extn, "extension");
  table.readCells(obj, "object");

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
	typename S::Match* m=nullptr;
	for (int j=0; j<matchExtns.size(); j++) {
	  if (matchExtns[j]<0) continue;  // Skip detections starved of their color
	  auto d = new typename S::Detection;  
	  d->catalogNumber = matchExtns[j];
	  d->objectNumber = matchObjs[j];
	  extensions[matchExtns[j]]->keepers.insert(std::pair<long,typename S::Detection*>
						    (matchObjs[j], d));
	  if (m)
	    m->add(d);
	  else {
	    // Make a new Match object from this detection.
	    // It might be a PMMatch if this is an appropriate catalog,
	    // in which case it will get the parallaxPrior too.
	    m = S::makeNewMatch(d, usePM, parallaxPrior);
	  }
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
		     const Exposure* e,
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
		     bool isTag) {
  d->xpix = getTableDouble(table, xKey, -1, xColumnIsDouble, irow);
  d->ypix = getTableDouble(table, yKey, -1, yColumnIsDouble, irow);

  // Get coordinates and transformation matrix
  startWcs->toWorld(d->xpix, d->ypix, d->xw, d->yw);  // no color in startWCS
  auto dwdp = startWcs->dWorlddPix(d->xpix, d->ypix);

  if (isTag) {
    d->invCov.setZero();
  } else {
    astrometry::Matrix22 cov(0.);
    if (xyErrKeys.size()==1) {
      // We have a single pixel error, diagonal
      double sigma = getTableDouble(table, xyErrKeys[0], -1, errorColumnIsDouble,irow);
      cov(0,0) = sigma*sigma;
      cov(1,1) = sigma*sigma;
    } else if (xyErrKeys.size()==3) {
      // We have three components of an ellipse, giving x^2, y^2, xy values:
      cov(0,0) = getTableDouble(table, xyErrKeys[0], -1, errorColumnIsDouble,irow);
      cov(1,1) = getTableDouble(table, xyErrKeys[1], -1, errorColumnIsDouble,irow);
      cov(0,1) = getTableDouble(table, xyErrKeys[2], -1, errorColumnIsDouble,irow);
      cov(1,0) = cov(0,1);
    } else {
      for (auto s : xyErrKeys)
	cerr << " --- " << s;
      cerr << " ---" << endl;
      throw AstrometryError("Invalid number of xyErrKeys passed to fillDetection");
    }
    cov = dwdp * cov * dwdp.transpose();  // Note this converts to world units (degrees)
    cov += e->astrometricCovariance; 
    d->invCov = cov.inverse();
    d->fitWeight = e->weight;
  }

  if (dynamic_cast<const PMMatch*>(d->itsMatch)) {
    // Build projection matrix if this Detection is being used in a PMMatch
    d->buildProjector(e->pmTDB, e->observatory, &fieldProjection);
  }

}

// This one reads a full 5d stellar solution
astrometry::PMDetection*
Astro::makePMDetection(astrometry::Detection* d, const Exposure* e,
		img::FTable& table, long irow,
		string xKey, string yKey,
		string pmRaKey, string pmDecKey, string parallaxKey, string pmCovKey,
		bool xColumnIsDouble, bool yColumnIsDouble,
		bool errorColumnIsDouble,
		const astrometry::PixelMap* startWcs) {

  auto out = new astrometry::PMDetection;
  out->catalogNumber = d->catalogNumber;
  out->objectNumber = d->objectNumber;
  
  out->xpix = getTableDouble(table, xKey, -1, xColumnIsDouble, irow);
  out->ypix = getTableDouble(table, yKey, -1, yColumnIsDouble, irow);
  // Read in the reference solution
  double pmRA, pmDec, parallax;
  table.readCell(pmRA, pmRaKey, irow);
  table.readCell(pmDec, pmDecKey, irow);
  table.readCell(parallax, parallaxKey, irow);

  // We will want to shift the inputs to move their reference time from
  // the exposure's (i.e. catalog's) reference to that of the field,
  // which is defined as 0 here.
  double epochShift = -e->pmTDB;

  // Shift the RA and Dec according to proper motion, including
  // factors for unit difference, and cos(dec) on the RA.
  double cosdec = cos(out->ypix * WCS_UNIT);
  out->xpix += epochShift * pmRA * PM_UNIT / WCS_UNIT / cosdec;
  out->ypix += epochShift * pmDec * PM_UNIT / WCS_UNIT;

  // Get coordinates and transformation matrix to world coordinates
  startWcs->toWorld(out->xpix, out->ypix, out->xw, out->yw);  // no color in startWCS
  auto dwdp = startWcs->dWorlddPix(out->xpix, out->ypix);
  // Add cos(dec) factors, since all error values are assumed for ra*cos(dec).
  // and we are assuming that PM catalogs are in RA/Dec, in WCS_UNIT (degrees, see Units.h)
  dwdp.col(0) /= cosdec;

  // Fill in the PM central values
  out->pmMean[astrometry::X0] = out->xw;
  out->pmMean[astrometry::Y0] = out->yw;
  out->pmMean[astrometry::VX] = pmRA;
  out->pmMean[astrometry::VY] = pmDec;
  out->pmMean[astrometry::PAR] = parallax;

  // Read the covariance matrix
  vector<double> tmp;
  table.readCell(tmp, pmCovKey, irow);
  astrometry::PMCovariance pmCov;  // Covariance from catalog
  astrometry::PMCovariance dwdp5(0.);  // 5d transformation matrix to world coords
  int k=0;
  for (int i=0; i<5; i++) {
    dwdp5(i,i) = 1.;
    for (int j=0; j<5; j++, k++) {
      pmCov(i,j) = tmp[k];
    }
  }

  if (epochShift != 0.) {
    // Transform covariance matrix for shift in reference time
    astrometry::PMCovariance shift(0.);
    for (int i=0; i<5; i++) shift(i,i) = 1.;
    shift(astrometry::X0, astrometry::VX) = epochShift * PM_UNIT / WCS_UNIT;
    shift(astrometry::Y0, astrometry::VY) = epochShift * PM_UNIT / WCS_UNIT;
    pmCov = shift * pmCov * shift.transpose();
  }
  
  // Now apply transformation from RA/Dec to world coordinate system.
  // Subtlety here: X0 and Y0 are in world coords.
  // The proper motion and parallax are kept in ICRS.
  // So only X0, Y0 are altered by spherical projection.
  dwdp5.subMatrix(0,2,0,2) = dwdp;
  astrometry::PMCovariance wCov = (dwdp5 * pmCov * dwdp5.transpose());
  // Save inverse
#ifdef USE_EIGEN
  pmCov.setIdentity();
  out->pmInvCov = wCov.ldlt().solve(pmCov);
#else
  out->pmInvCov = wCov.inverse();
#endif
  out->fitWeight = e->weight;

  // Now fill in the Detection 2d (inverse) covariance in case we end
  // up using the as a single observation at the field's reference epoch.
  astrometry::Matrix22 cov22 = pmCov.subMatrix(0,2,0,2).inverse();
  // Add any extra covariance associated with reference catalogs.
  cov22 += e->astrometricCovariance;
  out->invCov = cov22.inverse();

  return out;
}

void
Astro::handlePMDetection(astrometry::PMDetection* pmd, Astro::Detection* d) {
  auto mm = d->itsMatch;
  if (dynamic_cast<astrometry::PMMatch*>(mm)) {
    // If this Detection is being used in a PMMatch,
    // replace plain old Detection with this PMDetection.
    mm->remove(d);
    delete d;
    mm->add(pmd);
  } else {
    // PMDetection is being used in non-PM Match.  Slice it down
    // to pure Detection info, replace old detection with it.
    auto dd = new typename Astro::Detection(*pmd);
    mm->remove(d);
    delete d;
    mm->add(dd);
    delete pmd;
  }
}
 
astrometry::Match*
Astro::makeNewMatch(Astro::Detection* d, bool usePM, double parallaxPrior) {
  if (usePM) {
    // Make a PMMatch and give it parallax prior
    auto mpm = new astrometry::PMMatch(d);
    if (parallaxPrior>0.) 
      mpm->setParallaxPrior(parallaxPrior);
    return mpm;
  } else {
    // Plain old Match
    return new Match(d);
  }
}

// And a routine to get photometric information too
void
Photo::fillDetection(Photo::Detection* d,
		     const Exposure* e,
		     astrometry::SphericalCoords& fieldProjection, 
		     img::FTable& table, long irow,
		     string xKey, string yKey,
		     vector<string>& xyErrKeys,
		     string magKey, string magErrKey,
		     int magKeyElement, int magErrKeyElement,
		     bool xColumnIsDouble, bool yColumnIsDouble,
		     bool errorColumnIsDouble,
		     bool magColumnIsDouble, bool magErrColumnIsDouble,
		     double magshift,
		     const astrometry::PixelMap* startWcs,
		     bool isTag) {
  d->args.xDevice = getTableDouble(table, xKey, -1, xColumnIsDouble, irow);
  d->args.yDevice = getTableDouble(table, yKey, -1, yColumnIsDouble, irow);
  startWcs->toWorld(d->args.xDevice, d->args.yDevice,
		    d->args.xExposure, d->args.yExposure);

  // Get the mag input and its error
  d->magIn = getTableDouble(table, magKey, magKeyElement, magColumnIsDouble,irow)
    + magshift;
  double sigma = getTableDouble(table, magErrKey, magErrKeyElement, magErrColumnIsDouble,irow);

  if ( isnan(d->magIn) || isinf(d->magIn)) {
    cerr << "** NaN input at row " << irow << endl;
    d->isClipped = true;
  
    d->magIn = 0.;
  }
  // Map to output and estimate output error
  d->magOut = d->map->forward(d->magIn, d->args);

  sigma = std::sqrt(e->photometricVariance + sigma*sigma);
  d->sigma = sigma;
  sigma *= d->map->derivative(d->magIn, d->args);

  // no clips on tags
  double wt = isTag ? 0. : pow(sigma,-2.);
  d->clipsq = wt;
  d->wt = wt * e->magWeight;
}

// Read each Extension's objects' data from it FITS catalog
// and place into Detection structures.
template <class S>
void readObjects(const img::FTable& extensionTable,
		 const vector<Exposure*>& exposures,
		 vector<typename S::Extension*>& extensions,
		 vector<astrometry::SphericalCoords*> fieldProjections) {

  // Should be safe to multithread this loop as different threads write
  // only to distinct parts of memory.  Protect the FITS table read though.
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic,60) num_threads(CATALOG_THREADS)
#endif
  for (int iext = 0; iext < extensions.size(); iext++) {
    if (!extensions[iext]) continue; // Skip unused 

    // Relevant structures for this extension
    typename S::Extension& extn = *extensions[iext];
    Exposure& expo = *exposures[extn.exposure];
    if (extn.keepers.empty()) continue; // or useless

    bool pmCatalog = S::isAstro && expo.instrument==PM_INSTRUMENT;
    bool isTag = expo.instrument == TAG_INSTRUMENT;

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
    vector<string> xyErrKeys; // Holds name(s) of posn error keys
    string idKey;
    // For PM catalogs:
    string pmRaKey;
    string pmDecKey;
    string parallaxKey;
    string pmCovKey;
    // For magnitudes
    string magKey;
    string magErrKey;
    int magKeyElement;
    int magErrKeyElement;
    
    extensionTable.readCell(hduNumber, "extension", iext);
    extensionTable.readCell(xKey, "xKey", iext);
    extensionTable.readCell(yKey, "yKey", iext);
    extensionTable.readCell(idKey, "idKey", iext);

    // Keep track of what rows we need to read from this catalog
    vector<string> neededColumns;
    bool useRows = stringstuff::nocaseEqual(idKey, "_ROW");
    if (!useRows)
      neededColumns.push_back(idKey);
    neededColumns.push_back(xKey);
    neededColumns.push_back(yKey);

    if (S::isAstro) {
      if (!isTag && !pmCatalog) {
	// See if we have the keys for an error ellipse

	if (extensionTable.hasColumn("errXXKey") &&
	    extensionTable.hasColumn("errYYKey") &&
	    extensionTable.hasColumn("errXYKey")) {
	  string errXXKey, errYYKey, errXYKey;
	  extensionTable.readCell(errXXKey, "errXXKey", iext);
	  extensionTable.readCell(errYYKey, "errYYKey", iext);
	  extensionTable.readCell(errXYKey, "errXYKey", iext);
	  if (!(errXXKey.empty() || errYYKey.empty() || errXYKey.empty())) {
	    // We have all the keys for an error ellipse
	    xyErrKeys.push_back(errXXKey);
	    xyErrKeys.push_back(errYYKey);
	    xyErrKeys.push_back(errXYKey);
	  }
	}

	if (xyErrKeys.empty()) {
	  // We did not get an error ellipse.  We will need a
	  // scalar error value
	  string errKey;
	  extensionTable.readCell(errKey, "errKey", iext);
	  xyErrKeys.push_back(errKey);
	}

	for (auto s : xyErrKeys)
	  neededColumns.push_back(s);

      } else if (pmCatalog) {
	// Need PM information
	extensionTable.readCell(pmRaKey, "pmRaKey", iext);
	extensionTable.readCell(pmDecKey, "pmDecKey", iext);
	extensionTable.readCell(pmCovKey, "pmCovKey", iext);
	extensionTable.readCell(parallaxKey, "parallaxKey", iext);
	neededColumns.push_back(pmRaKey);
	neededColumns.push_back(pmDecKey);
	neededColumns.push_back(pmCovKey);
	neededColumns.push_back(parallaxKey);
      }
	
    } else {
      // For photometry:
      extensionTable.readCell(magKey, "magKey", iext);
      extensionTable.readCell(magErrKey, "magErrKey", iext);
      // Be willing to get an element of array-valued bintable cell
      // for mag or magerr.  Syntax would be
      // MAGAPER[4]  to get 4th (0-indexed) element of MAGAPER column
      magKeyElement = elementNumber(magKey);
      magErrKeyElement = elementNumber(magErrKey);
      neededColumns.push_back(magKey);
      neededColumns.push_back(magErrKey);
    }
    
    const typename S::SubMap* sm=extn.map;

    astrometry::Wcs* startWcs = extn.startWcs;

    if (!startWcs) {
      cerr << "Failed to find initial Wcs for exposure " << expo.name
	   << endl;
      exit(1);
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

    bool xColumnIsDouble = isDouble(ff, xKey, -1);
    bool yColumnIsDouble = isDouble(ff, yKey, -1);
    bool magColumnIsDouble;
    bool magErrColumnIsDouble;
    bool errorColumnIsDouble;
    double magshift = extn.magshift;

    // Get fieldProjection for this catalog
    astrometry::SphericalCoords* fieldProjection =
      fieldProjections[expo.field]->duplicate();
    
    if (S::isAstro) {
      errorColumnIsDouble = isDouble(ff, pmCatalog ? pmCovKey : xyErrKeys[0], -1);
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

      if (pmCatalog) {
	// Need to read data differently from catalog with
	// full proper motion solution.  Get PMDetection.
	auto pmd = S::makePMDetection(d, &expo,
				      ff, irow,
				      xKey, yKey,
				      pmRaKey, pmDecKey, parallaxKey, pmCovKey,
				      xColumnIsDouble,yColumnIsDouble,
				      errorColumnIsDouble, 
				      startWcs);

	// This routine will either replace the plain old Detection in
	// a PMMatch with this PMDetection; or if it is part of a plain
	// match, it will slice the PMDetection down to a Detection.
	S::handlePMDetection(pmd, d);
      } else {
	// Normal photometric or astrometric entry
	S::fillDetection(d, &expo, *fieldProjection,
			 ff, irow,
			 xKey, yKey, xyErrKeys,
			 magKey, magErrKey,
			 magKeyElement, magErrKeyElement,
			 xColumnIsDouble,yColumnIsDouble,
			 errorColumnIsDouble, magColumnIsDouble, magErrColumnIsDouble,
			 magshift,
			 startWcs, isTag);
      }
    } // End loop over catalog objects

    if (fieldProjection) delete fieldProjection;

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
    extensionTable.readCell(hduNumber, "extension", iext);
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
// Also deletes anything that is already marked as clipped
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
      // Keep it if error is small or if it's from a tag or reference "instrument"
      if ( (d->isClipped || d->getSigma() > maxError) &&
	   exposures[ extensions[d->catalogNumber]->exposure]->instrument >= 0) {
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
      if (color != astrometry::NODATA && (color < minColor || color > maxColor)) {
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
  ran::UniformDeviate<double> u;
  if (randomNumberSeed >0) u.seed(randomNumberSeed);

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
matchCensus(const list<typename S::Match*>& matches, ostream& os) {
  long dcount = 0;
  int dof = 0;
  double chi = 0.;
  double maxdev = 0.;
  for (auto mptr : matches) {
    mptr->countFit();
    dcount += mptr->fitSize();
    chi += mptr->chisq(dof, maxdev);
  }
  os << "# Using " << matches.size() 
       << " matches with " << dcount << " total detections." << endl;
  os << "#  chisq " << chi << " / " << dof << " dof maxdev " << maxdev << endl;
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
      
    double thresh = sqrt(chisq/dof) * clipThresh; // ??? expected chisq instead of dof?
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

// Save photometric fitting results (residual) to output FITS table.
void
Photo::saveResults(const list<Match*>& matches,
		   string outCatalog) {
  // Open the output bintable
  string tablename = "PhotoOut";
  FITS::FitsTable ft(outCatalog, FITS::ReadWrite + FITS::Create, tablename);
  img::FTable outTable = ft.use();;

  // Create vectors that will be put into output table
  // (create union of photo and astro columns)
  vector<int> sequence;
  outTable.addColumn(sequence, "sequenceNumber");
  vector<long> catalogNumber;
  outTable.addColumn(catalogNumber, "extension");
  vector<long> objectNumber;
  outTable.addColumn(objectNumber, "object");
  vector<bool> clip;
  outTable.addColumn(clip, "clip");
  vector<bool> reserve;
  outTable.addColumn(reserve, "reserve");
  vector<float> xpix;
  outTable.addColumn(xpix, "xPix");
  vector<float> ypix;
  outTable.addColumn(ypix, "yPix");
  vector<float> wtFrac;
  outTable.addColumn(wtFrac, "wtFrac");
  vector<float> color;
  outTable.addColumn(color, "color");

  vector<float> xExposure;
  outTable.addColumn(xExposure, "xExpo");
  vector<float> yExposure;
  outTable.addColumn(yExposure, "yExpo");
  vector<float> magOut;
  outTable.addColumn(magOut, "magOut");
  vector<float> sigMag;
  outTable.addColumn(sigMag, "sigMag");
  vector<float> magRes;
  outTable.addColumn(magRes, "magRes");

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

      outTable.writeCells(sequence, "sequenceNumber", pointCount);
      sequence.clear();
      outTable.writeCells(catalogNumber, "extension", pointCount);
      catalogNumber.clear();
      outTable.writeCells(objectNumber, "object", pointCount);
      objectNumber.clear();
      outTable.writeCells(clip, "clip", pointCount);
      clip.clear();
      outTable.writeCells(reserve, "reserve", pointCount);
      reserve.clear();
      outTable.writeCells(xpix, "xPix", pointCount);
      xpix.clear();
      outTable.writeCells(ypix, "yPix", pointCount);
      ypix.clear();
      outTable.writeCells(wtFrac, "wtFrac", pointCount);
      wtFrac.clear();
      outTable.writeCells(color, "color", pointCount);
      color.clear();

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

      pointCount += nAdded;
    }	// Done flushing the vectors to Table

    if (im==matches.end()) break;

    const Match& m = **im;
    
    // Get mean mag and weights for the whole match
    double mmag, wtot;
    m.getMean(mmag,wtot);

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
      color.push_back(detptr->args.color);

      // Photometry quantities:
      xpix.push_back(detptr->args.xDevice);
      ypix.push_back(detptr->args.yDevice);
      xExposure.push_back(detptr->args.xExposure);
      yExposure.push_back(detptr->args.yExposure);
      magOut.push_back(detptr->magOut);
      magRes.push_back(mmag==0. ? 0. : detptr->magOut - mmag);
      double sigma = detptr->clipsq;
      sigma = sigma>0. ? 1./sqrt(sigma) : 0.;
      sigMag.push_back(sigma);
      wtFrac.push_back(detptr->wt/wtot);
    } // End detection loop
  } // end match loop
}

/***** Astro version ***/
// Save fitting results (residuals) to output FITS tables.
void
Astro::saveResults(const list<Astro::Match*>& matches,
		   string outCatalog) {

  // Open the output bintable for pure position residuals
  string tableName = "WCSOut";
  string pmTableName = "PMOut";
  string starTableName = "StarCat";
  FITS::FitsTable ft(outCatalog, FITS::ReadWrite + FITS::Create, tableName);
  img::FTable outTable = ft.use();;

  // pmTable and starTable will be created after gathering all data

  // Create vectors that will be put into output tables
  vector<long> matchID;
  outTable.addColumn(matchID, "matchID");
  vector<long> catalogNumber;
  outTable.addColumn(catalogNumber, "extension");
  vector<long> objectNumber;
  outTable.addColumn(objectNumber, "object");
  vector<bool> clip;
  outTable.addColumn(clip, "clip");
  vector<bool> reserve;
  outTable.addColumn(reserve, "reserve");
  vector<bool> hasPM;  // Was this input a full PM estimate?
  outTable.addColumn(hasPM, "hasPM");
  vector<float> color;
  outTable.addColumn(color, "color");

  vector<float> xpix;
  outTable.addColumn(xpix, "xPix");
  vector<float> ypix;
  outTable.addColumn(ypix, "yPix");
  vector<float> sigpix;  // Circularized total error in pixels
  outTable.addColumn(sigpix, "sigPix");
  vector<float> xrespix;
  outTable.addColumn(xrespix, "xresPix");
  vector<float> yrespix;
  outTable.addColumn(yrespix, "yresPix");
  vector<float> xworld;
  outTable.addColumn(xworld, "xW");
  vector<float> yworld;
  outTable.addColumn(yworld, "yW");
  vector<float> xresw;
  outTable.addColumn(xresw, "xresW");
  vector<float> yresw;
  outTable.addColumn(yresw, "yresW");

  // Give full error ellipse
  vector<vector<float>> covTotalW;
  outTable.addColumn(covTotalW, "covTotalW");

  // The true chisq for this detection, and the expected value
  vector<float> chisq;
  outTable.addColumn(chisq, "chisq");
  vector<float> chisqExpected;
  outTable.addColumn(chisqExpected, "chisqExpected");

  // These are quantities we will want for PMDetections
  vector<long> pmMatchID;
  vector<long> pmCatalogNumber;
  vector<long> pmObjectNumber;
  vector<bool> pmClip;
  vector<bool> pmReserve;
  vector<vector<float>> pmMean;
  vector<vector<float>> pmInvCov;
  vector<float> pmChisq;
  vector<float> pmChisqExpected;
  // (Not putting residuals here since PM's will be findable
  //  by matchID in the star table)

  // These are quantities we will want to put into our output PM catalog
  vector<long> starMatchID;
  vector<bool> starReserve;  // Is this star reserved from fit?
  vector<float> starColor;    // Color for this star
  vector<int> starPMCount;   // How many PMDetections in it were fit?
  vector<int> starDetCount;  // How many non-PM Detections in it were fit?
  vector<int> starClipCount; // How many detections were clipped?
  vector<int> starDOF;       // Number of DOF in PM fit
  vector<float> starChisq;   // Total chisq
  vector<float> starX;       // The solution
  vector<float> starY;
  vector<float> starPMx;
  vector<float> starPMy;
  vector<float> starParallax;
  vector<vector<float>> starInvCov;  // flattened Fisher matrix of PM

  // Cumulative counter for rows written to table:
  long pointCount = 0;
  // Write vectors to table when this many rows accumulate:
  const long WriteChunk = 100000;

  // Write all matches to output catalog, deleting them along the way
  // and accumulating statistics of each exposure.
  //
  long matchCount = 0;
  for (auto im = matches.begin(); true; ++im, ++matchCount) {
    // First, write vectors to table if they've gotten big or we're at end
    if ( matchID.size() >= WriteChunk || im==matches.end()) {
      long nAdded = matchID.size();

      outTable.writeCells(matchID, "matchID", pointCount);
      matchID.clear();
      outTable.writeCells(catalogNumber, "extension", pointCount);
      catalogNumber.clear();
      outTable.writeCells(objectNumber, "object", pointCount);
      objectNumber.clear();
      outTable.writeCells(clip, "clip", pointCount);
      clip.clear();
      outTable.writeCells(reserve, "reserve", pointCount);
      reserve.clear();
      outTable.writeCells(xpix, "xPix", pointCount);
      xpix.clear();
      outTable.writeCells(ypix, "yPix", pointCount);
      ypix.clear();
      outTable.writeCells(color, "color", pointCount);
      color.clear();

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

      pointCount += nAdded;
    }	// Done flushing the vectors to Table

    if (im==matches.end()) break;

    const Match* m = *im;

    // Update all world coords and mean solutions
    m->remap(true);  // include remapping of clipped detections
    m->solve();

    // Don't report results for matches with no results
    if (m->getDOF() < 0) 
      continue;
    
    // Create a pointer to PMMatch if this is one:
    auto pmm = dynamic_cast<const astrometry::PMMatch*> (m); 

    astrometry::Vector2 xyMean;    // Match's mean position, in world coords
    if (!pmm)  {
      // We can use the same best-fit position for all detections
      xyMean = m->predict();
    }
    
    int detCount=0;
    int pmDetCount = 0;
    int clipCount = 0;
    
    double matchColor = astrometry::NODATA;

    for (auto detptr : *m) {

      // Get a pointer to a PMDetection if this is one:
      auto pmDetptr = dynamic_cast<const PMDetection*> (detptr);

      if (m->isFit(detptr)) {
	if (pmDetptr) {
	  pmDetCount++;
	} else {
	  detCount++;
	}
      } else {
	clipCount++;
      }
	
      // Set color if this is the first detection to have one
      if (matchColor!= astrometry::NODATA &&
	  detptr->color != astrometry::NODATA)
	matchColor = detptr->color;

      // Save some basics:
      matchID.push_back(matchCount);
      catalogNumber.push_back(detptr->catalogNumber);
      objectNumber.push_back(detptr->objectNumber);
      clip.push_back(detptr->isClipped);
      reserve.push_back(m->getReserved());
      hasPM.push_back( bool(pmDetptr));
      color.push_back(detptr->color);
      xpix.push_back(detptr->xpix);
      ypix.push_back(detptr->ypix);
      xworld.push_back(detptr->xw);
      yworld.push_back(detptr->yw);

      // Let's get the model prediction
      if (pmm)
	xyMean = pmm->predict(detptr);

      // Get world residuals (returned RESIDUAL_UNIT)
      astrometry::Vector2 residW = detptr->residWorld();
      xresw.push_back(residW[0]);
      yresw.push_back(residW[1]);

      // Get pixel residuals
      astrometry::Vector2 residP = detptr->residPix();
      xrespix.push_back(residP[0]);
      yrespix.push_back(residP[1]);

      astrometry::Matrix22 cov(0.);
      // Get error covariances, if we have any
      double chisqThis;
      if (detptr->invCov(0,0) > 0.) {
	cov = detptr->invCov.inverse() * pow(WCS_UNIT/RESIDUAL_UNIT,2.);
	vector<float> vv(3,0.);
	vv[0] = cov(0,0);
	vv[1] = cov(1,1);
	vv[2] = cov(0,1);
	covTotalW.push_back(vv);

	// Transform back to pixel errors to get circularized pixel error
	{
	  cov /= pow(WCS_UNIT/RESIDUAL_UNIT,2.); // Back to wcs units
	  auto dpdw = detptr->map->dPixdWorld(xyMean[0],xyMean[1]);
	  astrometry::Matrix22 covPix = dpdw * cov * dpdw.transpose();
	  double detCov = covPix(0,0)*covPix(1,1)-covPix(1,0)*covPix(0,1);
	  sigpix.push_back(detCov>0. ? pow(detCov,0.25) : 0.);
	}

	// Chisq and expected
	chisq.push_back(detptr->trueChisq());
	chisqExpected.push_back(detptr->expectedTrueChisq);
      } else {
	// Did not have a usable error matrix
	chisqThis = 0.;
	vector<float> vv(3,0.);
	covTotalW.push_back(vv);
	sigpix.push_back(0.);
	chisq.push_back(chisqThis);
	chisqExpected.push_back(detptr->expectedTrueChisq);
      }	

      if (pmDetptr) {
	// Add this object to PM output list
	pmMatchID.push_back(matchCount);
	pmCatalogNumber.push_back(detptr->catalogNumber);
	pmObjectNumber.push_back(detptr->objectNumber);
	pmClip.push_back(detptr->isClipped);
	pmReserve.push_back(m->getReserved());
	vector<float> vv(5);
	for (int i=0; i<5; i++)
	  vv[i] = pmDetptr->pmMean[i];
	pmMean.push_back(vv);
	vv.resize(25);
	int k=0;
	for (int i=0; i<5; i++)
	  for (int j=0; j<5; j++, k++)
	    vv[k] = pmDetptr->pmInvCov(i,j);
	pmInvCov.push_back(vv);
	pmChisq.push_back(detptr->trueChisq());
	pmChisqExpected.push_back(detptr->expectedTrueChisq);
      }
    } // End detection loop

    if (pmm) {
      // Add this object to the PM output catalog
      starMatchID.push_back(matchCount);
      starReserve.push_back(pmm->getReserved());
      starColor.push_back(matchColor);
      starPMCount.push_back(pmDetCount);
      starDetCount.push_back(detCount);
      starClipCount.push_back(clipCount);
      int dof=0.;
      double dummy;
      double chisqThis = pmm->chisq(dof,dummy);
      starDOF.push_back(dof);
      starChisq.push_back(chisqThis);
      {
	// Save the solution in a table
	auto pm = pmm->getPM();
	starX.push_back(pm[astrometry::X0]);
	starY.push_back(pm[astrometry::Y0]);
	starPMx.push_back(pm[astrometry::VX]);
	starPMy.push_back(pm[astrometry::VY]);
	starParallax.push_back(pm[astrometry::PAR]);
      }
      {
	// And the inverse covariance
	auto fisher = pmm->getInvCovPM();
	vector<float> vv(25);
	int k=0;
	for (int i=0; i<5; i++)
	  for (int j=0; j<5; j++, k++)
	    vv[k] = fisher(i,j);
	starInvCov.push_back(vv);
      }
    }  // Done adding a row to the star catalog.
  } // end match loop

  if (!pmMatchID.empty()) {
    // Create and fill the PMDetection output table if we have any
    FITS::FitsTable ft(outCatalog, FITS::ReadWrite + FITS::Create, pmTableName);
    auto pmDetTable = ft.use();

    pmDetTable.addColumn(pmMatchID,"matchID");
    pmDetTable.addColumn(pmCatalogNumber, "extension");
    pmDetTable.addColumn(pmObjectNumber, "object");
    pmDetTable.addColumn(pmClip, "clip");
    pmDetTable.addColumn(pmReserve, "reserve");
    pmDetTable.addColumn(pmMean, "pm");
    pmDetTable.addColumn(pmInvCov,"pmInvCov");
    pmDetTable.addColumn(pmChisq, "chisq");
    pmDetTable.addColumn(pmChisqExpected, "chisqExpected");
  }
    
  if (!starMatchID.empty()) {
    // Create and fill the star catalog table in the output file
    FITS::FitsTable ft(outCatalog, FITS::ReadWrite + FITS::Create, starTableName);
    auto starTable = ft.use();

    starTable.addColumn(starMatchID,"matchID");
    starTable.addColumn(starReserve, "reserve");
    starTable.addColumn(starPMCount, "pmCount");
    starTable.addColumn(starDetCount, "detCount");
    starTable.addColumn(starClipCount, "clipCount");
    starTable.addColumn(starChisq, "chisq");
    starTable.addColumn(starDOF, "dof");
    starTable.addColumn(starColor,"color");
    starTable.addColumn(starX, "xW");
    starTable.addColumn(starY, "yW");
    starTable.addColumn(starPMx, "pmX");
    starTable.addColumn(starPMy, "pmY");
    starTable.addColumn(starParallax, "parallax");
    starTable.addColumn(starInvCov,"pmInvCov");
  }
}

void
Photo::reportStatistics(const list<typename Photo::Match*>& matches,
			const vector<Exposure*>& exposures,
			const vector<typename Photo::Extension*>& extensions,
			ostream& os) {
  // Create Accum instances for fitted and reserved Detections on every
  // exposure, plus total accumulator for all reference
  // and all non-reference objects.
  typedef Accum<Photo> Acc;
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
    double magMean;
    double wtot;
    mptr->getMean(magMean, wtot);
    
    // Calculate number of DOF per Detection coordinate after
    // allowing for fit to centroid:
    int nFit = mptr->fitSize();
    double dofPerPt = (nFit > 1) ? 1. - 1./nFit : 0.;

    for (auto dptr : *mptr) {
      // Accumulate statistics for meaningful residuals
      int exposureNumber = extensions[dptr->catalogNumber]->exposure;
      Exposure* expo = exposures[exposureNumber];
      if (mptr->getReserved()) {
	if (expo->instrument==REF_INSTRUMENT ||
	    expo->instrument==PM_INSTRUMENT) {
	  // Treat PM like reference for photometry
	  refAccReserve.add(dptr, magMean, wtot, dofPerPt);
	  vaccReserve[exposureNumber].add(dptr, magMean, wtot, dofPerPt);
	} else if (expo->instrument==TAG_INSTRUMENT) {
	  // do nothing
	} else {
	  accReserve.add(dptr, magMean, wtot, dofPerPt);
	  vaccReserve[exposureNumber].add(dptr, magMean, wtot, dofPerPt);
	}
      } else {
	// Not a reserved object:
	if (expo->instrument==REF_INSTRUMENT ||
	    expo->instrument==PM_INSTRUMENT) {
	  // Treat PM like reference for photometry
	  refAccFit.add(dptr, magMean, wtot, dofPerPt);
	  vaccFit[exposureNumber].add(dptr, magMean, wtot, dofPerPt);
	} else if (expo->instrument==TAG_INSTRUMENT) {
	  // do nothing
	} else {
	  accFit.add(dptr, magMean, wtot, dofPerPt);
	  vaccFit[exposureNumber].add(dptr, magMean, wtot, dofPerPt);
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
  os << "#   Exposure           | " << Accum<Photo>::header() << endl;
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
  os << "Detection fit          | " << accFit.summary() << endl;
  if (accReserve.n>0)
    os << "Detection reserve      | " << accReserve.summary() << endl;
}

void
Astro::reportStatistics(const list<typename Astro::Match*>& matches,
			const vector<Exposure*>& exposures,
			const vector<typename Astro::Extension*>& extensions,
			ostream& os) {
  // Create Accum instances for fitted and reserved Detections on every
  // exposure, plus total accumulator for all reference
  // and all non-reference objects.
  typedef Accum<Astro> Acc;
  vector<Acc> vaccFit(exposures.size());
  vector<Acc> vaccReserve(exposures.size());
  Acc refAccFit;
  Acc refAccReserve;
  Acc accFit;
  Acc accReserve;

  // Accumulate stats over all detections.
  for (auto mptr : matches) {

    // Make sure match is ready, world coords for all detections done
    mptr->remap(true);
    mptr->solve();

    for (auto dptr : *mptr) {
      // Accumulate statistics for meaningful residuals
      int exposureNumber = extensions[dptr->catalogNumber]->exposure;
      Exposure* expo = exposures[exposureNumber];
      double magMean, wtot; // Dummy variables for accum::add() interface
      int dof;               // dummy also
      if (mptr->getReserved()) {
	if (expo->instrument==REF_INSTRUMENT ||
	    expo->instrument==PM_INSTRUMENT) {
	  refAccReserve.add(dptr, magMean, wtot, dof);
	  vaccReserve[exposureNumber].add(dptr, magMean, wtot, dof);
	} else if (expo->instrument==TAG_INSTRUMENT) {
	  // do nothing
	} else {
	  accReserve.add(dptr, magMean, wtot, dof);
	  vaccReserve[exposureNumber].add(dptr, magMean, wtot, dof);
	}
      } else {
	// Not a reserved object:
	if (expo->instrument==REF_INSTRUMENT ||
	    expo->instrument==PM_INSTRUMENT) {
	  refAccFit.add(dptr, magMean, wtot, dof);
	  vaccFit[exposureNumber].add(dptr, magMean, wtot, dof);
	} else if (expo->instrument==TAG_INSTRUMENT) {
	  // do nothing
	} else {
	  accFit.add(dptr, magMean, wtot, dof);
	  vaccFit[exposureNumber].add(dptr, magMean, wtot, dof);
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
  os << "#   Exposure           | " << Accum<Astro>::header() << endl;
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

  os << "Reference fit          | " << refAccFit.summary() << endl;
  if (refAccReserve.n>0)
    os << "Reference reserve      | " << refAccReserve.summary() << endl;

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
		    const vector<Exposure*>& exposures,		\
		    const vector<int>& exposureColorPriorities,	\
		    vector<AP::ColorExtension*>& colorExtensions,\
		    astrometry::YAMLCollector& inputYAML);       \
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
readObjects<AP>(const img::FTable& extensionTable, \
		const vector<Exposure*>& exposures, \
		vector<typename AP::Extension*>& extensions, \
		vector<astrometry::SphericalCoords*> fieldProjections); \
template void \
readMatches<AP>(img::FTable& table, \
		list<typename AP::Match*>& matches, \
		vector<typename AP::Extension*>& extensions, \
		vector<typename AP::ColorExtension*>& colorExtensions, \
		const ExtensionObjectSet& skipSet, \
		int minMatches, \
		bool usePM, \
		double parallaxPrior);  \
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
matchCensus<AP>(const list<AP::Match*>& matches, ostream& os);	\
template void  \
clipReserved<AP>(AP::Align& ca,  \
		 double clipThresh,  \
		 double minimumImprovement,  \
		 bool clipEntireMatch,  \
		 bool reportToCerr); \

INSTANTIATE(Astro)
INSTANTIATE(Photo)
