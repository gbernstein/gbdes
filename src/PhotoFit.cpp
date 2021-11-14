// Program to fit photometric solutions to detection catalogs already matched by WCSFoF.
#include <fstream>
#include <sstream>
#include <map>
#include <algorithm>

#include "Std.h"
#include "Astrometry.h"
#include "FitsTable.h"
#include "StringStuff.h"
#include "Pset.h"
#include "PixelMapCollection.h"
#include "YAMLCollector.h"
#include "PhotoMatch.h"
#include "PhotoMapCollection.h"
#include "Instrument.h"
#include "Units.h"

#include "FitSubroutines.h"
#include "MapDegeneracies.h"


using namespace std;
using namespace stringstuff;
using namespace photometry;
using img::FTable;

typedef Photo::ColorExtension ColorExtension;
typedef Photo::Extension Extension;

#define PROGRESS(val, msg) if (verbose>=val) cerr << "-->" <<  #msg << endl

string usage=
  "PhotoFit: Refine photometric solutions for a matched set of catalogs.\n"
  "usage: PhotoFit <match file> [parameter file] [parameter file...]\n"
  "   [-parameter[=]value...]\n"
  "      <match file>:  FITS file with binary tables produced by WCSFoF\n"
  "      Program parameters specified as command-line options or read from\n"
  "          parameter file(s) specified on cmd line";

// Temporary documentation:
// parameter renameInstruments is a string that has format like
// <instrument name>=<map name>, <instrument name> = <map name>, ...
// where the Instrument's PixelMap will be given the map name.  Regex allowed.
// Note that this is assuming that regexes do not include = or , characters.
// whitespace will be stripped from edges of each name.

// parameter existingMaps is a string with
// <mapName>@<file name>, ...
// which says that maps with the mapName should be initially be deserialized from file.
// Regexes ok; same caveats: no @ or , in names, whitespace stripped.

// parameter fixMaps is a string with
// <mapName>, ...
// where any given mapName should have its parameters fixed at initial values during
// the fitting.  Regexes allowed (no commas!).

// parameter canonicalExposures is a string with
// <exposureID>, ....
// where ares are exposures that will be given an identity exposure map.  There must be
// 0 or 1 of these specified for any instrument that has Instrument Map with free parameters.
// Default is to find an exposure that has data in all devices and use it.
// Will have an error if there is more than one constraint on any Instrument.

// skipFile parameter gives name of file holding extension and object numbers of detections to ignore.

// Note that PhotoMaps for devices within instrument will get names <instrument>/<device>.
// And PhotoMaps for individual exposures will get names <exposure>/<device>.

// magKey and magErrKey can be of form COLUMN[#] when COLUMN is an array column (float or double)
// and # is the element of this array that we want to use.

// The sysError command-line option will override values of photSysVar that are in input Exposure table,
// if it is >0.

int
main(int argc, char *argv[])
{
  double reserveFraction;
  int randomNumberSeed;
  string skipFile;

  double clipThresh;
  double maxMagError;
  double sysError;

  int minMatches;
  int minFitExposures;
  bool clipEntirePrior;
  bool clipEntireMatch;
  double priorClipThresh;
  double chisqTolerance;
  bool divideInPlace;
  bool purgeOutput;

  string inputMaps;
  string fixMaps;
  string priorFiles;
  string useInstruments;
  string skipExposures;

  string outCatalog;
  string outPhotFile;
  string outPriorFile;

  string colorExposures;
  double minColor;
  double maxColor;

  int verbose;
  
  Pset parameters;
  {
    const int def=PsetMember::hasDefault;
    const int low=PsetMember::hasLowerBound;
    const int up=PsetMember::hasUpperBound;
    const int lowopen = low | PsetMember::openLowerBound;
    const int upopen = up | PsetMember::openUpperBound;

    parameters.addMemberNoValue("INPUTS");
    parameters.addMember("maxMagError",&maxMagError, def | lowopen,
			 "cut objects with magnitude uncertainty above this (mmag)",
			 0.05, 0.);
    parameters.addMember("sysError",&sysError, def | low,
			 "systematic error to add to exposures (mmag)", 2., 0.);
    parameters.addMember("minMatch",&minMatches, def | low,
			 "minimum number of detections for usable match", 2, 2);
    parameters.addMember("minFitExposures",&minFitExposures, def | low,
			 "Minimum number of detections to fit exposure map", 200, 2);
    parameters.addMember("useInstruments",&useInstruments, def,
			 "the instruments to include in fit",".*");
    parameters.addMember("skipExposures",&skipExposures, def,
			 "exposures to ignore during fitting","");
    parameters.addMemberNoValue("COLORS");
    parameters.addMember("colorExposures",&colorExposures, def,
			 "exposures holding valid colors for stars","");
    parameters.addMember("minColor",&minColor, def,
			 "minimum value of color to be used",-10.);
    parameters.addMember("maxColor",&maxColor, def,
			 "maximum value of color to be used",+10.);
    parameters.addMemberNoValue("CLIPPING");
    parameters.addMember("clipThresh",&clipThresh, def | low,
			 "Clipping threshold (sigma)", 5., 2.);
    parameters.addMember("clipEntireMatch",&clipEntireMatch, def,
			 "Discard entire object if one outlier on later passes", false);
    parameters.addMember("priorClipThresh",&priorClipThresh, def | low,
			 "Clipping threshold (sigma)", 5., 2.);
    parameters.addMember("clipEntirePrior",&clipEntirePrior, def,
			 "Discard entire night's prior if one outlier", false);
    parameters.addMember("skipFile",&skipFile, def,
			 "optional file holding extension/object of detections to ignore","");
    parameters.addMemberNoValue("FITTING");
    parameters.addMember("priorFiles", &priorFiles, def,
			 "File(s) specifying any priors to apply to zeropoints and colors", "");
    parameters.addMember("reserveFraction",&reserveFraction, def | low,
			 "Fraction of matches reserved from re-fit", 0., 0.);
    parameters.addMember("seed",&randomNumberSeed, def,
			 "seed for reserving randomizer, <=0 to seed with time", 0);
    parameters.addMember("chisqTolerance",&chisqTolerance, def | lowopen,
			 "Fractional change in chisq for convergence", 0.001, 0.);
    parameters.addMember("inputMaps",&inputMaps, def,
			 "list of YAML files specifying maps","");
    parameters.addMember("fixMaps",&fixMaps, def,
			 "list of map components or instruments to hold fixed","");
    parameters.addMember("divideInPlace",&divideInPlace, def,
			 "use in-place Cholesky to save memory but lose debug of degeneracies",false);
    parameters.addMemberNoValue("OUTPUTS");
    parameters.addMember("purgeOutput",&purgeOutput, def,
			 "Purge un-fittable maps from output", false);
    parameters.addMember("outCatalog",&outCatalog, def,
			 "Output FITS binary catalog", "photo.fits");
    parameters.addMember("outPhotFile",&outPhotFile, def,
			 "Output serialized photometric solutions", "photfit.phot");
    parameters.addMember("outPriorFile",&outPriorFile, def,
			 "Output listing of zeropoints etc of exposures tied by priors","");
    parameters.addMember("verbose", &verbose, def,
			 "stderr detail level", 1);
}

  // Fractional reduction in RMS required to continue sigma-clipping:
  const double minimumImprovement=0.02;

  try {
    
    // Read all the command-line and parameter-file program parameters
    processParameters(parameters, usage, 1, argc, argv);
    string inputTables = argv[1];

    sysError *= MMAG;
    maxMagError *= MMAG;
    
    /////////////////////////////////////////////////////
    // Parse all the parameters 
    /////////////////////////////////////////////////////

    
    // Teach PhotoMapCollection about new kinds of PhotoMaps and PixelMaps
    loadPhotoMapParser();
    loadPixelMapParser();


    // This is list of regexes of PhotoMap names (or instrument names) that should
    // have their parameters held fixed.
    list<string> fixMapList = splitArgument(fixMaps);

    // The list of instruments that we will be matching together in this run:
    list<string> useInstrumentList = splitArgument(useInstruments);

    // The list of exposures that are considered valid sources of color information:
    list<string> useColorList = splitArgument(colorExposures);

    // Objects to ignore on input:
    ExtensionObjectSet skipSet(skipFile);

    // Exposures to ignore:
    list<string> skipExposureList = splitArgument(skipExposures);
    
    // Class that will build a starting YAML config for all extensions
    astrometry::YAMLCollector inputYAML(inputMaps, PhotoMapCollection::magicKey);
    // Make sure inputYAML knows about the Identity transformation:
    {
      istringstream iss("Identity:\n  Type:  Identity\n");
      inputYAML.addInput(iss);
    }
    /////////////////////////////////////////////////////
    //  Read in properties of all Instruments, Devices, Exposures
    /////////////////////////////////////////////////////

    // All the names will be stripped of leading/trailing white space, and internal
    // white space replaced with single underscore - this keeps PhotoMap parsing functional.


    PROGRESS(1,Reading fields);

    vector<double> fieldEpochs;
    vector<std::unique_ptr<astrometry::SphericalCoords>> fieldProjections;
    // Read the fields table & propagate the info,
    // and discard info as it is not used here.
    {
      NameIndex fieldNames;
      readFields(inputTables, outCatalog, fieldNames, fieldProjections, fieldEpochs);
      fieldProjections.clear();
    }

    // This flag will be set if we have already opened (and overwritten) the
    // output FITS catalog, as we have in readFields
    
    bool outputCatalogAlreadyOpen = true;

    // Let's figure out which of our FITS extensions are Instrument or MatchCatalog
    vector<int> instrumentHDUs;
    vector<int> catalogHDUs;
    inventoryFitsTables(inputTables, instrumentHDUs, catalogHDUs);
    
    if (verbose>=1)
      cerr << "Found " << instrumentHDUs.size() << " instrument HDUs" 
	   << " and " << catalogHDUs.size() << " catalog HDUs" << endl;


    PROGRESS(1,Reading instruments);

    // Read in all the instrument extensions and their device info from input
    // FITS file, save useful ones and write to output FITS file.
    vector<unique_ptr<Instrument>> instruments =
      readInstruments(instrumentHDUs, useInstrumentList, inputTables, outCatalog,
		      outputCatalogAlreadyOpen);
    
    PROGRESS(1,Reading exposures);

    // This vector will hold the color-priority value of each exposure.
    // -1 means an exposure that does not hold color info.
    vector<int> exposureColorPriorities;
    // Read in the table of exposures
    vector<unique_ptr<Exposure>> exposures =
      readExposures(instruments,
		    fieldEpochs,
		    exposureColorPriorities,
		    useColorList,
		    inputTables,
		    outCatalog,
		    skipExposureList,
		    false, // Do not use reference exposures for photometry
		    outputCatalogAlreadyOpen);

    if (sysError > 0.) {
      // Add global systematic error to exposuress
      double photometricVariance = sysError*sysError;
      for (auto const & e : exposures) {
	e->photometricVariance += photometricVariance;
      }
    }

    PROGRESS(1,Reading extensions);

    // Read info about all Extensions - we will keep the Table around.
    FTable extensionTable;
    {
      FITS::FitsTable ft(inputTables, FITS::ReadOnly, "Extensions");
      extensionTable = ft.extract();
      FITS::FitsTable out(outCatalog, FITS::ReadWrite+FITS::Create, "Extensions");
      out.copy(extensionTable);
    }

    vector<unique_ptr<ColorExtension>> colorExtensions;
    vector<unique_ptr<Extension>> extensions =
      readExtensions<Photo>(extensionTable,
			    instruments,
			    exposures,
			    exposureColorPriorities,
			    colorExtensions,
			    inputYAML);


    /////////////////////////////////////////////////////
    //  Create and initialize all magnitude maps
    /////////////////////////////////////////////////////

    PROGRESS(1, Building initial PhotoMapCollection);
    // Now build a preliminary set of photo maps from the configured YAML files
    unique_ptr<PhotoMapCollection> pmcInit(new PhotoMapCollection);
    pmcInit->learnMap(IdentityMap(), false, false);
    {
      istringstream iss(inputYAML.dump());
      if (!pmcInit->read(iss)) {
	cerr << "Failure parsing the initial YAML map specs" << endl;
	exit(1);
      }
    }

    
    /////////////////////////////////////////////////////
    // Check for degeneracies between compounded linear/poly/constant
    // maps. See if it is needed and sufficient to set some exposure
    // maps to Identity to break such degeneracies.
    /////////////////////////////////////////////////////

    // Now fix all map elements requested to be fixed
    fixMapComponents<Photo>(*pmcInit,
			    fixMapList,
			    instruments);

    PROGRESS(2,Checking for polynomial degeneracies);
    set<string> degenerateTypes={"Poly","Linear","Constant"};
    {
      MapDegeneracies<Photo> degen(extensions,
				   *pmcInit,
				   degenerateTypes,
				   false);  // Any non-fixed maps are examined
      // All exposure maps are candidates for setting to Identity
      // (the code will ignore those which already are Identity)
      set<string> exposureMapNames;
      for (auto const & expoPtr : exposures) {
	if (expoPtr && !expoPtr->name.empty())
	  exposureMapNames.insert(expoPtr->name);
      }

      auto replaceThese = degen.replaceWithIdentity(exposureMapNames);
      
      // Supersede their maps if there are any
      if (!replaceThese.empty()) {
	PhotoMapCollection pmcAltered;
	for (auto mapname : replaceThese) {
	  cerr << "..Setting map <" << mapname << "> to Identity" << endl;
	  pmcAltered.learnMap(IdentityMap(mapname));
	}	  
	// Add the altered maps specs to the input YAML specifications
	ostringstream oss;
	pmcAltered.write(oss);
	istringstream iss(oss.str());
	inputYAML.addInput(iss, "", true); // Prepend these specs to others
      }
    } // End of poly-degeneracy check/correction

    PROGRESS(1,Making final mapCollection);

    // Do not need the preliminary PMC any more.
    pmcInit.reset();
    // And clean out any old maps stored in the YAMLCollector
    inputYAML.clearMaps();

    // Make final map collection from the YAML
    PhotoMapCollection mapCollection;
    createMapCollection<Photo>(instruments,
			       exposures,
			       extensions,
			       inputYAML,
			       mapCollection);

    PROGRESS(2, Defining all maps);

    // Now construct a SubMap for every extension
    for (auto const & extnptr : extensions) {
      if (!extnptr) continue;  // not in use.
      Exposure& expo = *exposures[extnptr->exposure];
      if ( expo.instrument < 0) {
	// Reference exposures have no instruments, shouldn't be here
	cerr << "Trying to make PhotoMap for reference exposure "
	     << expo.name
	     << endl;
	exit(1);
      } 
      // Real instrument, make a map combining its exposure with its Device map:
      extnptr->map = mapCollection.issueMap(extnptr->mapName);
    } // Extension loop

    // Now fix all map elements requested to be fixed
    fixMapComponents<Photo>(mapCollection,
			    fixMapList,
			    instruments);
    
    // Recalculate all parameter indices - maps are ready to roll!
    mapCollection.rebuildParameterVector();
    
    cout << "Total number of free map elements " << mapCollection.nFreeMaps()
	 << " with " << mapCollection.nParams() << " free parameters."
	 << endl;


    //////////////////////////////////////////////////////////
    // Read in all the data
    //////////////////////////////////////////////////////////

    // List of all Matches - they will hold pointers to all Detections too.
    list<unique_ptr<Match>> matches;

    // Figure out which extensions' maps require a color entry to function
    whoNeedsColor<Photo>(extensions);
    
    // Start by reading all matched catalogs, creating Detection and Match arrays, and 
    // telling each Extension which objects it should retrieve from its catalog

    for (int icat = 0; icat < catalogHDUs.size(); icat++) {
      FITS::FitsTable ft(inputTables, FITS::ReadOnly, catalogHDUs[icat]);
      FTable ff = ft.use();
      {
	// Only do photometric fitting on the stellar objects:
	string affinity;
	if (!ff.getHdrValue("Affinity", affinity)) {
	  cerr << "Could not find affinity keyword in header of extension " 
	       << catalogHDUs[icat] 
	       << endl;
	  exit(1);
	}
	stripWhite(affinity);
	if (!stringstuff::nocaseEqual(affinity, stellarAffinity))
	  continue;
      }
      string dummy1, dummy2;
      ff.getHdrValue("Field", dummy1);
      ff.getHdrValue("Affinity", dummy2);
      if (verbose>=2)
	cerr << "-->Parsing catalog field " << dummy1 << " Affinity " << dummy2 << endl;
      
      readMatches<Photo>(ff, matches, extensions, colorExtensions, skipSet, minMatches);
      
    } // End loop over input matched catalogs

    if (verbose>=0) cerr << "Total match count: " << matches.size() << endl;

    // Now loop over all original catalog bintables, reading the desired rows
    // and collecting needed information into the Detection structures
    PROGRESS(1,Reading catalogs);
    readObjects<Photo>(extensionTable, exposures, extensions, fieldProjections);
    
    // Now loop again over all catalogs being used to supply colors,
    // and insert colors into all the PhotoArguments for Detections they match
    PROGRESS(1,Reading colors);
    readColors<Photo>(extensionTable, colorExtensions);

    PROGRESS(2,Purging defective detections and matches);
    // Get rid of Detections with errors too high or already clipped
    purgeNoisyDetections<Photo>(maxMagError, matches, exposures, extensions);
			 
    PROGRESS(2,Purging sparse matches);
    // Get rid of Matches with too few detections
    purgeSparseMatches<Photo>(minMatches, matches);

    PROGRESS(2,Purging out-of-range colors);
    // Get rid of Matches with color out of range (NODATA are kept)
    purgeBadColor<Photo>(minColor, maxColor, matches); 
    
    PROGRESS(2, Reserving matches);
    // Reserve desired fraction of matches
    if (reserveFraction>0.) 
      reserveMatches<Photo>(matches, reserveFraction, randomNumberSeed);

    PROGRESS(2,Purging underpopulated exposures);
    // Find exposures whose parameters are free but have too few
    // Detections being fit to the exposure model.
    auto badExposures = findUnderpopulatedExposures<Photo>(minFitExposures,
							   matches,
							   exposures,
							   extensions,
							   mapCollection);

    PROGRESS(2,Purging bad exposure parameters and Detections);
    // Freeze parameters of an exposure model and clip all
    // Detections that were going to use it.
    for (auto i : badExposures) {
      cout << "WARNING: Shutting down exposure map " << i.first
	   << " with only " << i.second
	   << " fitted detections "
	   << endl;
      freezeMap<Photo>(i.first, matches, extensions, mapCollection);
    } 

    if (purgeOutput) {
      PROGRESS(2,Purging unfittable maps);
      mapCollection.purgeInvalid();
    }

    PROGRESS(2,Match census);
    matchCensus<Photo>(matches, cout);

    ///////////////////////////////////////////////////////////
    // Construct priors
    ///////////////////////////////////////////////////////////

    PROGRESS(2,Constructing priors);
    list<PhotoPrior*> priors;
    {
      // List of the files to use 
      list<string> priorFileList = splitArgument(priorFiles);
      // Call subroutine that will read & initialize priors from a file
      for (list<string>::const_iterator i = priorFileList.begin();
	   i != priorFileList.end();
	   ++i) {
	list<PhotoPrior*> thisTime = readPriors(*i, instruments, exposures, extensions);
	//???						detectionsPerExposure);
	priors.splice(priors.end(), thisTime);
      }
    }

    // Remove any priors that have insufficient data to fit
    for (list<PhotoPrior*>::iterator i = priors.begin();
	 i != priors.end(); ) {
      if ((*i)->isDegenerate()) {
	cout << "WARNING: PhotoPrior " << (*i)->getName()
	     << " is degenerate and is being deleted " << endl;
	i = priors.erase(i);
      } else {
	++i;
      }
    }

    ///////////////////////////////////////////////////////////
    // Now do the re-fitting 
    ///////////////////////////////////////////////////////////

    PROGRESS(1,Begin fitting process);
    
    // make CoordAlign class
    PhotoAlign ca(mapCollection, matches, priors);

    int nclip;
    double oldthresh=0.;

    // Start off in a "coarse" mode so we are not fine-tuning the solution
    // until most of the outliers have been rejected:
    bool coarsePasses = true;
    ca.setRelTolerance(10.*chisqTolerance);
    // Here is the actual fitting loop 
    do {
      // Report number of active Matches / Detections in each iteration:
      {
	long int mcount=0;
	long int dcount=0;
	ca.count(mcount, dcount, false, 2);
	double maxdev=0.;
	int dof=0;
	double chi= ca.chisqDOF(dof, maxdev, false);
	if (verbose>=1)
	  cout << "Fitting " << mcount << " matches with " << dcount << " detections "
	       << " chisq " << chi << " / " << dof << " dof,  maxdev " << maxdev 
	       << " sigma" << endl;
      }

      // Do the fit here!!
      double chisq = ca.fitOnce(verbose>=1, divideInPlace);
      // Note that fitOnce() remaps *all* the matches, including reserved ones.
      double max;
      int dof;
      ca.chisqDOF(dof, max, false);	// Exclude reserved Matches
      double thresh = sqrt(chisq/dof) * clipThresh;
      if (verbose>=1)
	cout << "After fit: chisq " << chisq 
	     << " / " << dof << " dof, max deviation " << max
	     << "  new clip threshold at: " << thresh << " sigma"
	     << endl;
      if (thresh >= max || (oldthresh>0. && (1-thresh/oldthresh)<minimumImprovement)) {
	// Sigma clipping is no longer doing much.  Quit if we are at full precision,
	// else require full target precision and initiate another pass.
	if (coarsePasses) {
	  coarsePasses = false;
	  // This is the point at which we will clip aberrant exposures from the priors,
	  // after coarse passes are done, before we do fine passes.
	  nclip = 0;
	  PROGRESS(1,Starting prior clipping);
	  do {
	    if (nclip > 0) {
	      // Refit the data after clipping priors
	      ca.fitOnce(verbose>=1, divideInPlace);
	      if (verbose>=1)
		cout << "After prior clip: chisq " << chisq 
		     << " / " << dof << " dof, max deviation " << max
		     << endl;
	    }
	    // Clip up to one exposure per prior
	    nclip = ca.sigmaClipPrior(priorClipThresh, false);
	    cout << "Clipped " << nclip << " prior reference points" << endl;
	  } while (nclip > 0);
	  ca.setRelTolerance(chisqTolerance);
	  PROGRESS(1,Starting strict tolerance passes);
	  oldthresh = thresh;
	  nclip = ca.sigmaClip(thresh, false, true);
	  if (verbose>=0)
	    cout << "Clipped " << nclip
		 << " matches " << endl;
	  continue;
	} else {
	  // Done!
	  break;
	}
      }
      oldthresh = thresh;
      nclip = ca.sigmaClip(thresh, false, clipEntireMatch && !coarsePasses,
			   verbose>=1);
      if (nclip==0 && coarsePasses) {
	// Nothing being clipped; tighten tolerances and re-fit
	coarsePasses = false;
	ca.setRelTolerance(chisqTolerance);
	PROGRESS(1,Starting strict tolerance passes);
	continue;
      }
      if (verbose>=0) 
	cout << "# Clipped " << nclip
	     << " matches " << endl;
      
    } while (coarsePasses || nclip>0);
  
    // The re-fitting is now complete.  Serialize all the fitted magnitude solutions
    PROGRESS(1,Saving solution);
    {
      ofstream ofs(outPhotFile.c_str());
      if (!ofs) {
	cerr << "Error trying to open output file for fitted photometry: "
	     << outPhotFile << endl;
	// *** will not quit before dumping output ***
      } else {
	mapCollection.write(ofs, stringstuff::taggedCommandLine(argc,argv));
      }
    }

    // If there are reserved Matches, run sigma-clipping on them now.
    if (reserveFraction > 0.) {
      PROGRESS(1,Clipping reserved matches);
      //turn on cerr logging, no entire-match clipping
      clipReserved<Photo>(ca, clipThresh, minimumImprovement,
			  false, verbose>=1);  
    }

    //////////////////////////////////////
    // Output report on priors
    //////////////////////////////////////
    if (!outPriorFile.empty()) {
      PROGRESS(1,Saving prior results);
      ofstream ofs(outPriorFile.c_str());
      if (!ofs) {
	cerr << "Error trying to open output file for prior report: " << outPriorFile << endl;
	// *** will not quit before dumping output ***
      } else {
	PhotoPrior::reportHeader(ofs);
	for (list<PhotoPrior*>::const_iterator i = priors.begin();
	     i != priors.end();
	     ++i) {
	  (*i)->report(ofs);
	}
      }
    }

    //////////////////////////////////////
    // Output data and calculate some statistics
    //////////////////////////////////////

    PROGRESS(1,Saving photometric residuals);
    // Save the pointwise fitting results
    Photo::saveResults(matches, outCatalog);
    
    // Report summary of residuals to stdout
    Photo::reportStatistics(matches, exposures, extensions, cout);

  } catch (std::runtime_error& m) {
    quit(m,1);
  }
}

