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

#include "FitSubroutines.h"
#include "MapDegeneracies.h"


using namespace std;
using namespace stringstuff;
using namespace photometry;
using img::FTable;

typedef Photo::ColorExtension ColorExtension;
typedef Photo::Extension Extension;


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

int
main(int argc, char *argv[])
{
  double reserveFraction;
  int randomNumberSeed;
  string skipFile;

  double clipThresh;
  double maxMagError;
  string sysErrorColumn;
  double sysError;

  int minMatches;
  int minFitExposures;
  bool clipEntirePrior;
  bool clipEntireMatch;
  double priorClipThresh;
  double chisqTolerance;

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

  Pset parameters;
  {
    const int def=PsetMember::hasDefault;
    const int low=PsetMember::hasLowerBound;
    const int up=PsetMember::hasUpperBound;
    const int lowopen = low | PsetMember::openLowerBound;
    const int upopen = up | PsetMember::openUpperBound;

    parameters.addMemberNoValue("INPUTS");
    parameters.addMember("maxMagError",&maxMagError, def | lowopen,
			 "cut objects with magnitude uncertainty above this", 0.05, 0.);
    parameters.addMember("sysErrorColumn",&sysErrorColumn, def,
			 "Extension table column holding systematic error, if any", "");
    parameters.addMember("sysError",&sysError, def | low,
			 "additional systematic error for all detected mags", 0.002, 0.);
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
    parameters.addMemberNoValue("OUTPUTS");
    parameters.addMember("outCatalog",&outCatalog, def,
			 "Output FITS binary catalog", "photo.fits");
    parameters.addMember("outPhotFile",&outPhotFile, def,
			 "Output serialized photometric solutions", "photfit.phot");
    parameters.addMember("outPriorFile",&outPriorFile, def,
			 "Output listing of zeropoints etc of exposures tied by priors","");
}

  // Fractional reduction in RMS required to continue sigma-clipping:
  const double minimumImprovement=0.02;

  try {
    
    // Read all the command-line and parameter-file program parameters
    processParameters(parameters, usage, 1, argc, argv);
    string inputTables = argv[1];

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

    // Let's figure out which of our FITS extensions are Instrument or MatchCatalog
    vector<int> instrumentHDUs;
    vector<int> catalogHDUs;
    inventoryFitsTables(inputTables, instrumentHDUs, catalogHDUs);
    
    /**/cerr << "Found " << instrumentHDUs.size() << " instrument HDUs" 
	     << " and " << catalogHDUs.size() << " catalog HDUs" << endl;

    // Read the fields table & propagate the info,
    // and discard info as it is not used here.
    {
      NameIndex fieldNames;
      vector<astrometry::SphericalCoords*> fieldProjections;
      readFields(inputTables, outCatalog, fieldNames, fieldProjections);
      for (auto ptr : fieldProjections) delete ptr;
    }

    // This flag will be set if we have already opened (and overwritten) the
    // output FITS catalog, as we have in readFields
    
    bool outputCatalogAlreadyOpen = true;

    // Read in all the instrument extensions and their device info from input
    // FITS file, save useful ones and write to output FITS file.
    /**/cerr << "Reading instruments" << endl;
    vector<Instrument*> instruments =
      readInstruments(instrumentHDUs, useInstrumentList, inputTables, outCatalog,
		      outputCatalogAlreadyOpen);
    

    // This vector will hold the color-priority value of each exposure.  -1 means an exposure
    // that does not hold color info.
    vector<int> exposureColorPriorities;
    // Read in the table of exposures
    /**/cerr << "Reading exposures" << endl;
    vector<Exposure*> exposures =
      readExposures(instruments,
		    exposureColorPriorities,
		    useColorList,
		    inputTables,
		    outCatalog,
		    skipExposureList,
		    false, // Do not use reference exposures for photometry
		    outputCatalogAlreadyOpen);


    // Read info about all Extensions - we will keep the Table around.
    /**/cerr << "Reading extensions" << endl;
    FTable extensionTable;
    {
      FITS::FitsTable ft(inputTables, FITS::ReadOnly, "Extensions");
      extensionTable = ft.extract();
      FITS::FitsTable out(outCatalog, FITS::ReadWrite+FITS::Create, "Extensions");
      out.copy(extensionTable);
    }

    vector<ColorExtension*> colorExtensions;
    vector<Extension*> extensions =
      readExtensions<Photo>(extensionTable,
			    instruments,
			    exposures,
			    exposureColorPriorities,
			    colorExtensions,
			    inputYAML,
			    sysErrorColumn,
			    sysError);


    /////////////////////////////////////////////////////
    //  Create and initialize all magnitude maps
    /////////////////////////////////////////////////////

    /**/cerr << "Building initial PhotoMapCollection" << endl;
    // Now build a preliminary set of photo maps from the configured YAML files
    PhotoMapCollection* pmcInit = new PhotoMapCollection;
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

    /**/cerr << "Checking for polynomial degeneracies" << endl;
    set<string> degenerateTypes={"Poly","Linear","Constant"};
    {
      MapDegeneracies<Photo> degen(extensions,
				   *pmcInit,
				   degenerateTypes,
				   false);  // Any non-fixed maps are examined
      // All exposure maps are candidates for setting to Identity
      // (the code will ignore those which already are Identity)
      set<string> exposureMapNames;
      for (auto expoPtr : exposures) {
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

    /**/cerr << "Making final mapCollection" << endl;
    
    // Do not need the preliminary PMC any more.
    delete pmcInit;
    pmcInit = 0;
    // And clean out any old maps stored in the YAMLCollector
    inputYAML.clearMaps();

    // Make final map collection from the YAML
    PhotoMapCollection mapCollection;
    createMapCollection<Photo>(instruments,
			       exposures,
			       extensions,
			       inputYAML,
			       mapCollection);

    /**/cerr << "Defining all maps" << endl;

    // Now construct a SubMap for every extension
    for (auto extnptr : extensions) {
      if (!extnptr) continue;  // not in use.
      Exposure& expo = *exposures[extnptr->exposure];
      if ( expo.instrument < 0) {
	// Reference exposures have no instruments, but possible color term per exposure
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
    
    /**/cerr << "Total number of free map elements " << mapCollection.nFreeMaps()
	     << " with " << mapCollection.nParams() << " free parameters."
	     << endl;


    //////////////////////////////////////////////////////////
    // Read in all the data
    //////////////////////////////////////////////////////////

    // List of all Matches - they will hold pointers to all Detections too.
    list<Match*> matches;    

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
      /**/cerr << "Parsing catalog field " << dummy1 << " Affinity " << dummy2 << endl;
      
      readMatches<Photo>(ff, matches, extensions, colorExtensions, skipSet, minMatches);
      
    } // End loop over input matched catalogs

    /**/cerr << "Total match count: " << matches.size() << endl;

    // Now loop over all original catalog bintables, reading the desired rows
    // and collecting needed information into the Detection structures
    readObjects<Photo>(extensionTable, exposures, extensions);
    
    /**/cerr << "Done reading catalogs for magnitudes." << endl;

    // Now loop again over all catalogs being used to supply colors,
    // and insert colors into all the PhotoArguments for Detections they match
    readColors<Photo>(extensionTable, colorExtensions);

    cerr << "Done reading catalogs for colors." << endl;

    /**/cerr << "Done reading catalogs." << endl;

    // Get rid of Detections with errors too high or already clipped
    purgeNoisyDetections<Photo>(maxMagError, matches, exposures, extensions);
			 
    // Get rid of Matches with too few detections
    purgeSparseMatches<Photo>(minMatches, matches);

    // Get rid of Matches with color out of range (NODATA are kept)
    purgeBadColor<Photo>(minColor, maxColor, matches); 
    
    /**/cerr << "Done purging defective detections and matches" << endl;

    // Reserve desired fraction of matches
    if (reserveFraction>0.) 
      reserveMatches<Photo>(matches, reserveFraction, randomNumberSeed);

    // Find exposures whose parameters are free but have too few
    // Detections being fit to the exposure model.
    auto badExposures = findUnderpopulatedExposures<Photo>(minFitExposures,
							   matches,
							   exposures,
							   extensions,
							   mapCollection);

    // Freeze parameters of an exposure model and clip all
    // Detections that were going to use it.
    for (auto i : badExposures) {
      cout << "WARNING: Shutting down exposure map " << i.first
	   << " with only " << i.second
	   << " fitted detections "
	   << endl;
      freezeMap<Photo>(i.first, matches, extensions, mapCollection);
    } 

    matchCensus<Photo>(matches, cout);

    ///////////////////////////////////////////////////////////
    // Construct priors
    ///////////////////////////////////////////////////////////

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
	cout << "Fitting " << mcount << " matches with " << dcount << " detections "
	     << " chisq " << chi << " / " << dof << " dof,  maxdev " << maxdev 
	     << " sigma" << endl;
      }

      // Do the fit here!!
      double chisq = ca.fitOnce();
      // Note that fitOnce() remaps *all* the matches, including reserved ones.
      double max;
      int dof;
      ca.chisqDOF(dof, max, false);	// Exclude reserved Matches
      double thresh = sqrt(chisq/dof) * clipThresh;
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
	  do {
	    if (nclip > 0) {
	      // Refit the data after clipping priors
	      ca.fitOnce();
	      cout << "After prior clip: chisq " << chisq 
		   << " / " << dof << " dof, max deviation " << max
		   << endl;
	    }
	    // Clip up to one exposure per prior
	    nclip = ca.sigmaClipPrior(priorClipThresh, false);
	    cout << "Clipped " << nclip << " prior reference points" << endl;
	  } while (nclip > 0);
	  ca.setRelTolerance(chisqTolerance);
	  /**/cerr << "--Starting strict tolerance passes";
	  /**/if (clipEntireMatch) cerr << "; clipping full matches";
	  /**/cerr << endl;
	  oldthresh = thresh;
	  nclip = ca.sigmaClip(thresh, false, true);
	  cout << "Clipped " << nclip
	       << " matches " << endl;
	  continue;
	} else {
	  // Done!
	  break;
	}
      }
      oldthresh = thresh;
      nclip = ca.sigmaClip(thresh, false, clipEntireMatch && !coarsePasses);
      if (nclip==0 && coarsePasses) {
	// Nothing being clipped; tighten tolerances and re-fit
	coarsePasses = false;
	ca.setRelTolerance(chisqTolerance);
	/**/cerr << "--Starting strict tolerance passes";
	/**/if (clipEntireMatch) cerr << "; clipping full matches";
	/**/cerr << endl;
	continue;
      }
      cout << "Clipped " << nclip
	   << " matches " << endl;
      long int dcount=0;
      long int mcount=0;
      
    } while (coarsePasses || nclip>0);
  
    // The re-fitting is now complete.  Serialize all the fitted magnitude solutions
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
      cout << "** Clipping reserved matches: " << endl;
      //turn on cerr logging, no entire-match clipping
      clipReserved<Photo>(ca, clipThresh, minimumImprovement,
			  false, true);  
    }

    //////////////////////////////////////
    // Output report on priors
    //////////////////////////////////////
    if (!outPriorFile.empty()) {
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
      /**/cerr << "Wrote prior results to " << outPriorFile << endl;
    }

    //////////////////////////////////////
    // Output data and calculate some statistics
    //////////////////////////////////////

    // Save the pointwise fitting results
    saveResults<Photo>(matches, outCatalog);
    
    /**/cerr << "Saved results to FITS table" << endl;
    
    // Report summary of residuals to stdout
    reportStatistics<Photo>(matches, exposures, extensions, cout);

    //////////////////////////////////////
    // Cleanup:
    //////////////////////////////////////

    // Get rid of matches:
    for (auto im = matches.begin(); im!=matches.end(); ) {
      (*im)->clear(true);  // deletes detections
      // And get rid of match itself.
      im = matches.erase(im);
    }
    // Get rid of extensions
    for (int i=0; i<extensions.size(); i++)
      if (extensions[i]) delete extensions[i];
    // Get rid of exposures
    for (int i=0; i<exposures.size(); i++)
      if (exposures[i]) delete exposures[i];
    // Get rid of instruments
    for (int i=0; i<instruments.size(); i++)
      if (instruments[i]) delete instruments[i];

  } catch (std::runtime_error& m) {
    quit(m,1);
  }
}

