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
#include "Random.h"
#include "YAMLCollector.h"
#include "PhotoMatch.h"
#include "PhotoMapCollection.h"
#include "Instrument.h"

#include "FitSubroutines.h"
#include "PhotoSubs.h"

using namespace std;
using namespace stringstuff;
using namespace photometry;
using img::FTable;


string usage=
  "PhotoFit: Refine photometric solutions for a matched set of catalogs.\n"
  "usage: WCSFit <match file> [parameter file] [parameter file...]\n"
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
  double magSysError;

  int minMatches;
  bool clipEntirePrior;
  bool clipEntireMatch;
  double priorClipThresh;
  double chisqTolerance;

  string inputMaps;
  string fixMaps;
  string priorFiles;
  string useInstruments;

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
    parameters.addMember("magSysError",&magSysError, def | low,
			 "additional systematic error for all detected mags", 0.002, 0.);
    parameters.addMember("minMatch",&minMatches, def | low,
			 "minimum number of detections for usable match", 2, 2);
    parameters.addMember("useInstruments",&useInstruments, def,
			 "the instruments to include in fit",".*");
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
			 "Discard entire object if one outlier", false);
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
    // have their parameters held fixed.
    list<string> useInstrumentList = splitArgument(useInstruments);

    // The list of exposures that are considered valid sources of color information:
    list<string> useColorList = splitArgument(colorExposures);

    // Objects to ignore on input:
    ExtensionObjectSet skipSet(skipFile);

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

    // This flag will be set if we have already opened (and overwritten) the
    // output FITS catalog.
    bool outputCatalogAlreadyOpen = false;

    // Read in all the instrument extensions and their device info from input
    // FITS file, save useful ones and write to output FITS file.
    vector<Instrument*> instruments =
      readInstruments(instrumentHDUs, useInstrumentList, inputTables, outCatalog,
		      outputCatalogAlreadyOpen);
    

    // This vector will hold the color-priority value of each exposure.  -1 means an exposure
    // that does not hold color info.
    vector<int> exposureColorPriorities;
    // Read in the table of exposures
    vector<Exposure*> exposures =
      readExposures(instruments,
		    exposureColorPriorities,
		    useColorList,
		    inputTables,
		    outCatalog,
		    false, // Do not use reference exposures for photometry
		    outputCatalogAlreadyOpen);

    // Read info about all Extensions - we will keep the Table around.
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
			    inputYAML);

    /////////////////////////////////////////////////////
    //  Create and initialize all magnitude maps
    /////////////////////////////////////////////////////

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

    /**/cerr << "Successfully built initial PixelMapCollection" << endl;
    
    /////////////////////////////////////////////////////
    // Now for each instrument, determine if we need a canonical
    // exposure to break an exposure/instrument degeneracy.  If
    // so, choose one, and replace its exposure map with Identity.
    /////////////////////////////////////////////////////

    // Here's where we'll collect the map definitions that override the
    // inputs in order to eliminate degeneracies.
    PhotoMapCollection pmcAltered;  

    for (int iInst=0; iInst < instruments.size(); iInst++) {
      if (!instruments[iInst]) continue;  // Not using instrument
      auto& instr = *instruments[iInst];

      int canonicalExposure =
	findCanonical<Photo>(instr, iInst, exposures, extensions, *pmcInit);

      if (canonicalExposure >= 0)
	// Make a new map spec for the canonical exposure
	pmcAltered.learnMap(IdentityMap(exposures[canonicalExposure]->name));
    } // End instrument loop

    // Add the altered maps specs to the input YAML specifications
    {
      ostringstream oss;
      pmcAltered.write(oss);
      istringstream iss(oss.str());
      inputYAML.addInput(iss, "", true); // Prepend these specs to others
    }
    
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

    /**/cerr << "Done making final mapCollection!" << endl;

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
    readObjects<Photo>(extensionTable, exposures, extensions, magSysError);
    
    /**/cerr << "Done reading catalogs for magnitudes." << endl;

    // Now loop again over all catalogs being used to supply colors,
    // and insert colors into all the PhotoArguments for Detections they match
    for (int iext = 0; iext < colorExtensions.size(); iext++) {
      if (!colorExtensions[iext]) continue; // Skip unused 
      ColorExtension& extn = *colorExtensions[iext];
      if (extn.keepers.empty()) continue; // Not using any colors from this catalog
      string filename;
      extensionTable.readCell(filename, "Filename", iext);
      /**/ cerr << "# Reading color catalog " << iext
				<< "/" << colorExtensions.size()
				<< " from " << filename << endl;
      int hduNumber;
      //**      extensionTable.readCell(hduNumber, "HDUNumber", iext);
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
      FTable ff = ft.use();

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
	map<long, Match*>::iterator pr = extn.keepers.find(id[irow]);
	if (pr == extn.keepers.end()) continue; // Not a desired object

	// Have a desired object. Put the color into everything it matches
	Match* m = pr->second;
	extn.keepers.erase(pr);

	for ( Match::iterator j = m->begin();
	      j != m->end();
	      ++j) 
	  (*j)->args.color = color[irow];

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
	
    cerr << "Done reading catalogs for colors." << endl;

    // Make a pass through all matches to reserve as needed and purge 
    // those not wanted.  

    // Also count fitted detections per exposure to freeze parameters of those without any detections
    vector<long> detectionsPerExposure(exposures.size(), 0);

    {
      ran::UniformDeviate u;
      if (randomNumberSeed >0) u.Seed(randomNumberSeed);
      long dcount=0;
      int dof=0;
      double chi=0.;
      double maxdev=0.;

      /**/cerr << "Before reserving, match count: " << matches.size() << endl;
      list<Match*>::iterator im = matches.begin();
      while (im != matches.end() ) {
	// Decide whether to reserve this match
	(*im)->countFit();
	if (reserveFraction > 0.)
	  (*im)->setReserved( u < reserveFraction );

	// Remove Detections with too-large errors:
	Match::iterator j=(*im)->begin();
	while (j != (*im)->end()) {
	  Detection* d = *j;
	  // Keep it if mag error is small or if it's from a tag or reference "instrument"
	  if ( d->sigmaMag > maxMagError
	       && exposures[ extensions[d->catalogNumber]->exposure]->instrument >= 0) {
	    j = (*im)->erase(j, true);  // will delete the detection too.
	  } else {
	    ++j;
	  }
	}
	int nFit = (*im)->fitSize();

	if ((*im)->size() > 0) {
	  // See if color is in range, using color from first Detection (they should
	  // all have the same color).
	  double color = (*(*im)->begin())->args.color;
	  if (color != PhotoArguments::NODATA && (color < minColor || color > maxColor))
	    nFit = 0; // Kill the whole Match below.
	}

	if ( nFit < minMatches) {
	  // Remove entire match if it's too small, and kill its Detections too
	  (*im)->clear(true);
	  im = matches.erase(im);
	} else {
	  // Still a good match.

	  dcount += nFit;	// Count total good detections
	  chi += (*im)->chisq(dof, maxdev);

	  // Count up fitted Detections in each exposure from this Match
	  if (!(*im)->getReserved()) {
	    for (Match::iterator j = (*im)->begin();
		 j != (*im)->end();
		 ++j) {
	      if ((*j)->isClipped || (*j)->wt<=0.) continue; // Not a fitted object
	      int iExtn = (*j)->catalogNumber;
	      int iExpo = extensions[iExtn]->exposure;
	      detectionsPerExposure[iExpo]++;
	    }
	  }

	  ++im;
	}
      } // End loop over all matches

      cout << "Using " << matches.size() 
	   << " matches with " << dcount << " total detections." << endl;
      cout << " chisq " << chi << " / " << dof << " dof maxdev " << maxdev << endl;

    } // End Match-culling section.

    // Check for unconstrained exposures:
    for (int iExposure = 0; iExposure < detectionsPerExposure.size(); iExposure++) {
      if (detectionsPerExposure[iExposure] <= 0 && exposures[iExposure]) {
	Exposure& expo = *exposures[iExposure];
	cerr << "WARNING: Exposure " << expo.name << " has no fitted detections.  Freezing its parameters."
	     << endl;
	mapCollection.setFixed(expo.name, true);
      }
    }

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
	list<PhotoPrior*> thisTime = readPriors(*i, instruments, exposures, extensions, 
						detectionsPerExposure);
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
	  cout << "--Starting strict tolerance passes, clipping full matches" << endl;
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
      nclip = ca.sigmaClip(thresh, false, clipEntireMatch || !coarsePasses);
      if (nclip==0 && coarsePasses) {
	// Nothing being clipped; tighten tolerances and re-fit
	coarsePasses = false;
	ca.setRelTolerance(chisqTolerance);
	cout << "No clipping--Starting strict tolerance passes, clipping full matches" << endl;
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
	cerr << "Error trying to open output file for fitted photometry: " << outPhotFile << endl;
	// *** will not quit before dumping output ***
      } else {
	mapCollection.write(ofs, stringstuff::taggedCommandLine(argc,argv));
      }
    }

    // If there are reserved Matches, run sigma-clipping on them now.
    if (reserveFraction > 0.) {
      cout << "** Clipping reserved matches: " << endl;
      oldthresh = 0.;
      do {
	// Report number of active Matches / Detections in each iteration:
	long int mcount=0;
	long int dcount=0;
	ca.count(mcount, dcount, true, 2);
	double max;
	int dof=0;
	double chisq= ca.chisqDOF(dof, max, true);
	cout << "Clipping " << mcount << " matches with " << dcount << " detections "
	     << " chisq " << chisq << " / " << dof << " dof,  maxdev " << max 
	     << " sigma" << endl;
      
	double thresh = sqrt(chisq/dof) * clipThresh;
	cout << "  new clip threshold: " << thresh << " sigma"
	     << endl;
	if (thresh >= max) break;
	if (oldthresh>0. && (1-thresh/oldthresh)<minimumImprovement) break;
	oldthresh = thresh;
	nclip = ca.sigmaClip(thresh, true, clipEntireMatch);
	cout << "Clipped " << nclip
	     << " matches " << endl;
      } while (nclip>0);
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
    }

    //////////////////////////////////////
    // Output data and calculate some statistics
    //////////////////////////////////////

    // Create Accum instances for fitted and reserved Detections on every
    // exposure, plus total accumulator for all reference
    // and all non-reference objects.
    vector<Accum> vaccFit(exposures.size());
    vector<Accum> vaccReserve(exposures.size());
    Accum refAccFit;
    Accum refAccReserve;
    Accum accFit;
    Accum accReserve;

    // Open the output bintable
    FITS::FitsTable ft(outCatalog, FITS::ReadWrite + FITS::Create, "PhotoOut");
    FTable outTable = ft.use();;

    // Create vectors that will be put into output table
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
    vector<float> xPix;
    outTable.addColumn(xPix, "xPix");
    vector<float> yPix;
    outTable.addColumn(yPix, "yPix");
    vector<float> xExposure;
    outTable.addColumn(xExposure, "xExpo");
    vector<float> yExposure;
    outTable.addColumn(yExposure, "yExpo");
    vector<float> color;
    outTable.addColumn(color, "Color");

    vector<float> magOut;
    outTable.addColumn(magOut, "magOut");
    vector<float> sigMag;
    outTable.addColumn(sigMag, "sigMag");
    vector<float> magRes;
    outTable.addColumn(magRes, "magRes");

    vector<float> wtFrac;
    outTable.addColumn(wtFrac, "wtFrac");

    // Cumulative counter for rows written to table:
    long pointCount = 0;
    // Write vectors to table when this many rows accumulate:
    const long WriteChunk = 100000;

    // Write all matches to output catalog, deleting them along the way
    // and accumulating statistics of each exposure.
    // 
    list<Match*>::iterator im = matches.begin();
    while ( im != matches.end()) {
      // First, write current vectors to table if they've gotten big
      if ( sequence.size() > WriteChunk) {
	outTable.writeCells(sequence, "SequenceNumber", pointCount);
	outTable.writeCells(catalogNumber, "Extension", pointCount);
	outTable.writeCells(objectNumber, "Object", pointCount);
	outTable.writeCells(clip, "Clip", pointCount);
	outTable.writeCells(reserve, "Reserve", pointCount);
	outTable.writeCells(xPix, "xPix", pointCount);
	outTable.writeCells(yPix, "yPix", pointCount);
	outTable.writeCells(xExposure, "xExpo", pointCount);
	outTable.writeCells(yExposure, "yExpo", pointCount);
	outTable.writeCells(color, "Color", pointCount);
	outTable.writeCells(magOut, "magOut", pointCount);
	outTable.writeCells(magRes, "magRes", pointCount);
	outTable.writeCells(sigMag, "sigMag", pointCount);
	outTable.writeCells(wtFrac, "wtFrac", pointCount);

	pointCount += sequence.size();

	sequence.clear();
	catalogNumber.clear();
	objectNumber.clear();
	clip.clear();
	reserve.clear();
	xPix.clear();
	yPix.clear();
	xExposure.clear();
	yExposure.clear();
	color.clear();
	magOut.clear();
	magRes.clear();
	sigMag.clear();
	wtFrac.clear();
      }	// Done flushing the vectors to Table

      Match* m = *im;
      double mean;
      double wtot;
      m->getMean(mean, wtot);
      // Calculate number of DOF per Detection coordinate after
      // allowing for fit to centroid:
      int nFit = m->fitSize();
      double dofPerPt = (nFit > 1) ? 1. - 1./nFit : 0.;
	    
      int detcount=0;
      for (Match::const_iterator idet=m->begin();
	   idet != m->end();
	   ++idet, ++detcount) {
	const Detection* d = *idet;

	// Prepare output quantities
	sequence.push_back(detcount);
	catalogNumber.push_back(d->catalogNumber);
	objectNumber.push_back(d->objectNumber);
	clip.push_back(d->isClipped);
	reserve.push_back(m->getReserved());
	wtFrac.push_back( d->isClipped ? 0. : d->wt / wtot);

	xPix.push_back(d->args.xDevice);
	yPix.push_back(d->args.yDevice);
	xExposure.push_back(d->args.xExposure);
	yExposure.push_back(d->args.yExposure);
	color.push_back(d->args.color);

	magOut.push_back(d->magOut);
	// ?it's updated already, right??
	magRes.push_back( mean==0. ? 0. : d->magOut - mean);
	double sig=d->clipsq;
	sig = (sig > 0.) ? 1./sqrt(sig) : 0.;
	sigMag.push_back(sig);

	dofPerPt = 1. - wtFrac.back();
	
	if (mean!=0.) {
	  // Accumulate statistics for meaningful residuals
	  if (dofPerPt > 0. && !d->isClipped) {
	    int exposureNumber = extensions[d->catalogNumber]->exposure;
	    Exposure* expo = exposures[exposureNumber];
	    if (m->getReserved()) {
	      if (expo->instrument==REF_INSTRUMENT) {
		refAccReserve.add(d, mean, dofPerPt);
		vaccReserve[exposureNumber].add(d, mean, dofPerPt);
	      } else if (expo->instrument==TAG_INSTRUMENT) {
		// do nothing
	      } else {
		accReserve.add(d, mean, dofPerPt);
		vaccReserve[exposureNumber].add(d, mean, dofPerPt);
	      }
	    } else {
	      // Not a reserved object:
	      if (expo->instrument==REF_INSTRUMENT) {
		refAccFit.add(d, mean, dofPerPt);
		vaccFit[exposureNumber].add(d, mean, dofPerPt);
	      } else if (expo->instrument==TAG_INSTRUMENT) {
		// do nothing
	      } else {
		accFit.add(d, mean, dofPerPt);
		vaccFit[exposureNumber].add(d, mean, dofPerPt);
	      }
	    }
	  } // end statistics accumulation

	} // End residuals calculation

      } // End detection loop

      // Done with this match, delete it with its Detections
      m->clear(true);
      // And get rid of match itself.
      im = matches.erase(im);
    } // end match loop

    // Write remaining results to output table:
    outTable.writeCells(sequence, "SequenceNumber", pointCount);
    outTable.writeCells(catalogNumber, "Extension", pointCount);
    outTable.writeCells(objectNumber, "Object", pointCount);
    outTable.writeCells(clip, "Clip", pointCount);
    outTable.writeCells(reserve, "Reserve", pointCount);
    outTable.writeCells(xPix, "xPix", pointCount);
    outTable.writeCells(yPix, "yPix", pointCount);
    outTable.writeCells(xExposure, "xExpo", pointCount);
    outTable.writeCells(yExposure, "yExpo", pointCount);
    outTable.writeCells(color, "Color", pointCount);
    outTable.writeCells(magOut, "magOut", pointCount);
    outTable.writeCells(magRes, "magRes", pointCount);
    outTable.writeCells(sigMag, "sigMag", pointCount);
    outTable.writeCells(wtFrac, "wtFrac", pointCount);

    // Now for standard output, an exposure-by-exposure table.

    // Print summary statistics for each exposure
    cout << "#    Exp         N    RMS   chi_red   xExposure    yExposure \n"
	 << "#                                    |.......degrees......|"
	 << endl;
    // Sort the exposures by name lexical order
    vector<int> ii;
    for (int i=0; i<exposures.size(); i++)
      if (exposures[i])
	ii.push_back(i);
    std::sort(ii.begin(), ii.end(), NameSorter<Exposure>(exposures));
      
    for (int i=0; i<ii.size(); i++) {
      int iexp = ii[i];
      if (!exposures[iexp]) continue;
      cout << "Fit     " << setw(10) << exposures[iexp]->name
	   << " " << vaccFit[iexp].summary()
	   << endl;
      if (reserveFraction > 0. && vaccReserve[iexp].n > 0) 
	  cout << "Reserve " << setw(10) << exposures[iexp]->name
	       << " " << vaccReserve[iexp].summary()
	       << endl;
    } // exposure summary loop
    
    // Output summary data for reference catalog and detections
    cout << "# " << endl;
    cout << "#                   N    DOF  RMS chi_red \n"
	 << "#                                           "
	 << endl;
    cout << "Reference fit     " << refAccFit.summary() << endl;
    if (reserveFraction>0. && refAccReserve.n>0)
      cout << "Reference reserve " << refAccReserve.summary() << endl;

    cout << "Detection fit     " << accFit.summary() << endl;
    if (reserveFraction>0. && accReserve.n>0)
      cout << "Detection reserve " << accReserve.summary() << endl;

    // Cleanup:

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

