// Program to fit astrometric solutions to detection catalogs already matched by WCSFoF.
#include <fstream>
#include <sstream>
#include <map>

#include "Std.h"
#include "Astrometry.h"
#include "Match.h"
#include "FitsTable.h"
#include "StringStuff.h"
#include "Pset.h"
#include "SCAMPMap.h"
#include "Random.h"

using namespace std;
using namespace astrometry;
using namespace stringstuff;
using img::FTable;

string usage=
  "WCSFit: Refine coordinate solutions for a matched set of catalogs.\n"
  "usage: WCSFit <match file> [parameter file] [parameter file...]\n"
  "      <match file>:  FITS file with binary tables produced by WCSFoF\n"
  "     stdin is read as parameters.";


#include "Instrument.h"

class MapMarq {
public:
  MapMarq(const list<Detection*>& dets_, PixelMap* pm_): dets(dets_), pm(pm_) {}
  void operator()(const DVector& a, double& chisq, 
		  DVector& beta, tmv::SymMatrix<double>& alpha);
private:
  const list<Detection*>& dets;
  PixelMap* pm;
};

class Accum {
public:
  double sumx;
  double sumy;
  double sumxx;
  double sumyy;
  double sumxxw;
  double sumyyw;
  double sumxw;
  double sumyw;
  double sumwx;
  double sumwy;
  int n;
  double xpix;
  double ypix;
  double xw;
  double yw;
  double sumdof;
  Accum(): sumxw(0.), sumyw(0.), 
	   sumx(0.), sumy(0.), sumwx(0.), sumwy(0.), 
	   sumxx(0.), sumyy(0.), sumxxw(0.), sumyyw(0.),
	   sumdof(0.), xpix(0.), ypix(0.), xw(0.), yw(0.),
	   n(0) {}
  void add(const Detection* d, double xoff=0., double yoff=0., double dof=1.) {
    sumx += (d->xw-xoff);
    sumy += (d->yw-yoff);
    sumxw += (d->xw-xoff)*d->wtx;
    sumyw += (d->yw-yoff)*d->wty;
    sumxxw += (d->xw-xoff)*(d->xw-xoff)*d->wtx;
    sumyyw += (d->yw-yoff)*(d->yw-yoff)*d->wty;
    sumxx += (d->xw-xoff)*(d->xw-xoff);
    sumyy += (d->yw-yoff)*(d->yw-yoff);
    sumwx += d->wtx;
    sumwy += d->wty;
    sumdof += dof;
    ++n;
  }
  double rms() const {
    return n > 0 ? sqrt( (sumxx+sumyy)/(2.*n)) : 0.;
  }
  double reducedChisq() const {
    return sumdof>0 ? (sumxxw+sumyyw)/(2.*sumdof) : 0.;
  }
  string summary() const {
    ostringstream oss;
    double dx = 0., dy = 0., sigx=0., sigy=0.;
    if (n>0) {
      dx = sumxw / sumwx;
      sigx = 1./sqrt(sumwx);
      dy = sumyw / sumwy;
      sigy = 1./sqrt(sumwy);
    }
    oss << setw(4) << n 
	<< fixed << setprecision(1)
	<< " " << setw(6) << sumdof
	<< " " << setw(5) << dx*1000.*DEGREE/ARCSEC 
	<< " " << setw(5) << sigx*1000.*DEGREE/ARCSEC
	<< " " << setw(5) << dy*1000.*DEGREE/ARCSEC 
	<< " " << setw(5) << sigy*1000.*DEGREE/ARCSEC
	<< " " << setw(5) << rms()*1000.*DEGREE/ARCSEC
	<< setprecision(2) 
	<< " " << setw(5) << reducedChisq()
	<< setprecision(0) << noshowpoint
	<< " " << setw(5) << xpix 
	<< " " << setw(5) << ypix
	<< setprecision(5) << showpoint << showpos
	<< " " << setw(9) << xw 
	<< " " << setw(9) << yw ;
    return oss.str();
  }
};


int
main(int argc, char *argv[])
{
  bool reFit;	// until we separate the statistics part ???
  double reserveFraction;
  double clipThresh;
  double minPixError;
  double maxPixError;
  double minReferenceError;
  int minMatches;
  bool clipEntireMatch;
  int exposureOrder;
  int deviceOrder;
  string outCatalog;
  string catalogSuffix;
  string newHeadSuffix;
  double chisqTolerance;
  Pset parameters;
  {
    const int def=PsetMember::hasDefault;
    const int low=PsetMember::hasLowerBound;
    const int up=PsetMember::hasUpperBound;
    const int lowopen = low | PsetMember::openLowerBound;
    const int upopen = up | PsetMember::openUpperBound;

    parameters.addMemberNoValue("INPUTS");
    parameters.addMember("minPixError",&minPixError, def | low,
			 "Floor for pixel posn uncertainty", 0.01, 0.);
    parameters.addMember("maxPixError",&maxPixError, def | lowopen,
			 "Cut objects with pixel posn uncertainty above this", 0.1, 0.);
    parameters.addMember("minReferenceError",&minReferenceError, def | low,
			 "Floor for ref. posn uncertainty (arcsec)", 0.003, 0.);
    parameters.addMember("minMatch",&minMatches, def | low,
			 "Minimum number of detections for usable match", 2, 2);
    parameters.addMemberNoValue("CLIPPING");
    parameters.addMember("clipThresh",&clipThresh, def | low,
			 "Clipping threshold (sigma)", 5., 2.);
    parameters.addMember("clipEntireMatch",&clipEntireMatch, def,
			 "Discard entire object if one outlier", false);
    parameters.addMemberNoValue("FITTING");
    parameters.addMember("reserveFraction",&reserveFraction, def | low,
			 "Fraction of matches reserved from re-fit", 0., 0.);
    parameters.addMember("exposureOrder",&exposureOrder, def | low,
			 "Order of per-exposure map", 1, 0);
    parameters.addMember("deviceOrder",&deviceOrder, def | low,
			 "Order of per-extension map", 2, 0);
    parameters.addMember("chisqTolerance",&chisqTolerance, def | lowopen,
			 "Fractional change in chisq for convergence", 0.001, 0.);
    parameters.addMemberNoValue("FILES");
    parameters.addMember("outCatalog",&outCatalog, def,
			 "Output FITS binary catalog", "wcscat.fits");
  }

  // Positional accuracy (in degrees) demanded of numerical solutions for inversion of 
  // pixel maps: 
  const double worldTolerance = 0.001*ARCSEC/DEGREE;
  // Fractional reduction in RMS required to continue sigma-clipping:
  const double minimumImprovement=0.02;
  const int REF_INSTRUMENT=-1;	// Instrument for reference objects (no WCS fitting)
  const int TAG_INSTRUMENT=-2;	// Exposure number for tag objects (no WCS fitting nor contrib to stats)

  try {
    
    // Read parameters
    if (argc<2) {
      cerr << usage << endl;
      exit(1);
    }
    string inputTables = argv[1];

    parameters.setDefault();
    for (int i=2; i<argc; i++) {
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

    minReferenceError *= ARCSEC/DEGREE;

    // ??? log parameters to output file somehow?

    /////////////////////////////////////////////////////
    //  Read in properties of all Fields, Instruments, Devices, Exposures
    /////////////////////////////////////////////////////

    // All we care about fields are names and orientations:
    NameIndex fieldNames;
    vector<Orientation> fieldOrients;

    {
      FITS::FitsTable in(inputTables, FITS::ReadOnly, "Fields");
      FITS::FitsTable out(outCatalog, FITS::ReadWrite + FITS::OverwriteFile + FITS::Create, "Fields");
      FTable ft = in.extract();
      out.adopt(ft);
      vector<double> ra;
      vector<double> dec;
      vector<string> name;
      ft.readCells(name, "Name");
      ft.readCells(ra, "RA");
      ft.readCells(dec, "Dec");
      for (int i=0; i<name.size(); i++) {
	fieldNames.append(name[i]);
	fieldOrients.push_back( Orientation(SphericalICRS(ra[i]*DEGREE, dec[i]*DEGREE)));
      }
    }

    // Let's figure out which of our FITS extensions are Instrument or MatchCatalog
    vector<int> instrumentHDUs;
    vector<int> catalogHDUs;
    {
      // Find which extensions are instrument tables
      FITS::FitsFile ff(inputTables);
      for (int i=1; i<ff.HDUCount(); i++) {
	FITS::Hdu h(inputTables, FITS::HDUAny, i);
	if (stringstuff::nocaseEqual(h.getName(), "Instrument"))
	  instrumentHDUs.push_back(i);
	else if (stringstuff::nocaseEqual(h.getName(), "MatchCatalog"))
	  catalogHDUs.push_back(i);
      }
    }
    
    // Read in all the instrument extensions and their device info.
    vector<Instrument*> instruments(instrumentHDUs.size());
    for (int i=0; i<instrumentHDUs.size(); i++) {
      FITS::FitsTable ft(inputTables, FITS::ReadOnly, instrumentHDUs[i]);
      FITS::FitsTable out(outCatalog, FITS::ReadWrite+FITS::Create, "Instrument");
      Assert( stringstuff::nocaseEqual(ft.getName(), "Instrument"));
      FTable ff=ft.extract();
      out.adopt(ff);
      string instrumentName;
      int instrumentNumber;
      if (!ff.header()->getValue("Name", instrumentName)
	  || !ff.header()->getValue("Number", instrumentNumber)) {
	cerr << "Could not read name and/or number of instrument at extension "
	     << instrumentHDUs[i] << endl;
      }
      Assert(instrumentNumber < instruments.size());
      Instrument* instptr = new Instrument(instrumentName);
      instruments[instrumentNumber] = instptr;
      
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
	instptr->addDevice( devnames[j],
			    Bounds<double>( vxmin[j], vxmax[j], vymin[j], vymax[j]));
      }
    }

    // Read in the table of exposures
    vector<Exposure*> exposures;
    {
      FITS::FitsTable ft(inputTables, FITS::ReadOnly, "Exposures");
      FITS::FitsTable out(outCatalog, FITS::ReadWrite+FITS::Create, "Exposures");
      FTable ff = ft.extract();
      out.adopt(ff);
      vector<string> names;
      vector<double> ra;
      vector<double> dec;
      vector<int> fieldNumber;
      vector<int> instrumentNumber;
      ff.readCells(names, "Name");
      ff.readCells(ra, "RA");
      ff.readCells(dec, "Dec");
      ff.readCells(fieldNumber, "fieldNumber");
      ff.readCells(instrumentNumber, "InstrumentNumber");
      for (int i=0; i<names.size(); i++) {
	Exposure* expo = new Exposure(names[i],
				      Orientation(SphericalICRS(ra[i]*DEGREE, 
								dec[i]*DEGREE)));
	expo->field = fieldNumber[i];
	expo->instrument = instrumentNumber[i];
	exposures.push_back(expo);
      }
    }

    // Read and keep the table giving instructions for all input files:
    FTable fileTable;
    {
      FITS::FitsTable ft(inputTables, FITS::ReadOnly, "Files");
      FITS::FitsTable out(outCatalog, FITS::ReadWrite+FITS::Create, "Files");
      fileTable = ft.extract();
      out.copy(fileTable);
    }

    // Read info about all Extensions - we will keep the Table around.
    FTable extensionTable;
    vector<Extension*> extensions;
    {
      FITS::FitsTable ft(inputTables, FITS::ReadOnly, "Extensions");
      extensionTable = ft.extract();
      FITS::FitsTable out(outCatalog, FITS::ReadWrite+FITS::Create, "Extensions");
      out.copy(extensionTable);
    }
    for (int i=0; i<extensionTable.nrows(); i++) {
      extensions.push_back(new Extension);
      Extension* extn = extensions.back();
      int j;
      extensionTable.readCell(j, "ExposureNumber", i);
      extn->exposure = j;
      extensionTable.readCell(j, "DeviceNumber", i);
      extn->device = j;
      string s;
      extensionTable.readCell(s, "WCS", i);
      if (stringstuff::nocaseEqual(s, "ICRS")) {
	// leave startpm null to signal RA/Dec coordinates
      } else {
	istringstream iss(s);
	img::Header h;
	iss >> h;
	extn->startpm = new SCAMPMap(h, &(fieldOrients[exposures[extn->exposure]->field]));
      }
      string wcsOutFile;
      extensionTable.readCell(wcsOutFile, "WCSOut", i);
      extn->wcsOutFile = wcsOutFile;
    }

    /////////////////////////////////////////////////////
    //  Create and initialize all coordinate maps
    /////////////////////////////////////////////////////

    // Create a PixelMapCollection to hold all the components of the PixelMaps.
    PixelMapCollection mapCollection;

    // Create an identity map used by all reference catalogs
    PixelMapKey identityKey = mapCollection.add(new IdentityMap, "Identity");
    PixelMapChain nullChain;	// Can use empty chain to also be identity "map"
    PixelMapChain identityChain; 
    identityChain.append(identityKey);
    SubMap* identitySubMap = mapCollection.issue(identityChain);

    // Set up coordinate maps for each instrument & exposure
    for (int iinst=0; iinst<instruments.size(); iinst++) {
      Instrument* inst = instruments[iinst];
      // First job is to find an exposure that has all devices for this instrument
      int iexp1;
      for (iexp1 = 0; iexp1 < exposures.size(); iexp1++) {
	if (exposures[iexp1]->instrument != iinst) continue; 
	// Filter the extension table for this exposure:
	set<long> vexp;
	vexp.insert(iexp1);
	FTable exts = extensionTable.extractRows("ExposureNumber", vexp);
	// Stop here if this exposure has right number of extensions (exactly 1 per device):
	if (exts.nrows() == inst->nDevices) break; 
      }

      if (iexp1 >= exposures.size()) {
	cerr << "Could not find exposure with all devices for intrument "
	     << inst->name << endl;
	exit(1);
      }

      // Used this exposure to initialize the map for each extension
      Exposure* exp1 = exposures[iexp1];
      // First exposure of instrument will have identity map
      exp1->warp = nullChain;

      // Save reprojection map for first exposure with this instrument:
      // Map from exposure's projection to field's projection
      {
	astrometry::Orientation& fOrient = fieldOrients[exp1->field];
	exp1->reproject = mapCollection.add(new ReprojectionMap( TangentPlane(exp1->orient.getPole(),
									      exp1->orient),
								 TangentPlane(fOrient.getPole(),
									      fOrient),
								 DEGREE),
					    exp1->name + "_reprojection");
      }

      // Create new PolyMap for each device and fit to the starting map of the reference exposure
      const int nGridPoints=400;	// Number of test points for map initialization

      // Get all Extensions that are part of this exposure
      list<Extension*> thisExposure;
      for (int i=0; i<extensions.size(); i++)
	if (extensions[i]->exposure == iexp1)
	  thisExposure.push_back(extensions[i]);
      Assert(thisExposure.size()==inst->nDevices);

      for (int idev=0; idev<inst->nDevices; idev++) {
	// Create coordinate map for each Device, initialize its parameter by fitting to
	// the input map for this device in the reference exposure of this Instrument.
	Bounds<double> b=inst->domains[idev];
	double step = sqrt(b.area()/nGridPoints);
	int nx = static_cast<int> (ceil((b.getXMax()-b.getXMin())/step));
	int ny = static_cast<int> (ceil((b.getYMax()-b.getYMin())/step));
	double xstep = (b.getXMax()-b.getXMin())/nx;
	double ystep = (b.getYMax()-b.getYMin())/ny;
	list<Detection*> testPoints;

	Extension* extn=0;
	for (list<Extension*>::iterator i = thisExposure.begin();
	     i != thisExposure.end();
	     ++i) 
	  if ((*i)->device==idev) {
	    extn = *i;
	    break;
	  }
	if (!extn) {
	  cerr << "Could not find Extension for device " << inst->deviceNames.nameOf(idev)
	       << " in exposure " << exp1->name << endl;
	  exit(1);
	}

	PixelMap* startpm=extn->startpm;
	if (!startpm) {
	  cerr << "Did not read starting WCS for device " << inst->deviceNames.nameOf(idev)
		 << " in exposure " << exp1->name
		 << " with instrument " << inst->name
		 << endl;
	  exit(1);
	}
	const PixelMap* reproject = mapCollection.element(exp1->reproject);

	for (int ix=0; ix<=nx; ix++) {
	  for (int iy=0; iy<=ny; iy++) {
	    double xpix = b.getXMin() + ix*xstep;
	    double ypix = b.getYMin() + iy*ystep;
	    double xw, yw;
	    startpm->toWorld(xpix, ypix, xw, yw);
	    // Undo the change of projections:
	    reproject->toPix(xw,yw,xw,yw);
	    Detection* d = new Detection;
	    d->xpix = xpix;
	    d->ypix = ypix;
	    d->xw = xw;
	    d->yw = yw;
	    // Note MapMarq does not use positional uncertainties.
	    testPoints.push_back(d);
	  }
	}
	PixelMap* pm = new PolyMap(deviceOrder, worldTolerance);
	MapMarq mm(testPoints, pm);
	// Since polynomial fit is linear, should just take one iteration:
	double var=0.;
	DVector beta(pm->nParams(), 0.);
	DVector params(pm->nParams(), 0.);
	tmv::SymMatrix<double> alpha(pm->nParams(), 0.);
	mm(params, var, beta, alpha);
	beta /= alpha;
	params = beta;
	beta.setZero();
	alpha.setZero();
	var = 0.;
	mm(params, var, beta, alpha);
	beta /= alpha;
	params += beta;
	pm->setParams(params);

	// Clear out the testPoints:
	for (list<Detection*>::iterator i = testPoints.begin();
	     i != testPoints.end();
	     ++i)
	  delete *i;

	// Save this map into the collection
	PixelMapKey key = mapCollection.add(pm, inst->name + "_" + inst->deviceNames.nameOf(idev));
	inst->maps[idev].append(key);
	inst->pixelMaps[idev] = mapCollection.issue(inst->maps[idev]);
      } // device loop

      // Now create an exposure map for all other exposures with this instrument,
      // and set initial parameters by fitting to input WCS map.

      int pointsPerDevice = nGridPoints / inst->nDevices + 1;
      for (int iexp=0; iexp<exposures.size(); iexp++) {
	if ( exposures[iexp]->instrument != iinst) continue;
	if (iexp == iexp1) continue;  // Skip the reference exposure

	// Each additional exposure needs a warping tranformation
	Exposure* expo = exposures[iexp];
	PixelMap* warp;
	if (exposureOrder==1) 
	  warp = new LinearMap;
	else
	  warp = new PolyMap(exposureOrder, worldTolerance);

	// Set numerical-derivative step to 1 arcsec, polynomial is working in degrees:
	warp->setPixelStep(ARCSEC/DEGREE);

	// And a reprojection map:
	astrometry::Orientation& fOrient = fieldOrients[expo->field];
	expo->reproject = mapCollection.add(new ReprojectionMap( TangentPlane(expo->orient.getPole(),
									      expo->orient),
								 TangentPlane(fOrient.getPole(),
									      fOrient),
								 DEGREE),
					    expo->name + "_reprojection");
	const PixelMap* reproject = mapCollection.element(expo->reproject);

	// Get all Extensions that are part of this exposure
	thisExposure.clear();
	for (int i=0; i<extensions.size(); i++)
	  if (extensions[i]->exposure == iexp)
	    thisExposure.push_back(extensions[i]);

	// Build new test points
	list<Detection*> testPoints;
	for (list<Extension*>::iterator j=thisExposure.begin();
	     j != thisExposure.end();
	     ++j) {
	  Extension* extn = *j;
	  int idev = extn->device;
	  Bounds<double> b=inst->domains[idev];
	  double step = sqrt(b.area()/pointsPerDevice);
	  int nx = static_cast<int> (ceil((b.getXMax()-b.getXMin())/step));
	  int ny = static_cast<int> (ceil((b.getYMax()-b.getYMin())/step));
	  double xstep = (b.getXMax()-b.getXMin())/nx;
	  double ystep = (b.getYMax()-b.getYMin())/ny;
	  // starting pixel map for the first exposure with this instrument
	  PixelMap* startpm=extn->startpm;
	  if (!startpm) {
	    cerr << "Failed to get starting PixelMap for device "
		 << inst->deviceNames.nameOf(idev)
		 << " in exposure " << expo->name
		 << endl;
	    exit(1);
	  }
	  for (int ix=0; ix<=nx; ix++) {
	    for (int iy=0; iy<=ny; iy++) {
	      double xpix = b.getXMin() + ix*xstep;
	      double ypix = b.getYMin() + iy*ystep;
	      double xw, yw;
	      startpm->toWorld(xpix, ypix, xw, yw);
	      // Undo the change of projections:
	      reproject->toPix(xw,yw,xw,yw);
	      // And process the pixel coordinates through Instrument map:
	      inst->pixelMaps[idev]->toWorld(xpix,ypix,xpix,ypix);
	      Detection* d = new Detection;
	      d->xpix = xpix;
	      d->ypix = ypix;
	      d->xw = xw;
	      d->yw = yw;
	      testPoints.push_back(d);
	    }
	  }
	} // finish device loop for test points

	MapMarq mm(testPoints, warp);
	// Since polynomial fit is linear, should just take one iteration:
	double var=0.;
	DVector beta(warp->nParams(), 0.);
	DVector params(warp->nParams(), 0.);
	tmv::SymMatrix<double> alpha(warp->nParams(), 0.);
	mm(params, var, beta, alpha);
	beta /= alpha;
	params = beta;
	beta.setZero();
	alpha.setZero();
	var = 0.;

	mm(params, var, beta, alpha);
	beta /= alpha;
	params += beta;
	warp->setParams(params);
	
	// Save this warp into the collection and save key for exposure
	expo->warp.append(mapCollection.add(warp, expo->name));

	// Clear out the testPoints:
	for (list<Detection*>::iterator i = testPoints.begin();
	     i != testPoints.end();
	     ++i)
	  delete *i;
      } // End exposure loop for this instrument
    } // end instrument loop.

    // Now construct compound maps for every Extension
    // and save SubMap for it
    for (int iext=0; iext<extensions.size(); iext++) {
      Extension* extn = extensions[iext];
      Exposure* expo = exposures[extn->exposure];
      if ( expo->instrument < 0) {
	// Tag & reference exposures have no instruments and no fitting
	// being done.  Coordinates are fixed to xpix = xw.
	extn->extensionMap = identitySubMap;
      } else {
	// Real instrument, make its compounded map:
	// Start with the instrument map for device:
	PixelMapChain chain = instruments[expo->instrument]->maps[extn->device];
	// Then the exposure's warp
	chain.append(expo->warp);
	// And exposure's reprojection:
	chain.append(expo->reproject);
	extn->extensionMap = mapCollection.issue(chain);
      }
    } // Extension loop

    /**/cerr << "Total number of map elements " << mapCollection.nMaps()
	     << " with " << mapCollection.nParams() << " free parameters."
	     << endl;


    //////////////////////////////////////////////////////////
    // Read in all the data
    //////////////////////////////////////////////////////////

    // List of all Matches - they will hold pointers to all Detections too.
    list<Match*> matches;    // ??? Store objects, not pointers??

    // Start by reading all matched catalogs, creating Detection and Match arrays, and 
    // telling each Extension which objects it should retrieve from its catalog

    for (int icat = 0; icat < catalogHDUs.size(); icat++) {
      /**/cerr << "# Reading catalog extension " << catalogHDUs[icat] << endl;
      FITS::FitsTable ft(inputTables, FITS::ReadOnly, catalogHDUs[icat]);
      FTable ff = ft.use();
      vector<int> seq;
      vector<LONGLONG> extn;
      vector<LONGLONG> obj;
      ff.readCells(seq, "SequenceNumber");
      ff.readCells(extn, "Extension");
      ff.readCells(obj, "Object");
      // Smaller collections for each match
      vector<long> extns;
      vector<long> objs;
      for (int i=0; i<seq.size(); i++) {
	if (seq[i]==0 && !extns.empty()) {
	  // Make a match from previous few entries
	  Detection* d = new Detection;
	  d->catalogNumber = extns[0];
	  d->objectNumber = objs[0];
	  matches.push_back(new Match(d));
	  extensions[extns[0]]->keepers.insert(std::pair<long, Detection*>(objs[0], d));
	  for (int j=1; j<extns.size(); j++) {
	    d = new Detection;
	    d->catalogNumber = extns[j];
	    d->objectNumber = objs[j];
	    matches.back()->add(d);
	    extensions[extns[j]]->keepers.insert(std::pair<long, Detection*>(objs[j], d));
	  }
	  extns.clear();
	  objs.clear();
	} // Finished processing previous match
	extns.push_back(extn[i]);
	objs.push_back(obj[i]);
      } // End loop of catalog entries
      
      // Make a last match out of final entries
      if (!extns.empty()) {
	Detection* d = new Detection;
	d->catalogNumber = extns[0];
	d->objectNumber = objs[0];
	matches.push_back(new Match(d));
	extensions[extns[0]]->keepers.insert(std::pair<long,Detection*>(objs[0], d));
	for (int j=1; j<extns.size(); j++) {
	  Detection* d = new Detection;
	  d->catalogNumber = extns[j];
	  d->objectNumber = objs[j];
	  matches.back()->add(d);
	  extensions[extns[j]]->keepers.insert(std::pair<long, Detection*>(objs[j], d));
	}
      }
    } // End loop over input matched catalogs

    // Now loop over all original catalog bintables, reading the desired rows
    // and collecting needed information into the Detection structures
    for (int iext = 0; iext < extensions.size(); iext++) {
      string filename;
      extensionTable.readCell(filename, "Filename", iext);
      /**/cerr << "# Reading object catalog " << iext
	       << "/" << extensions.size()
	       << " from " << filename << endl;
      int hduNumber;
      extensionTable.readCell(hduNumber, "HDUNumber", iext);
      string xKey;
      extensionTable.readCell(xKey, "xKey", iext);
      string yKey;
      extensionTable.readCell(yKey, "yKey", iext);
      string idKey;
      extensionTable.readCell(idKey, "idKey", iext);
      string errKey;
      extensionTable.readCell(errKey, "errKey", iext);
      double weight;
      extensionTable.readCell(weight, "Weight", iext);

      // Relevant structures for this extension
      Extension* extn = extensions[iext];
      Exposure* expo = exposures[extn->exposure];
      astrometry::Orientation& fOrient = fieldOrients[expo->field];
      Instrument* instr = 0;
      const SubMap* pixmap=extn->extensionMap;
      const PixelMap* startMap = extn->startpm;
      bool isReference = false;
      bool isTag = false;
      if (expo->instrument == REF_INSTRUMENT) {
	isReference = true;
      } else if (expo->instrument == TAG_INSTRUMENT) {
	isTag = true;
      }

      bool useRows = stringstuff::nocaseEqual(idKey, "_ROW");
      // What we need to read from the FitsTable:
      vector<string> neededColumns;
      if (!useRows)
	neededColumns.push_back(idKey);
      neededColumns.push_back(xKey);
      neededColumns.push_back(yKey);
      neededColumns.push_back(errKey);

      FITS::FitsTable ft(filename, FITS::ReadOnly, hduNumber);
      FTable ff = ft.extract(0, -1, neededColumns);
      vector<long> id;
      if (useRows) {
	id.resize(ff.nrows());
	for (long i=0; i<id.size(); i++)
	  id[i] = i;
      } else {
	ff.readCells(id, idKey);
      }
      Assert(id.size() == ff.nrows());

      bool errorColumnIsDouble = true;
      try {
	double d;
	ff.readCell(d, errKey, 0);
      } catch (img::FTableError& e) {
	errorColumnIsDouble = false;
      }
      for (long irow = 0; irow < ff.nrows(); irow++) {
	map<long, Detection*>::iterator pr = extn->keepers.find(id[irow]);
	if (pr == extn->keepers.end()) continue; // Not a desired object

	// Have a desired object now.  Fill its Detection structure
	Detection* d = pr->second;
	extn->keepers.erase(pr);

	d->map = pixmap;
	ff.readCell(d->xpix, xKey, irow);
	ff.readCell(d->ypix, yKey, irow);
	double sigma;
	if (errorColumnIsDouble) {
	  ff.readCell(sigma, errKey, irow);
	} else {
	  float f;
	  ff.readCell(f, errKey, irow);
	  sigma = f;
	}

	if (startMap) {
	  sigma = std::max(minPixError, sigma);
	  d->sigmaPix = sigma;
	  startMap->toWorld(d->xpix, d->ypix, d->xw, d->yw);
	  Matrix22 dwdp = startMap->dWorlddPix(d->xpix, d->ypix);
	  d->clipsqx = pow(d->sigmaPix,-2.) / (dwdp(0,0)*dwdp(0,0)+dwdp(0,1)*dwdp(0,1));
	  d->clipsqy = pow(d->sigmaPix,-2.) / (dwdp(1,0)*dwdp(1,0)+dwdp(1,1)*dwdp(1,1));
	  d->wtx = d->clipsqx * (isTag ? 0. : weight);
	  d->wty = d->clipsqy * (isTag ? 0. : weight);
	    
	} else {
	  // Input coords were degrees RA/Dec; project into Field's tangent plane
	  TangentPlane tp( SphericalICRS(d->xpix*DEGREE, d->ypix*DEGREE), fOrient);
	  double xw, yw;
	  tp.getLonLat(xw,yw);
	  xw /= DEGREE;
	  yw /= DEGREE;  // Back to our WCS system of degrees.

	  // Anything that's coming in as RA/Dec gets "reference" floor on errors
	  sigma = std::max(minReferenceError, sigma);

	  d->xpix = d->xw = xw;
	  d->ypix = d->yw = yw;
	  d->sigmaPix = sigma;
	  // Tags have no weight in fits and don't get clipped:
	  double wt = isTag ? 0. : pow(sigma,-2.);
	  d->clipsqx = wt;
	  d->clipsqy = wt;
	  wt *= weight;  // optional increase in weight for fitting:
	  d->wtx = wt;
	  d->wty = wt;
	}
      } // End loop over catalog objects

      if (!extn->keepers.empty()) {
	cerr << "Did not find all desired objects in catalog " << filename
	     << " extension " << hduNumber
	     << endl;
	exit(1);
      }
    } // end loop over catalogs to read
	
    cerr << "Done reading catalogs." << endl;

    // Make a pass through all matches to reserve as needed and purge 
    // those not wanted.  

    {
      ran::UniformDeviate u;
      long dcount=0;
      int dof=0;
      double chi=0.;
      double maxdev=0.;

      list<Match*>::iterator im = matches.begin();
      while (im != matches.end() ) {
	// Remove Detections with too-large errors:
	Match::iterator j=(*im)->begin();
	while (j != (*im)->end()) {
	  Detection* d = *j;
	  // Keep it if pixel error is small or if it's from a tag or reference "instrument"
	  if ( d->sigmaPix > maxPixError
	       && exposures[ extensions[d->catalogNumber]->exposure]->instrument >= 0) {
	    j = (*im)->erase(j, true);  // will delete the detection too.
	  } else {
	    ++j;
	  }
	}
	int nFit = (*im)->countFit();
	if ( nFit < minMatches) {
	  // Remove entire match if it's too small, and kill its Detections too
	  (*im)->clear(true);
	  im = matches.erase(im);
	} else {
	  // Still a good match.

	  dcount += nFit;	// Count total good detections
	  chi += (*im)->chisq(dof, maxdev);

	  // Decide whether to reserve this match
	  if (reserveFraction > 0.)
	    (*im)->setReserved( u < reserveFraction );

	  ++im;
	}
      } // End loop over all matches

      cout << "Using " << matches.size() 
	   << " matches with " << dcount << " total detections." << endl;
      cout << " chisq " << chi << " / " << dof << " dof maxdev " << sqrt(maxdev) << endl;

    } // End Match-culling section.


    ///////////////////////////////////////////////////////////
    // Now do the re-fitting 
    ///////////////////////////////////////////////////////////

    // make CoordAlign class
    CoordAlign ca(mapCollection, matches);

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
	     << " chisq " << chi << " / " << dof << " dof,  maxdev " << sqrt(maxdev) 
	     << " sigma" << endl;
      }

      // Do the fit here!!
      double chisq = ca.fitOnce();
      // Note that fitOnce() remaps *all* the matches, including reserved ones.
      double max;
      int dof;
      ca.chisqDOF(dof, max, false);	// Exclude reserved Matches
      double thresh = sqrt(chisq/dof) * clipThresh;
      cout << "Final chisq: " << chisq 
	   << " / " << dof << " dof, max deviation " << max
	   << "  new clip threshold at: " << thresh << " sigma"
	   << endl;
      if (thresh >= max || (oldthresh>0. && (1-thresh/oldthresh)<minimumImprovement)) {
	// Sigma clipping is no longer doing much.  Quit if we are at full precision,
	// else require full target precision and initiate another pass.
	if (coarsePasses) {
	  coarsePasses = false;
	  ca.setRelTolerance(chisqTolerance);
	  cout << "--Starting strict tolerance passes" << endl;
	  oldthresh = thresh;
	  nclip = ca.sigmaClip(thresh, false, clipEntireMatch);
	  cout << "Clipped " << nclip
	       << " matches " << endl;
	  continue;
	} else {
	  // Done!
	  break;
	}
      }
      oldthresh = thresh;
      nclip = ca.sigmaClip(thresh, false, clipEntireMatch);
      if (nclip==0 && coarsePasses) {
	// Nothing being clipped; tighten tolerances and re-fit
	coarsePasses = false;
	ca.setRelTolerance(chisqTolerance);
	cout << "No clipping--Starting strict tolerance passes" << endl;
	continue;
      }
      cout << "Clipped " << nclip
	   << " matches " << endl;
      long int dcount=0;
      long int mcount=0;
      
    } while (coarsePasses || nclip>0);
  
    // The re-fitting is now complete.  Save the new header files
    
    // Keep track of those we've already written to:
    set<string> alreadyOpened;

    const double SCAMPTolerance=0.0001*ARCSEC/DEGREE;  // accuracy of SCAMP-format fit to real solution
    for (int iext=0; iext<extensions.size(); iext++) {
      Extension& extn = *extensions[iext];
      Exposure& expo = *exposures[extn.exposure];
      if (expo.instrument < 0) continue;
      Instrument& inst = *instruments[expo.instrument];
      // Open file for new ASCII headers  ??? Change the way these are done!!!
      string newHeadName = extn.wcsOutFile;
      bool overwrite = false;
      if ( alreadyOpened.find(newHeadName) == alreadyOpened.end()) {
	overwrite = true;
	alreadyOpened.insert(newHeadName);
      }
      ofstream newHeads(newHeadName.c_str(), overwrite ? ios::out : (ios::out | ios::ate));
      if (!newHeads) {
	cerr << "WARNING: could not open new header file " << newHeadName << endl;
	cerr << "...Continuing with program" << endl;
	continue;
      }
      PixelMap* thisMap = extn.extensionMap;
      img::Header h = FitSCAMP(inst.domains[extn.device],
			       *thisMap,
			       fieldOrients[expo.field],
			       expo.orient.getPole(),
			       SCAMPTolerance);
      
      newHeads << h;
    }

    // If there are reserved Matches, run sigma-clipping on them now.
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
	   << " chisq " << chisq << " / " << dof << " dof,  maxdev " << sqrt(max) 
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
    FITS::FitsTable ft(outCatalog, FITS::ReadWrite + FITS::Create, "WCSOut");
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
    vector<float> xpix;
    outTable.addColumn(xpix, "xPix");
    vector<float> ypix;
    outTable.addColumn(ypix, "yPix");
    vector<float> sigpix;
    outTable.addColumn(sigpix, "sigPix");
    vector<float> xrespix;
    outTable.addColumn(xrespix, "xresPix");
    vector<float> yrespix;
    outTable.addColumn(yrespix, "yresPix");

    vector<float> xw;
    outTable.addColumn(xw, "xW");
    vector<float> yw;
    outTable.addColumn(yw, "yW");
    vector<float> sigw;
    outTable.addColumn(sigw, "sigW");
    vector<float> xresw;
    outTable.addColumn(xresw, "xresW");
    vector<float> yresw;
    outTable.addColumn(yresw, "yresW");

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
	outTable.writeCells(xpix, "xPix", pointCount);
	outTable.writeCells(ypix, "yPix", pointCount);
	outTable.writeCells(xrespix, "xresPix", pointCount);
	outTable.writeCells(yrespix, "yresPix", pointCount);
	outTable.writeCells(xw, "xW", pointCount);
	outTable.writeCells(yw, "yW", pointCount);
	outTable.writeCells(xresw, "xresW", pointCount);
	outTable.writeCells(yresw, "yresW", pointCount);
	outTable.writeCells(sigpix, "sigPix", pointCount);
	outTable.writeCells(sigw, "sigW", pointCount);
	outTable.writeCells(wtFrac, "wtFrac", pointCount);

	pointCount += sequence.size();

	sequence.clear();
	catalogNumber.clear();
	objectNumber.clear();
	clip.clear();
	reserve.clear();
	xpix.clear();
	ypix.clear();
	xrespix.clear();
	yrespix.clear();
	xw.clear();
	yw.clear();
	xresw.clear();
	yresw.clear();
	sigpix.clear();
	sigw.clear();
	wtFrac.clear();
      }	// Done flushing the vectors to Table

      Match* m = *im;
      double xc, yc;
      double wtotx, wtoty;
      m->centroid(xc, yc, wtotx, wtoty);
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
	wtFrac.push_back( d->isClipped ? 0. : (d->wtx + d->wty) / (wtotx+wtoty));  //Just take mean of x & y

	xpix.push_back(d->xpix);
	ypix.push_back(d->ypix);
	sigpix.push_back(d->sigmaPix);

	xw.push_back(d->xw);
	yw.push_back(d->yw);
	double sigmaWorld=0.5*(d->clipsqx + d->clipsqy);
	sigmaWorld = (sigmaWorld > 0.) ? 1./sqrt(sigmaWorld) : 0.;
	sigw.push_back(sigmaWorld);

	// Calculate residuals if we have a centroid for the match:
	double xcpix=0., ycpix=0.;
	double xerrw=0., yerrw=0.;
	double xerrpix=0., yerrpix=0.;

	if (xc!=0. || yc!=0.) {
	  xcpix = d->xpix;
	  ycpix = d->ypix;
	  Assert (d->map);
	  d->map->toPix(xc, yc, xcpix, ycpix);

	  xerrw = d->xw - xc;
	  yerrw = d->yw - yc;
	  xerrpix = d->xpix - xcpix;
	  yerrpix = d->ypix - ycpix;

	  // Accumulate statistics for meaningful residuals
	  if (dofPerPt >= 0. && !d->isClipped) {
	    int exposureNumber = extensions[d->catalogNumber]->exposure;
	    Exposure* expo = exposures[exposureNumber];
	    if (m->getReserved()) {
	      if (expo->instrument==REF_INSTRUMENT) {
		refAccReserve.add(d, xc, yc, dofPerPt);
		vaccReserve[exposureNumber].add(d, xc, yc, dofPerPt);
	      } else if (expo->instrument==TAG_INSTRUMENT) {
		// do nothing
	      } else {
		accReserve.add(d, xc, yc, dofPerPt);
		vaccReserve[exposureNumber].add(d, xc, yc, dofPerPt);
	      }
	    } else {
	      // Not a reserved object:
	      if (expo->instrument==REF_INSTRUMENT) {
		refAccFit.add(d, xc, yc, dofPerPt);
		vaccFit[exposureNumber].add(d, xc, yc, dofPerPt);
	      } else if (expo->instrument==TAG_INSTRUMENT) {
		// do nothing
	      } else {
		accFit.add(d, xc, yc, dofPerPt);
		vaccFit[exposureNumber].add(d, xc, yc, dofPerPt);
	      }
	    }
	  } // end statistics accumulation

	} // End residuals calculation

	// Put world residuals into milliarcsec
	xrespix.push_back(xerrpix);
	yrespix.push_back(yerrpix);
	xresw.push_back(xerrw*1000.*DEGREE/ARCSEC);
	yresw.push_back(yerrw*1000.*DEGREE/ARCSEC);

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
    outTable.writeCells(xpix, "xPix", pointCount);
    outTable.writeCells(ypix, "yPix", pointCount);
    outTable.writeCells(xrespix, "xresPix", pointCount);
    outTable.writeCells(yrespix, "yresPix", pointCount);
    outTable.writeCells(xw, "xW", pointCount);
    outTable.writeCells(yw, "yW", pointCount);
    outTable.writeCells(xresw, "xresW", pointCount);
    outTable.writeCells(yresw, "yresW", pointCount);
    outTable.writeCells(sigpix, "sigPix", pointCount);
    outTable.writeCells(sigw, "sigW", pointCount);
    outTable.writeCells(wtFrac, "wtFrac", pointCount);

    // Print summary statistics for each exposure
    cout << "#    Exp    N    DOF    dx    +-    dy    +-   RMS chi_red  "
      "xpix  ypix      xw       yw \n"
	 << "#                     |.....milliarcsec...........|                     |....degrees....|"
	 << endl;
    for (int iexp=0; iexp<exposures.size(); iexp++) {
      cout << "Fit     " << setw(3) << iexp
	   << " " << vaccFit[iexp].summary()
	   << endl;
      if (reserveFraction > 0. && vaccReserve[iexp].n > 0) 
	  cout << "Reserve " << setw(3) << iexp
	       << " " << vaccReserve[iexp].summary()
	       << endl;
    } // exposure summary loop
    
    // Output summary data for reference catalog and detections
    cout << "# " << endl;
    cout << "#                   N    DOF    dx    +-    dy    +-   RMS chi_red  \n"
	 << "#                             |.....milliarcsec...........|         "
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
      delete extensions[i];
    // Get rid of exposures
    for (int i=0; i<exposures.size(); i++)
      delete exposures[i];
    // Get rid of instruments
    for (int i=0; i<instruments.size(); i++)
      delete instruments[i];

  } catch (std::runtime_error& m) {
    quit(m,1);
  }
}
    
void
MapMarq::operator()(const DVector& a, 
		    double& chisq, 
		    DVector& beta, 
		    tmv::SymMatrix<double>& alpha) {
  chisq = 0;
  beta.setZero();
  alpha.setZero();
  int np = a.size();
  Assert(pm->nParams()==np);
  Assert(beta.size()==np);
  Assert(alpha.nrows()==np);
  DMatrix derivs(2,np);
  pm->setParams(a);

  for (list<Detection*>::const_iterator i = dets.begin();
       i != dets.end();
       ++i) {
    const Detection* d = *i;
    if (d->isClipped) continue;
    double xmod, ymod;
    pm->toWorld(d->xpix, d->ypix, xmod, ymod, derivs);
    xmod = d->xw - xmod;
    ymod = d->yw - ymod;

    chisq += xmod*xmod + ymod*ymod;
    for (int i=0; i<np; i++) {
      beta[i] += xmod*derivs(0,i) + ymod*derivs(1,i);
      for (int j=0; j<=i; j++) 
	alpha(i,j)+=derivs(0,i)*derivs(0,j) + derivs(1,i)*derivs(1,j);
    }
  }
}
