// Program to fit astrometric solutions to detection catalogs already matched by WCSFoF.
#include <fstream>
#include <sstream>
#include <map>

#include "Std.h"
#include "Astrometry.h"
#include "Match.h"
#include "Image.h"
#include "StringStuff.h"
#include "Pset.h"
#include "HeaderFromStream.h"
#include "SCAMPMap.h"
#include "Random.h"

using namespace std;
using namespace astrometry;
using namespace stringstuff;

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
			 "Catalog of output points, if desired", "test.out");
    parameters.addMember("catalogSuffix",&catalogSuffix, def,
			 "Suffix for input ASCII catalog files", ".txt");
    parameters.addMember("newHeadSuffix",&newHeadSuffix, def,
			 "Suffix for new ASCII header files", ".head2");
  }

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

    // Open up output catalog file if desired ???
    ofstream ocatfile;
    if (!outCatalog.empty()) {
      ocatfile.open(outCatalog.c_str());
      if (!ocatfile) {
	cerr << "Could not open output catalog " << outCatalog << endl;
	exit(1);
      }
      ocatfile << stringstuff::taggedCommandLine(argc, argv) << endl;
      parameters.dump(ocatfile);
    }


    /////////////////////////////////////////////////////
    //  Read in properties of all Fields, Instruments, Devices, Exposures
    /////////////////////////////////////////////////////

    // All we care about fields are names and orientations:
    NameIndex fieldNames;
    vector<Orientation> fieldOrients;

    {
      FitsTable in(inputTables, FITS::ReadOnly, "Fields");
      FTable ft = in.use();
      vector<double> ra;
      vector<double> dec;
      vector<string> name;
      ft.readCells(name, "Name");
      ft.readCells(ra, "RA");
      ft.readCells(dec, "Dec");
      for (int i=0; i<name.size(); i++) {
	fieldNames.append(name[i]);
	fieldOrients.append( Orientation(SphericalICRS(ra[i]*DEGREE, dec[i]*DEGREE)));
      }
    }

    // Let's figure out which of our FITS extensions are Instrument or MatchCatalog
    vector<int> instrumentHDUs;
    vector<int> catalogHDUs;
    {
      // Find which extensions are instrument tables
      FitsFile ff(inputTables);
      for (int i=1; i<ff.HDUCount(); i++) {
	Hdu h(inputTables, FITS::HduAny, i);
	if (stringstuff::nocaseEqual(h.getName(), "Instrument"))
	  instrumentHDUs.push_back(i);
	else if (stringstuff::nocaseEqual(h.getName(), "MatchCatalog"))
	  catalogHDUs.push_back(i);
      }
    }
    
    // Read in all the instrument extensions and their device info.
    vector<Instrument*> instruments(instrumentHDUs.size(),
				    new Instrument("dummy"));
    for (int i=0; i<instrumentHDUs.size(); i++) {
      FitsTable ft(inputTables, FITS::ReadOnly, instrumentHDUs[i]);
      Assert( stringstuff::nocaseEqual(ft.getName(), "Instrument"));
      FTable ff=ft.use();
      string instrumentName;
      int instrumentNumber;
      if (!ff.header()->getValue(instrumentName, "Name")
	  || !ff.header()->getValue(instrumentNumber, "Number")) {
	cerr << "Could not read name and/or number of instrument at extension "
	     << instrumentHDUs[i] << endl;
      }
      Assert(instrumentNumber < instruments.size());
      Instrument* instptr = instruments[instrumentNumber];
      instptr->name = instrumentName;
      
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
      for (int j=0; j<devnames.size(); j++)
	instptr->addDevice( devnames[j],
			    Bounds<double>( vxmin[j], vxmax[j], vymin[j], vymax[j]));
    }

    // Read in the table of exposures
    vector<Exposure*> exposures;
    {
      FitsTable ft(inputTables, FITS::ReadOnly, "Exposures");
      FTable ff = ft.use();
      vector<string> names;
      vector<double> ra;
      vector<double> dec;
      vector<int> fieldNumber;
      vector<int> instrumentNumber;
      ff.readCells(names, "Name");
      ff.readCells(ra, "RA");
      ff.readCells(dec, "Dec");
      ff.readCells(fieldNumber, "fieldNumber");
      ft.readCells(instrumentNumber, "InstrumentNumber");
      for (int i=0; i<names.size(); i++) {
	Exposure* expo = new Exposure(names[i],
				      Orientation(SphericalICRS(ra[i]*DEGREE, 
								dec[i]*DEGREE)));
	expo->field = fieldNumber[i];
	expo->instrument = instrumentNumber[j];
	exposures.push_back(expo);
      }
    }

    // Read and keep the table giving instructions for all input files:
    FTable fileTable;
    {
      FitsTable ft(inputTables; FITS::ReadOnly, "Files");
      fileTable = ft.extract();
    }

    // Read info about all Extensions - we will keep the Table around.
    FTable extensionTable;
    vector<Extension*> extensions;
    {
      FitsTable ft(inputTables, FITS::ReadOnly, "Extensions");
      extensions = ft.extract();
    }
    for (int i=0; i<extensions.nrows(); i++) {
      extensions.push_back(new Extension);
      Extension* extn = extensions.back();
      int j;
      ft.readCell(j, "ExposureNumber");
      extn->exposure = j;
      ft.readCell(j, "DeviceNumber");
      extn->device = j;
      string s;
      ft.readCell(s, "WCS");
      if (stringstuff::nocaseEqual(s, "ICRS")) {
	// leave startpm null to signal RA/Dec coordinates
      } else {
	istringstream iss(s);
	img::Header h;
	h << iss;
	extn->startpm = new SCAMPMap(h, &(fieldOrients[extn->exposure->field]));
      }
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
	vector<int> vexp(1,iexp1);
	FTable exts = extensionTable.extract("ExposureNumber", vexp);
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
	const Field& fp = fields[exp1->field];
	exp1->reproject = mapCollection.add(new ReprojectionMap( TangentPlane(exp1->orient->getPole(),
									      *exp1->orient),
								 TangentPlane(fp.orient->getPole(),
									      *fp.orient),
								 DEGREE),
					    exp1->name + "_reprojection");
      }

      // Create new PolyMap for each device and fit to the starting map of the reference exposure
      const int nGridPoints=400;	// Number of test points for map initialization

      // Get all Extensions that are part of this exposure
      list<Extension*> thisExposure;
      for (int i=0; i<exposures.size(); i++)
	if (exposures[i]->exposure = iexp1)
	  thisExposure.push_back(i);
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
	  cerr << "Did not read starting WCS for device " << inst->deviceNames.nameOf(iext)
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
	PixelMapKey key = mapCollection.add(pm, inst->name + "_" + inst->deviceNames.nameOf(iext));
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
	const Field& fp = fields[expo->field];
	expo->reproject = mapCollection.add(new ReprojectionMap( TangentPlane(expo->orient.getPole(),
									      expo->orient),
								 TangentPlane(fp.orient.getPole(),
									      fp.orient),
								 DEGREE),
					    expo->name + "_reprojection");
	const PixelMap* reproject = mapCollection.element(expo->reproject);

	// Get all Extensions that are part of this exposure
	thisExposure.clear();
	for (int i=0; i<exposures.size(); i++)
	  if (exposures[i]->exposure = iexp)
	    thisExposure.push_back(i);

	// Build new test points
	list<Detection*> testPoints;
	for (list<Extenstion*>::iterator j=thisExposure.begin();
	     j != thisExposure.end();
	     ++j) {
	  int iext = *j;
	  Extension* extn = extensions[iext];
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
      if ( expo->instrument == REF_INSTRUMENT) {
	// This is a reference catalog, use identity map:
	extn->extensionMap = identitySubMap;
      } else if (exp->instrument == TAG_INSTRUMENT) {
	// Tag catalog - will not be fit, signal with null SubMap
	extn->extensionMap = 0;
      } else {
	// Real instrument, make its compounded map:
	// Start with the instrument map for device:
	PixelMapChain chain = inst->maps[extn->device];
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
      FitsTable ft(inputTables, FITS::ReadOnly, catalogHDUs[icat]);
      FTable ff = ft.use();
      vector<int> sequence;
      vector<long> extension;
      vector<long> object;
      ff.readCells(sequence, "Sequence");
      ff.readCells(extension, "Extension");
      ff.readCells(object, "Object");
      // Smaller collections for each match
      vector<long> extns;
      vector<long> objs;
      for (int i=0; i<sequence.size(); i++) {
	if (sequence[i]==0 && !extns.empty()) {
	  // Make a match from previous few entries
	  Detection* d = new Detection;
	  d->catalogNumber = extns[0];
	  d->objectNumber = objs[0];
	  matches.push_back(new Match(d));
	  extensions[extns[0]].keepers.insert(std::pair(objs[0], d));
	  for (int j=1; j<extns.size(); j++) {
	    d = new Detection;
	    d->catalogNumber = extns[j];
	    d->objectNumber = objs[j];
	    matches.back()->add(d);
	    extensions[extns[j]].keepers.insert(std::pair(objs[j], d));
	  }
	  extns.clear();
	  objs.clear();
	} // Finished processing previous match
	extns.push_back(extension[i]);
	objs.push_back(object[i]);
      } // End loop of catalog entries
      
      // Make a last match out of final entries
      if (!extns.empty()) {
	detections.push_back(Detection());
	Detection* d = &(detections.back());
	d->catalogNumber = extns[0];
	d->objectNumber = objs[0];
	matches.push_back(new Match(d));
	extensions[extns[0]].keepers.insert(std::pair(objs[0], d));
	for (int j=1; j<extns.size(); j++) {
	  detections.push_back(Detection());
	  Detection* d = &(detections.back());
	  d->catalogNumber = extns[j];
	  d->objectNumber = objs[j];
	  matches.back()->add(d);
	  extensions[extns[j]].keepers.insert(std::pair(objs[j], d));
	}
      }
    } // End loop over input matched catalogs

    // Now loop over all original catalog bintables, reading the desired rows
    // and collecting needed information into the Detection structures
    for (int iext = 0; iext < extensions.size(); iext++) {
      string filename;
      extensionTable.readCell(filename, "Filename", iext);
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
      const Field& fld = fields[expo->field];
      Instrument* instr = 0;
      const SubMap* map=extn->extensionMap;
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

      FitsTable ft(filename, FITS::ReadOnly, hduNumber);
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

      for (long irow = 0; irow < ff.nrows(); irow++) {
	map<long, Detection*>::iterator pr = extn->keepers.find(id[irow]);
	if (pr == extn->keepers.end()) continue; // Not a desired object

	// Have a desired object now.  Fill its Detection structure
	Detection* d = pr->second;
	extn->keepers.erase(pr);

	d->map = map;
	ff.readCell(d->xpix, xKey, irow);
	ff.readCell(d->ypix, yKey, irow);
	double sigma;
	ff.readCell(sigma, idKey, irow);
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
	  TangentPlane tp( SphericalICRS(d->xpix*DEGREE, d->ypix*DEGREE), fld.orient);
	  double xw, yw;
	  tp.getLonLat(xw,yw);
	  xw /= DEGREE;
	  yw /= DEGREE;  // Back to our WCS system of degrees.

	  // Anything that's coming in as RA/Dec gets "reference" floor on errors
	  sigma = std::max(minReferenceError, sigma);

	  d->xpix = d->xw = xw;
	  d->ypix = d->yw = yw;
	  d->sigmaPix = sigma;
	  double wt = pow(sigma,-2.);
	  d->clipsqx = wt;
	  d->clipsqy = wt;
	  wt *= isTag ? 0. : weight;   // Tags have no weight in fits.
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

      list<Match*>::iterator i = matches.begin();
      while (i != matches.end() ) {
	// Remove Detections with too-large errors:
	int ngood = 0;
	list<Detection*> kill;
	for (Match::iterator j=(*i)->begin();
	     j != (*i)->end();
	     ++j) {
	  Detection* d = *j;
	  // Keep it if pixel error is small or if it's from a tag or reference "instrument"
	  if ( d->sigmaPix > maxPixError
	       && exposures[ extensions[d->catalogNumber]->exposure]->instrument >= 0) {
	    kill.push_back(d);
	  } else {
	    // count Detections with non-zero weights
	    if ( d->wtx > 0.) ngood++;
	  }
	}
	if ( ngood < minMatches) {
	  // Remove entire match if it's too small, and kill its Detections too
	  (*i)->clear(true);
	  i = matches.erase(i);
	} else {
	  // Still a good match.  Kill the detections above maxPixError
	  for (list<Detection*>::iterator j=kill.begin();
	       j != kill.end();
	       ++j) {
	    (*i)->erase(*j);
	    delete *j;
	  }

	  dcount += ngood;	// Count total good detections
	  chi += (*i)->chisq(dof, maxdev);

	  // Decide whether to reserve this match
	  if (reserveFraction > 0.)
	    (*i)->setReserved( u < reserveFraction );

	  ++i;
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
      
      // ??? check on count()
      ca.count(mcount, dcount, false, 2, -2, -1);
      cout << "-1 / 0: " << dcount << endl;
        for (int iexp=0; iexp<exposures.size(); iexp++) {
          int iinst = exposures[iexp]->instrument;
          for (int iext=0; iext<instruments[iinst]->nExtensions; iext++) {
                ca.count(mcount, dcount, false, 2, -2, iexp, iext);
                cout << iexp << " / " << iext << ": " << dcount << endl;
          }
        }
    } while (coarsePasses || nclip>0);

    // The re-fitting is now complete.  Save the new header files
    const double SCAMPTolerance=0.0001*ARCSEC/DEGREE;  // accuracy of SCAMP-format fit to real solution
      for (int iexp=0; iexp<exposures.size(); iexp++) {
	int iinst = exposures[iexp]->instrument;
	for (int iext=0; iext<instruments[iinst]->nExtensions; iext++) {
	  // Open file for new ASCII headers
	    string newHeadName = exposures[iexp]->name + "_" 
	      + instruments[iinst]->extensionNames.nameOf(iext)
	      + newHeadSuffix;
	    ofstream newHeads(newHeadName.c_str());
	    if (!newHeads) {
	      cerr << "WARNING: could not open new header file " << newHeadName << endl;
	      cerr << "...Continuing with program" << endl;
	      continue;
	    }	    
	    PixelMap* thisMap = exposures[iexp]->extensionMaps[iext];
	    img::ImageHeader h = FitSCAMP(instruments[iinst]->domains[iext],
					  *thisMap,
					  *(fields[exposures[iexp]->field].orient),
					  exposures[iexp]->orient->getPole(),
					  SCAMPTolerance);
	    newHeads << h;
	}
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

    // First tell output file the codes for exposures
    if (ocatfile) {
      ocatfile << "# Exposure code guide: " << endl;
      ocatfile << "# Code   name    field   instrument " << endl;
      for (int iexp=0; iexp<exposures.size(); iexp++) {
	ocatfile << "# " << setw(3) << iexp 
		 << " " << setw(10) << exposures[iexp]->name
		 << " " << setw(10) << fieldNames.nameOf(exposures[iexp]->field)
		 << " " << setw(10) << instrumentNames.nameOf(exposures[iexp]->instrument)
		 << endl;
      }
    } // Done with code number guide


    // Create Accum instances for fitted and reserved Detections on every
    // extension of every exposure, plus total accumulator for all reference
    // and all non-reference objects.
    vector< vector<Accum> > vvaccFit(exposures.size());
    vector< vector<Accum> > vvaccReserve(exposures.size());
    Accum refAccFit;
    Accum refAccReserve;
    Accum accFit;
    Accum accReserve;

    for (int iexp=0; iexp<exposures.size(); iexp++) {
      int iinst = exposures[iexp]->instrument;
      vvaccFit[iexp] = vector<Accum>(instruments[iinst]->nExtensions);
      vvaccReserve[iexp] = vector<Accum>(instruments[iinst]->nExtensions);
      for (int iext=0; iext<instruments[iinst]->nExtensions; iext++) {
	// Skip if this exposure has no map (= no objects) in this extension
	if (iext >= exposures[iexp]->startpm.size()
	    || !(exposures[iexp]->startpm[iext])) continue;
	Position<double> p=instruments[iinst]->domains[iext].center();
	vvaccFit[iexp][iext].xpix = p.x;
	vvaccFit[iexp][iext].ypix = p.y;
	double xw, yw;
	exposures[iexp]->startpm[iext]->toWorld(p.x,p.y,xw,yw);
	vvaccFit[iexp][iext].xw = xw;
	vvaccFit[iexp][iext].yw = yw;
	vvaccReserve[iexp][iext] = vvaccFit[iexp][iext];
      }
    }

    // Column header for output file
    if (ocatfile)
      ocatfile << "#(1)(2)(3)  (4) (5)  (6)   (7)   (8)    (9)     (10)    (11)  (12) "
	"     (13)      (14)      (15) (16) (17)\n"
	       << "#ct/res/clip ID Expo Extn wtfrac xpix  ypix     xres   yres  sigpix"
	"     xw         yw       xres yres sigw\n"
	       << "#                               |............pixels................|"
	"........degrees.......|..milliarcsec...|"
	       << endl;

    // Total up deviations for all extensions, cycling through matches
    for (int ifield=0; ifield<fields.size(); ifield++) {
      for (map<int, MCat>::iterator iaff = fields[ifield].affinities.begin();
	   iaff != fields[ifield].affinities.end();
	   ++iaff) {
	MCat& mc = iaff->second;
	for (MCat::iterator icat=mc.begin();
	     icat!=mc.end();
	     ++icat) {
	  Match& m = *(*icat);
	  double xc, yc;
	  double wtotx, wtoty;
	  m.centroid(xc, yc, wtotx, wtoty);
	  int detcount=0;
	  // Calculate number of DOF per Detection coordinate after
	  // allowing for fit to centroid:
	  int ngood = m.goodSize();
	  double dofPerPt = (ngood > 1) ? 1. - 1./ngood : 0.;
	    
	  for (Match::const_iterator idet=m.begin();
	       idet != m.end();
	       ++idet, ++detcount) {
	    const Detection* d = *idet;

	    if (ocatfile) {
	      // Prepare output quantities
	      double xcpix=0., ycpix=0.;
	      double xresid=0., yresid=0.;
	      double xresidpix=0., yresidpix=0.;
	      double wtfrac = (d->isClipped) ? 0. : (d->wtx + d->wty) / (wtotx+wtoty);  //Just take mean of x & y
	      double sigmaWorld=0.5*(d->clipsqx + d->clipsqy);
	      sigmaWorld = (sigmaWorld > 0.) ? 1./sqrt(sigmaWorld) : 0.;
	      // Get extension name, unless it's a reference object:
	      string eName = (d->exposure < 0) ?
		"00" :
		instruments[exposures[d->exposure]->instrument]->extensionNames.nameOf(d->extension);
	      if (xc!=0. || yc!=0.) {
		// Have a valid centroid, map back to pixel coordinates
		xcpix = d->xpix;
		ycpix = d->ypix;
		if (d->map) {
		  // Note that reference-catalog objects always come here
		  d->map->toPix(xc, yc, xcpix, ycpix);
		} else if (exposures[d->exposure]->startpm[d->extension]) {
		  // use input map for this extension
		  exposures[d->exposure]->startpm[d->extension]->toPix(xc, yc, xcpix, ycpix);
		} else {
		  cerr << "No valid pixel map for detection " << d->id
		       << " exposure " << d->exposure
		       << " extension " << eName
		       << endl;
		  continue;
		}
		xresid = d->xw - xc;
		yresid = d->yw - yc;
		xresidpix = d->xpix - xcpix;
		yresidpix = d->ypix - ycpix;
	      }
	      ocatfile << right
		       << setw(2) << detcount
		       << " " << setw(1) << m.getReserved()
		       << " " << setw(1) << d->isClipped
		       << " " << setw(8) << d->id
		       << " " << setw(3) << d->exposure
		       << " " << setw(3) << eName
		       << " " << setprecision(2) << fixed << setw(4) << wtfrac
		// pixel coords, residuals, and sigma
		       << " " << setprecision(3) << fixed 
		       << setw(8) << d->xpix
		       << " " << setw(8) << d->ypix
		       << showpos
		       << " " << setw(5) << xresidpix
		       << " " << setw(5) << yresidpix
		       << noshowpos
		       << " " << setw(5) << d->sigmaPix
		// world coords (usually degrees), residuals, and sigma (mas)
		       << " " << setprecision(7) << showpos
		       << setw(11) << d->xw
		       << " " << setw(11) << d->yw
		       << setprecision(0) 
		       << " " << setw(4) << xresid*1000.*DEGREE/ARCSEC
		       << " " << setw(4) << yresid*1000.*DEGREE/ARCSEC
		       << noshowpos
		       << " " << setw(3) << sigmaWorld*1000.*DEGREE/ARCSEC
		// ancillary info
		       << " " << d->ancillary
		       << endl;
	    } // end output file
	      // Accumulate data for this detection if it's not clipped or a singlet
	    if (d->isClipped || dofPerPt <= 0.) continue;
	    if (m.getReserved()) {
	      if (d->exposure==REF_EXPOSURE) refAccReserve.add(d, xc, yc, dofPerPt);
	      else if (d->exposure==TAG_EXPOSURE) {
		// do nothing
	      } else {
		accReserve.add(d, xc, yc, dofPerPt);
		vvaccReserve[d->exposure][d->extension].add(d, xc, yc, dofPerPt);
	      }
	    } else {
	      // Not a reserved object:
	      if (d->exposure==REF_EXPOSURE) refAccFit.add(d, xc, yc, dofPerPt);
	      else if (d->exposure==TAG_EXPOSURE) {
		// do nothing
	      } else {
		accFit.add(d, xc, yc, dofPerPt);
		vvaccFit[d->exposure][d->extension].add(d, xc, yc, dofPerPt);
	      }
	    }
	  } // end detection loop
	} // end match loop
      } // end affinity loop
    } // end field loop

      // Done accumulating statistics and making output catalog.
      // Now print summary statistics to stdout for each chip of each exposure

    cout << "#    Expo Extn    N    DOF    dx    +-    dy    +-   RMS chi_red  "
      "xpix  ypix      xw       yw \n"
	 << "#                           |.....milliarcsec...........|                     |....degrees....|"
	 << endl;
    for (int iexp=0; iexp<exposures.size(); iexp++) {
      int iinst = exposures[iexp]->instrument;
      for (int iext=0; iext<instruments[iinst]->nExtensions; iext++) {
	string eName = instruments[iinst]->extensionNames.nameOf(iext);

	{
	  // Output for the fitted points:
	  const Accum& a = vvaccFit[iexp][iext];
	  double dx = 0., dy = 0., sigx=0., sigy=0.;
	  if (a.n>0) {
	    dx = a.sumxw / a.sumwx;
	    sigx = 1./sqrt(a.sumwx);
	    dy = a.sumyw / a.sumwy;
	    sigy = 1./sqrt(a.sumwy);
	  }
	  cout << "Fit     " << setw(3) << iexp
	       << " " << setw(3) << iext
	       << " " << a.summary()
	       << endl;
	}
	const Accum& a = vvaccReserve[iexp][iext];
	if (reserveFraction > 0. && a.n > 0) {
	  cout << "Reserve " << setw(3) << iexp
	       << " " << setw(3) << iext
	       << " " << a.summary()
	       << endl;
	}
      } // extension summary loop
    } // exposure summary loop
    
    // Output summary data for reference catalog
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
    // Get rid of Detections and Matches
    for (int ifield=0; ifield<fields.size(); ifield++)
      for (map<int, MCat>::iterator iaff = fields[ifield].affinities.begin();
	   iaff != fields[ifield].affinities.end();
	   ++iaff) 
	for (MSet::iterator im=iaff->second.begin();
	     im!=iaff->second.end();
	     ++im) {
	  (*im)->clear(true);
	  delete *im;
	}

    // Get rid of exposures
    for (int i=0; i<exposures.size(); i++)
      delete exposures[i];
    // Get rid of instruments
    for (int i=0; i<instruments.size(); i++)
      delete instruments[i];
    // ??? kill the orients of the Fields

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
