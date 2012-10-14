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


#include "Exposure.h"

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
  int extensionOrder;
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
    parameters.addMember("reFit",&reFit, def,
			 "Re-fit coordinate maps?", true);
    parameters.addMember("reserveFraction",&reserveFraction, def | low,
			 "Fraction of matches reserved from re-fit", 0., 0.);
    parameters.addMember("exposureOrder",&exposureOrder, def | low,
			 "Order of per-exposure map", 1, 0);
    parameters.addMember("extensionOrder",&extensionOrder, def | low,
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

    // Open up output catalog file if desired
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

    MCat::setDefaultMatchTolerance(matchtol);

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
      ff.readCells(names, "Name");
      ff.readCells(ra, "RA");
      ff.readCells(dec, "Dec");
      ff.readCells(fieldNumber, "fieldNumber");
      for (int i=0; i<names.size(); i++) {
	Exposure* expo = new Exposure(names[i],
				      Orientation(SphericalICRS(ra[i]*DEGREE, 
								dec[i]*DEGREE)));
	expo->fieldNumber = fieldNumber[i];
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
      //** want: exposure, instrument, device
      //** wcs (startpm)
      //** column keys for x, y, error, object number
      //** filename and extension of bintable
    }

    // A list of all Detections
    list<Detection> detections;
    // ...and a list of all Matches.
    list<Match*> matches;    // ??? Store objects, not pointers??

    
    // *** create astrometric models for each instrument and exposure
    // *** initialize instrument model by finding a reference exposure that has all devices
    // *** initialize exposure models

    // *** Read in all matched catalogs, create Detections, purge minMatches, map to world
    // *** summarize read-in statistics?

    // *** re-fit

    // ** write out put catalog
    // ** statistics out (separate code??)


    // Create a PixelMapCollection to hold all the components of the PixelMaps.
    PixelMapCollection mapCollection;

    // Create an identity map used by all reference catalogs
    PixelMapKey identityKey = mapCollection.add(new IdentityMap, "Identity");
    PixelMapChain nullChain;	// Can use empty chain to also be identity "map"
    PixelMapChain identityChain; 
    identityChain.append(identityKey);
    SubMap* identitySubMap = mapCollection.issue(identityChain);

	  f.orient = new Orientation(pole);
	  fields.push_back(f);
	  fieldNames.append(name);
	  Assert(fields.size()==fieldNames.size());
	}
	Field& fp = fields.back();
	ifield = fields.size()-1;

	// Open reference catalog for each
	ifstream refs(refcat.c_str());
	if (!refs) {
	  cerr << "Cannot open reference catalog " << refcat << endl;
	  exit(1);
	}

	cerr << "Reading reference catalog " << refcat
	     << " for field " << fp.name << endl;

	// Read each entry in reference catalog
	while (getlineNoComment(refs, buffer)) {
	  SphericalICRS posn;
	  double sigma;
	  string name;
	  istringstream iss(buffer);
	  // input RA in hours 
	  if (!(iss >> name >> posn >> sigma)) {
	    cerr << "Bad reference catalog entry <" << buffer << ">" << endl;
	    cerr << "In catalog " << refcat << endl;
	    exit(1);
	  }

	  Detection* d = new Detection;

	  // Project into this field's tangent plane
	  double xw, yw;
	  TangentPlane tp(posn, *fp.orient);
	  tp.getLonLat(xw, yw);
	  // set "pixel" and world coords
	  d->id = name;
	  d->xpix = xw/DEGREE;  
	  d->ypix = yw/DEGREE;
	  sigma = MAX(minReferenceError,sigma);
	  d->sigmaPix = sigma*ARCSEC/DEGREE; // Sigma comes in arcsec
	  d->xw = xw/DEGREE;
	  d->yw = yw/DEGREE;
	  d->clipsqx = pow(d->sigmaPix, -2.);
	  d->clipsqy = pow(d->sigmaPix, -2.);
	  // Reference objects can be weighted more/less than their
	  // uncertainties indicate:
	  d->wtx = d->clipsqx * referenceWeight;
	  d->wty = d->clipsqy * referenceWeight;
	  d->map = identitySubMap;
	  d->field = ifield;
	  d->exposure = REF_EXPOSURE;
	  d->extension = 0;	// means nothing here
	  d->isStar = true;	// *** note that all reference objects assumed stars
	  d->affinity = 0;	// *** and given affinity=0

#ifdef DUMP
	  cout << d->field
	       << " " << d->exposure
	       << " 0"
	       << " " << d->id
	       << " " << d->xw
	       << " " << d->yw
	       << endl;
#endif

	  // Match object into catalog for its affinity
	  // Map will create a new MCat if this is a new affinity
	  (fp.affinities[d->affinity]).add(d);
	} // end reference catalog object loop


	// Open tag catalog if there is one:
	if (!tagcat.empty()) {
	  ifstream tags(tagcat.c_str());
	  if (!tags) {
	    cerr << "Cannot open tag catalog " << tagcat << endl;
	    exit(1);
	  }

	  cerr << "Reading tag catalog " << tagcat
	       << " for field " << fp.name << endl;

	  // Read each entry in tag catalog
	  while (getlineNoComment(tags, buffer)) {
	    SphericalICRS posn;
	    double sigma;
	    string name;
	    istringstream iss(buffer);
	    // input RA in hours 
	    if (!(iss >> name >> posn >> sigma)) {
	      cerr << "Bad tag catalog entry <" << buffer << ">" << endl;
	      cerr << "In catalog " << refcat << endl;
	      exit(1);
	    }

	    // get ancillary info (tags!)
	    iss >> ws;
	    string ancillary;
	    getline(iss, ancillary);

	    Detection* d = new Detection;

	    // Project into this field's tangent plane
	    double xw, yw;
	    TangentPlane tp(posn, *fp.orient);
	    tp.getLonLat(xw, yw);
	    // set "pixel" and world coords
	    d->id = name;
	    d->xpix = xw/DEGREE;  
	    d->ypix = yw/DEGREE;
	    sigma = MAX(minReferenceError,sigma);
	    d->sigmaPix = sigma*ARCSEC/DEGREE; // Sigma comes in arcsec
	    d->xw = xw/DEGREE;
	    d->yw = yw/DEGREE;
	    d->clipsqx = pow(d->sigmaPix, -2.);
	    d->clipsqy = pow(d->sigmaPix, -2.);
	    // Tag catalog objects get ZERO WEIGHT:
	    d->wtx = 0.;
	    d->wty = 0.;
	    d->map = identitySubMap;
	    d->field = ifield;
	    d->exposure = TAG_EXPOSURE;	
	    d->extension = 0;	// means nothing here
	    d->isStar = true;	// *** note that all reference objects assumed stars
	    d->affinity = 0;	// *** and given affinity=0
	    d->ancillary = ancillary;
	    
	    // Match object into catalog for its affinity
	    // Map will create a new MCat if this is a new affinity
	    (fp.affinities[d->affinity]).add(d);
	  } // end tag catalog object loop
	} // end tag catalog read
      } // end field loop
    } // end use of fieldSpec file


    {
      // Read Exposures:
      ifstream ifs(exposureSpecs.c_str());
      if (!ifs) {
	cerr << "Can't open exposure specification file " << exposureSpecs << endl;
	exit(1);
      }

      string buffer;
      while (stringstuff::getlineNoComment(ifs, buffer)) {
	string basename;
	string field;
	string instName;
	string filter;
	SphericalICRS pointing;
	double positionAngle=0.; // ???
	{
	  istringstream iss(buffer);
	  if (!(iss >> basename >> pointing 
		>> field >> filter >> instName)) {
	    cerr << "Bad exposure spec <" << buffer << ">" << endl;
	    exit(1);
	  }
	}
	  
	// Check for duplicate exposure name
	if (exposureNames.has(basename)) {
	  cerr << "Duplicate exposure basename: " << basename << endl;
	  exit(1);
	}

	// Find the field that this exposure is in
	int ifield = fieldNames.indexOf(field);
	if (ifield < 0) {
	  cerr << "Exposure " << basename << " has field name <" << field 
	       << "> not given in field specification file" << endl;
	  exit(1);
	}
	Field& fp = fields[ifield];

	// Add this filter name to the index list if needed
	if (filterNames.indexOf(filter) < 0) filterNames.append(filter);
	// Choose affinity numbers for stars and galaxies in this exposure
	int starAffinity = 0;
	// If useAffinities=true, then galaxies are only matched if filters match:
	int galaxyAffinity = useAffinities ? filterNames.indexOf(filter)+1 : 0;

	// Find the instrument that is in use
	int iinst = instrumentNames.indexOf(instName);
	bool firstExposureForInstrument = false;
	if (iinst < 0) {
	  // First exposure with this instrument
	  firstExposureForInstrument = true;
	  instruments.push_back(new Instrument(instName));
	  instrumentNames.append(instName);
	  iinst = instrumentNames.indexOf(instName);
	  Assert(instruments.size()==instrumentNames.size());
	  Assert(iinst>=0);
	}
	Instrument& inst = *instruments[iinst];

	// Add this exposure to the list of ID's
	exposureNames.append(basename);
	int iexp = exposureNames.indexOf(basename);
	Exposure* expo = new Exposure;
	expo->name = basename;
	expo->field = ifield;
	expo->instrument = iinst;
	// Exposure will delete this Orientation when it is destroyed.
	// Must stick around for all the TangentPlane coordinates to work!
	expo->orient = new Orientation(pointing, positionAngle);

	exposures.push_back(expo);
	Assert(exposures.size()==exposureNames.size());

	// Tell each field and instrument about all of its exposures
	fp.exposures.push_back(iexp);
	inst.exposures.push_back(iexp);

	// Read catalog
	string catfile = exposures.back()->name + catalogSuffix;
	cerr << "Reading catalog " << exposures.size()-1 << " " << catfile << endl;
	ifstream catfs(catfile.c_str());
	if (!catfs) {
	  cerr << "Could not open catalog file " << catfile << endl;
	  exit(1);
	}
	while (getlineNoComment(catfs, buffer)) {
	  string name;
	  string extName;
	  double xpix;
	  double ypix;
	  double sigmaPix;
	  istringstream iss(buffer);
	  if (!(iss >> name >> extName >> xpix >> ypix >> sigmaPix)) {
	    cerr << "Bad catalog entry <" << buffer
		 << "> in catalog " << catfile
		 << endl;
	    exit(1);
	  }
	  // See if there is a s/g specifier
	  iss >> ws;	// strip leading whitespace
	  string ancillary;
	  getline(iss, ancillary);
	  bool isStar = true;   // Set true if not using affinities
	  if (ancillary.empty()) {
	    if (useAffinities) {
	      // Need a s/g spec, don't have it
	      cerr << "Using affinities but catalog entry does not give s|g"
		" in catalog " << catfile << endl;
	      cerr << buffer << endl;
	      exit(1);
	    }
	  } else if (ancillary[0]=='s') {
	    // It's a star
	    isStar = true;
	    // Strip first character and additional whitespace
	    ancillary.erase(0,1);
	    while (!ancillary.empty() && std::isspace(ancillary[0])) ancillary.erase(0,1);
	  } else if (ancillary[0]=='g') {
	    // It's a galaxy
	    isStar = false;
	    // Strip first character and additional whitespace
	    ancillary.erase(0,1);
	    while (!ancillary.empty() && std::isspace(ancillary[0])) ancillary.erase(0,1);
	  // and the rest of the string is "ancillary"
	  } else if (useAffinities) {
	    // Did not get s or g but needed it for affinities:
	    cerr << "Using affinities but catalog entry does not give s|g"
	      " in catalog " << catfile << endl;
	    cerr << buffer << endl;
	    exit(1);
	  }

	  // Done reading the detection's fields. 

	  int iext = inst.extensionNames.indexOf(extName);
	  if (iext < 0) {
	    // New extension for this instrument.
	    inst.addExtension(extName);
	    iext = inst.extensionNames.indexOf(extName);
	  }

	  // First pad the Exposure's startpm array to reach this extension number:
	  while (expo->startpm.size() < iext+1) expo->startpm.push_back(0);


	  // Get the WCS map if this is first Detection in this extension
	  // for this Exposure
	  if (!expo->startpm[iext]) {
	    string mapfile = basename + "_" + extName + oldHeadSuffix;
	    ifstream mapfs(mapfile.c_str());
	    if (!mapfs) {
	      cerr << "Could not open map spec file " << mapfile << endl;
	      exit(1);
	    } else {
	      cerr << "Opened map spec file " << mapfile << endl;
	    }

	    img::ImageHeader h = img::HeaderFromStream(mapfs);
	    expo->startpm[iext] = new SCAMPMap(h, fp.orient);
	    if (expo->startpm[iext]==0) {
	      cerr << "Error reading initial WCS map for extension " << iext
		   << " in file " << mapfile
		   << endl;
	      exit(1);
	    }
	  }

	  // Don't use things with too-large errors
	  if (sigmaPix > maxPixError) continue;

	  Detection* d = new Detection;
	  d->id = name;
	  d->xpix = xpix;
	  d->ypix = ypix;
	  d->sigmaPix = MAX(minPixError,sigmaPix);
	  d->field = ifield;
	  d->exposure = iexp;
	  d->extension = iext;
	  d->affinity = d->isStar ? starAffinity : galaxyAffinity;
	  d->isStar = isStar;
	  d->ancillary = ancillary;

	  //   map pixel coordinates to world from original map
	  double xw, yw;
	  exposures[iexp]->startpm[iext]->toWorld(xpix, ypix, xw, yw);
	  //**/cerr << iexp << " " << iext << " " << xw << " " << yw << endl;
	  d->xw = xw;
	  d->yw = yw;
	  Matrix22 dwdp = exposures[iexp]->startpm[iext]->dWorlddPix(xpix, ypix);
	  d->wtx = pow(d->sigmaPix,-2.) / (dwdp(0,0)*dwdp(0,0)+dwdp(0,1)*dwdp(0,1));
	  d->wty = pow(d->sigmaPix,-2.) / (dwdp(1,0)*dwdp(1,0)+dwdp(1,1)*dwdp(1,1));
	  d->clipsqx = d->wtx;
	  d->clipsqy = d->wty;

#ifdef DUMP
	  cout << d->field
		   << " " << d->exposure
		   << " " << d->extension
		   << " " << d->id
		   << " " << d->xw
		   << " " << d->yw
		   << endl;
#endif

	  // Tell Instrument about this detection to complete its range
	  inst.addDetection(*d);

	  // Match object into its catalog of its affinity
	  (fp.affinities[d->affinity]).add(d);
	
	} // end detection loop
      } // end exposure loop
    }// close exposure specification file
    
    /**/cerr << "done reading exposures" << endl;
    for (int ifield=0; ifield<fields.size(); ifield++) {
      for (map<int, MCat>::iterator iaff = fields[ifield].affinities.begin();
	   iaff != fields[ifield].affinities.end();
	   ++iaff) {
	MCat& mc = iaff->second;
	int affinity = iaff->first;
	// Purge Matches with multiple entries from single exposure or <minMatches:
	mc.minMatches(minMatches, true);
	cout << "Field " << fields[ifield].name
	     << " affinity " << affinity
	     << " has " << mc.size()
	     << " matches." 
	     << endl;
	long int mcount=0;
	int dof=0;
	double chi=0.;
	double maxdev=0.;
	for (MCat::const_iterator i=mc.begin();
	     i != mc.end();
	     ++i) {
	  mcount += (*i)->size();
	  chi += (*i)->chisq(dof, maxdev);
	}
	cout << " with " << mcount << " total detections." << endl;
	cout << " chisq " << chi << " / " << dof << " dof maxdev " << sqrt(maxdev) << endl;
	mc.purgeSelfMatches(true);
	cout << "After purgeSelfMatches: " << mc.size()
	     << " matches." 
	     << endl;
	mcount=0;
	dof = 0;
	chi = maxdev = 0.;
	for (MCat::const_iterator i=mc.begin();
	     i != mc.end();
	     ++i) {
	  mcount += (*i)->size();
	  chi += (*i)->chisq(dof, maxdev);
	}
	cout << " with " << mcount << " total detections." << endl;
	cout << " chisq " << chi << " / " << dof << " dof maxdev " << sqrt(maxdev) << endl;
      }
    }

    // Now everyone is matched up and cleaned up.

    ///////////////////////////////////////////////////////////
    // Now do the re-fitting if desired
    ///////////////////////////////////////////////////////////
    if (reFit) {


      // Set up coordinate maps for each instrument & exposure
      for (int iinst=0; iinst<instruments.size(); iinst++) {
	Instrument* inst = instruments[iinst];
	// First exposure for this instrument gets identity exposure map,
	// Used to initialize the map for each extension
	Exposure* exp1 = exposures[inst->exposures.front()];
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
	const int nGridPoints=400;	// Number of test points for map initialization

	for (int iext=0; iext<inst->nExtensions; iext++) {
	  Bounds<double> b=inst->domains[iext];
	  double step = sqrt(b.area()/nGridPoints);
	  int nx = static_cast<int> (ceil((b.getXMax()-b.getXMin())/step));
	  int ny = static_cast<int> (ceil((b.getYMax()-b.getYMin())/step));
	  double xstep = (b.getXMax()-b.getXMin())/nx;
	  double ystep = (b.getYMax()-b.getYMin())/ny;
	  list<Detection*> testPoints;
	  // starting pixel map for the first exposure with this instrument
	  PixelMap* startpm=exp1->startpm[iext];
	  if (!startpm) {
	    cerr << "Did not read starting WCS for extension " << inst->extensionNames.nameOf(iext)
		 << " in first exposure " << exp1->name
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
	  PixelMap* pm = new PolyMap(extensionOrder, worldTolerance);
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
	  PixelMapKey key = mapCollection.add(pm, inst->name + "_" + inst->extensionNames.nameOf(iext));
	  inst->maps[iext].append(key);
	  inst->pixelMaps[iext] = mapCollection.issue(inst->maps[iext]);
	} // extension loop

	int pointsPerExtension = nGridPoints / inst->nExtensions + 1;
	for (int iexp=1; iexp<inst->exposures.size(); iexp++) {
	  // Each additional exposure needs a warping tranformation
	  Exposure* expo = exposures[inst->exposures[iexp]];
	  PixelMap* warp;
	  if (exposureOrder==1) 
	    warp = new LinearMap;
	  else
	    warp = new PolyMap(exposureOrder, worldTolerance);
	  // Set numerical-derivative step to 1 arcsec, polynomial is working in degrees:
	  warp->setPixelStep(ARCSEC/DEGREE);

	  // And a reprojection map:
	  const Field& fp = fields[expo->field];
	  expo->reproject = mapCollection.add(new ReprojectionMap( TangentPlane(expo->orient->getPole(),
									      *expo->orient),
								   TangentPlane(fp.orient->getPole(),
										*fp.orient),
								   DEGREE),
					      expo->name + "_reprojection");
	  const PixelMap* reproject = mapCollection.element(expo->reproject);

	  // Build new test points
	  list<Detection*> testPoints;
	  for (int iext=0; iext<inst->nExtensions; iext++) {
	    // Skip extensions that did not have data in this exposure,
	    // hence did not read their initial pixel maps:
	    if (iext >= expo->startpm.size() || !expo->startpm[iext]) continue;
	    Bounds<double> b=inst->domains[iext];
	    double step = sqrt(b.area()/pointsPerExtension);
	    int nx = static_cast<int> (ceil((b.getXMax()-b.getXMin())/step));
	    int ny = static_cast<int> (ceil((b.getYMax()-b.getYMin())/step));
	    double xstep = (b.getXMax()-b.getXMin())/nx;
	    double ystep = (b.getYMax()-b.getYMin())/ny;
	    // starting pixel map for the first exposure with this instrument
	    PixelMap* startpm=expo->startpm[iext];
	    for (int ix=0; ix<=nx; ix++) {
	      for (int iy=0; iy<=ny; iy++) {
		double xpix = b.getXMin() + ix*xstep;
		double ypix = b.getYMin() + iy*ystep;
		double xw, yw;
		startpm->toWorld(xpix, ypix, xw, yw);
		// Undo the change of projections:
		reproject->toPix(xw,yw,xw,yw);
		// And process the pixel coordinates through Instrument map:
		inst->pixelMaps[iext]->toWorld(xpix,ypix,xpix,ypix);
		Detection* d = new Detection;
		d->xpix = xpix;
		d->ypix = ypix;
		d->xw = xw;
		d->yw = yw;
		testPoints.push_back(d);
	      }
	    }
	  } // finish extension loop for test points

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
	}

	// Now construct compound maps for every extension of
	// every exposure and save SubMap for it
	for (int iexp=0; iexp<inst->exposures.size(); iexp++) {
	  Exposure* expo = exposures[inst->exposures[iexp]];
	  for (int iext=0; iext<instruments[iinst]->nExtensions; iext++) {
	    // Start with the instrument map for extension:
	    PixelMapChain chain = inst->maps[iext];
	    // Then the exposure's warp
	    chain.append(expo->warp);
	    // And exposure's reprojection:
	    chain.append(expo->reproject);
	    expo->extensionMaps.push_back(mapCollection.issue(chain));

	    Assert(expo->extensionMaps.size()==iext+1);

	  } // extension loop for compound maps
	} // exposure loop for compound maps
      } // instrument loop

      /**/cerr << "Total number of map elements " << mapCollection.nMaps()
	       << " with " << mapCollection.nParams() << " free parameters."
	       << endl;

      // Assign maps to every detection
      for (int ifield=0; ifield<fields.size(); ifield++) {
	for (map<int, MCat>::iterator iaff = fields[ifield].affinities.begin();
	     iaff != fields[ifield].affinities.end();
	     ++iaff) {
	  MCat& mc = iaff->second;
	  for (MCat::iterator icat=mc.begin();
	       icat!=mc.end();
	       ++icat) {
	    Match& m = *(*icat);
	    for (Match::const_iterator idet=m.begin();
		 idet != m.end();
		 ++idet) {
	      Detection* d = *idet;
	      if (d->exposure < 0) 
		d->map = identitySubMap;
	      else 
		d->map = exposures[d->exposure]->extensionMaps[d->extension];
	    }
	  }
	}
      }

      // Mark a selected fraction of all matches as being reserved from fit
      if (reserveFraction > 0.) {
	ran::UniformDeviate u;
	for (int ifield=0; ifield<fields.size(); ifield++) {
	  for (map<int, MCat>::iterator iaff = fields[ifield].affinities.begin();
	       iaff != fields[ifield].affinities.end();
	       ++iaff) {
	    MCat& mc = iaff->second;
	    for (MCat::iterator icat=mc.begin();
		 icat!=mc.end();
		 ++icat) {
	      Match& m = *(*icat);
	      m.setReserved( u < reserveFraction );
	    }
	  }
	}
      } // Done reserving matches

      // make CoordAlign class
      CoordAlign ca(mapCollection);
      for (int ifield=0; ifield<fields.size(); ifield++)
	for (map<int, MCat>::iterator iaff = fields[ifield].affinities.begin();
	     iaff != fields[ifield].affinities.end();
	     ++iaff) 
	  ca.add(iaff->second);
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
      for (int ifield=0; ifield<fields.size(); ifield++)
{
      for (map<int, MCat>::iterator iaff = fields[ifield].affinities.begin();
           iaff != fields[ifield].affinities.end();
           ++iaff) {
        MCat& mc = iaff->second;
          for (int iext=0; iext<instruments[0]->nExtensions; iext++) {
             mc.detailed_info(iext);
          }
          }
         
	} }
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

    } else { 
      // reFit = false.  Do sigma-clipping on input world coordinates.
      // Create a CoordAlign class used just for clipping
      PixelMapCollection nullCollection;	// an empty list
      CoordAlign ca(nullCollection);
      for (int ifield=0; ifield<fields.size(); ifield++)
	for (map<int, MCat>::iterator iaff = fields[ifield].affinities.begin();
	     iaff != fields[ifield].affinities.end();
	     ++iaff) 
	  ca.add(iaff->second);
      int nclip;
      double oldthresh=0.;
      do {
	// Report number of active Matches / Detections in each iteration:
	long int mcount=0;
	long int dcount=0;
	ca.count(mcount, dcount, false, 2);
	double max;
	int dof=0;
	double chisq= ca.chisqDOF(dof, max, false);
	cout << "Clipping " << mcount << " matches with " << dcount << " detections "
	     << " chisq " << chisq << " / " << dof << " dof,  maxdev " << sqrt(max) 
	     << " sigma" << endl;

	double thresh = sqrt(chisq/dof) * clipThresh;
	cout << "  new clip threshold: " << thresh << " sigma"
	     << endl;
	if (thresh >= max) break;
	if (oldthresh>0. && (1-thresh/oldthresh)<minimumImprovement) break;
	oldthresh = thresh;
	nclip = ca.sigmaClip(thresh, false, clipEntireMatch);
	cout << "Clipped " << nclip
	     << " matches " << endl;
      } while (nclip>0);
    }

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
