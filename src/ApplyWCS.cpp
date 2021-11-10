#include <fstream>
#include <sstream>
#include <map>

#include "Std.h"
#include "StringStuff.h"
#include "Astrometry.h"
#include "PixelMapCollection.h"
#include "TemplateMap.h"
#include "PiecewiseMap.h"
#include "Header.h"
#include "TPVMap.h"
using namespace std;
using namespace astrometry;

string usage="ApplyWCS: use WCS systems specified by YAML file or FITS keywords to map from one\n"
             "   pixel system to another, or from pixel system to/from world coords\n"
             "usage: ApplyWCS <WCS in> <WCS out> \n"
	     "      <WCS in, out>:  each is either:\n"
             "             (a) name of ASCII file holding FITS-style WCS spec,\n"
             "             (b) <wcsname>@<filename>, reads WCS with given name serialized in the file\n"
             "             (c) the word \"gnomonic\" followed by RA and Dec (in hours, degrees)\n"
             "                 that will be pole of a projection to/from (xi,eta) in degrees\n"
             "             (d) either \"icrs\" or \"ecliptic\" to use true sky coords in degrees\n"
             "             (e) a dash - to get input/output in the native projection of the Wcs.\n"
             "     stdin: each line is <xin> <yin> [color] which are either input pixel \n"
             "            or input WCS coords.  Optional color defaults to zero.\n"
             "     stdout: each line is <xout> <yout> for either output pixel system or WCS";

int
main(int argc, char *argv[])
{
  try {
    // Read parameters
    if (argc<3) {
      cerr << usage << endl;
      exit(1);
    }

    // put all new kinds of PixelMaps we might be deserializing here:
    PixelMapCollection::registerMapType<astrometry::TemplateMap>();
    PixelMapCollection::registerMapType<astrometry::PiecewiseMap>();

    int iarg=1;
    vector<Wcs*> wcs(2, (Wcs*) 0);
    vector<unique_ptr<SphericalCoords>> projection(2);
    vector<bool> useNative(2,false);

    // Save pmc to avoid reading it twice if it's repeated
    PixelMapCollection* pmc = 0;
    string pmcFile = "";
    
    for (int i=0; i<2; i++) {
      // Choose input (i=0) / output (i=1) system:
      if (iarg >= argc) {
	cerr << "Not enough command line arguments."  << endl;
	cerr << usage << endl;
	exit(1);
      }

      string wcsfile = argv[iarg++];

      if (stringstuff::nocaseEqual(wcsfile, "gnomonic")) {
	if (iarg+1 >= argc) {
	  cerr << "Not enough command line arguments after gnomonic input" << endl;
	  exit(1);
	}
	string radecString = argv[iarg++];
	radecString += " ";
	radecString += argv[iarg++];
	SphericalICRS pole;
	istringstream iss(radecString);
	if (!(iss >> pole)) {
	  cerr << "Bad RA/Dec: " << radecString << endl;
	  exit(1);
	}
	Orientation orient(pole);
	projection[i] = unique_ptr<SphericalCoords>(new Gnomonic(orient));
      } else if (stringstuff::nocaseEqual(wcsfile, "icrs")) {
	projection[i] = unique_ptr<SphericalCoords>(new SphericalICRS);
      } else if (stringstuff::nocaseEqual(wcsfile, "ecliptic")) {
	projection[i] = unique_ptr<SphericalCoords>(new SphericalEcliptic);
      } else if (wcsfile=="-") {
	useNative[i] = true;
      } else if (wcsfile.find('@')!=string::npos) {
	// Deserialize WCS from a file
	size_t at = wcsfile.find('@');
	string wcsname = wcsfile.substr(0,at);
	at++;
	if (at >= wcsfile.size()) {
	  cerr << "Missing filename in wcsname@filename argument: " << wcsfile << endl;
	  exit(1);
	}
	string filename = wcsfile.substr(at);
	if (pmc==0 || filename!=pmcFile)  {
	  if (pmc) delete pmc;
	  ifstream ifs(filename.c_str());
	  if (!ifs) {
	    cerr << "Could not open serialized WCS file " << filename << endl;
	    exit(1);
	  }
	  pmc = new PixelMapCollection;
	  if (!pmc->read(ifs)) {
	    cerr << "File <" << filename << "> is not serialized WCS file" << endl;
	    exit(1);
	  }
	  pmcFile = filename;
	}
	wcs[i] = pmc->cloneWcs(wcsname);
      } else {
	// Get WCS as TPV in FITS-header-style file
	ifstream mapfs(wcsfile.c_str());
	if (!mapfs) {
	  cerr << "Could not open TPV spec file " << wcsfile << endl;
	  exit(1);
	}
	img::Header h; 
	mapfs >> h;
	wcs[i] = readTPV(h);
	if (!wcs[i]) {
	  cerr << "Error reading input TPV map from file " << wcsfile << endl;
	  exit(1);
	}
      }
    }

    if (pmc) delete pmc;

    if ( useNative[0] && !wcs[1] || useNative[1] && !wcs[0]) {
      cerr << "If you ask for native projection with \"-\" argument, other arg must be a WCS" << endl;
      exit(1);
    }

    // If we are outputting degrees, use 7 decimal places.
    if (projection[1] || useNative[1])
      cout << fixed << setprecision(7);

    // Begin reading data
    string buffer;
    while (stringstuff::getlineNoComment(cin, buffer)) {
      double xin;
      double yin;
      double color;
      double xout = 0.;
      double yout = 0.;
      istringstream iss(buffer);
      if (!(iss >> xin >> yin)) {
	cerr << "Bad input entry <" << buffer
	     << endl;
	exit(1);
      }
      if (!(iss >> color))
	color = 0.;
      // Input coordinates to world system:
      SphericalICRS icrs;
      if (wcs[0]) {
	if (useNative[1]) {
	  // Just map input to native system of the input Wcs
	  wcs[0]->getMap()->toWorld(xin, yin, xout, yout, color);
	  cout << xout << " " << yout << endl;
	  continue;
	} else {
	  icrs.convertFrom(wcs[0]->toSky(xin, yin, color));
	}
      } else if (projection[0]) {
	projection[0]->setLonLat(xin*DEGREE, yin*DEGREE);
	icrs.convertFrom(*projection[0]);
      } else if (useNative[0]) {
	// go straight to output system's pixels
	wcs[1]->getMap()->toPix(xin, yin, xout, yout, color);
	cout << xout << " " << yout << endl;
	continue;
      } else {
	cerr << "Logic error." << endl;
	exit(1);
      }

      // Now convert to output
      if (wcs[1]) {
	wcs[1]->fromSky(icrs, xout, yout, color);
      } else if (projection[1]) {
	projection[1]->convertFrom(icrs);
	projection[1]->getLonLat(xout, yout);
	xout /= DEGREE; yout /= DEGREE;
      } else {
	cerr << "Logic error." << endl;
	exit(1);
      }
      cout << xout << " " << yout << endl;
    }

    // Clean up
    for (int i=0; i<2; i++) {
      if (wcs[i]) delete wcs[i];
    }
  } catch (std::runtime_error &m) {
    quit(m,1);
  }
}
