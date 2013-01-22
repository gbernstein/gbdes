// 	$Id: ApplyWCS.cpp,v 1.3 2012/01/27 01:38:17 garyb Exp $	
#include <fstream>
#include <sstream>
#include <map>

#include "Std.h"
#include "Astrometry.h"
#include "TPVMap.h"
#include "Header.h"
#include "StringStuff.h"

using namespace std;
using namespace astrometry;

string usage="ApplyWCS: use WCS systems specified by FITS keywords to map from one\n"
             "   pixel system to another, or from pixel system to/from world coords\n"
             "usage: ApplyWCS <WCS in> <WCS out> [ra dec]\n"
	     "      <WCS in, out>:  each is either (a) name of file holding FITS-style WCS spec,\n"
             "             (b) the word \"gnomonic\" followed by RA and Dec (in hours, degrees)\n"
             "                 that will be pole of a projection to/from (xi,eta) in degrees\n"
             "             (c) a dash - to use RA and Dec directly (in hours, degrees).\n"
             "     stdin: each line is <xin> <yin> which are either input pixel or input WCS coords\n"
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
    Wcs* wcsIn = 0;
    int iarg=1;
    string wcsfilein = argv[iarg];
    iarg++;
    SphericalCoords* inProjection = 0;

    // Choose input method:
    if (stringstuff::nocaseEqual(wcsfilein, "gnomonic")) {
      if (argc < iarg+3) {
	cerr << "Not enough command line arguments after gnomonic input" << endl;
	exit(1);
      }
      string radecString = argv[iarg];
      radecString += " ";
      radecString += argv[iarg+1];
      iarg += 2;
      SphericalICRS pole;
      istringstream iss(radecString);
      if (!(iss >> pole)) {
	cerr << "Bad RA/Dec: " << radecString << endl;
	exit(1);
      }
      Orientation orient(pole);
      inProjection = new Gnomonic(orient);
    } else if (wcsfilein=="-") {
      // Nothing to do, will just read RA and Dec
    } else {
      // Get WCS from file:
      ifstream mapfs(wcsfilein.c_str());
      if (!mapfs) {
	cerr << "Could not open map spec file " << wcsfilein << endl;
	exit(1);
      }
      img::Header h; 
      mapfs >> h;
      wcsIn = readTPV(h);
      if (!wcsIn) {
	cerr << "Error reading input WCS map from file " << wcsfilein
	     << endl;
	exit(1);
      }
    }

    if (argc < iarg+1) {
      cerr << "Not enough command line arguments" << endl;
      exit(1);
    }
    Wcs* wcsOut = 0;
    string wcsfileout = argv[iarg];
    iarg++;
    SphericalCoords* outProjection = 0;

    // Choose output method:
    if (stringstuff::nocaseEqual(wcsfileout, "gnomonic")) {
      if (argc < iarg+2) {
	cerr << "Not enough command line arguments after gnomonic input" << endl;
	exit(1);
      }
      string radecString = argv[iarg];
      radecString += " ";
      radecString += argv[iarg+1];
      iarg += 2;
      SphericalICRS pole;
      istringstream iss(radecString);
      if (!(iss >> pole)) {
	cerr << "Bad RA/Dec: " << radecString << endl;
	exit(1);
      }
      Orientation orient(pole);
      outProjection = new Gnomonic(orient);
    } else if (wcsfileout=="-") {
      // Nothing to do, will just read RA and Dec
    } else {
      // Get WCS from file:
      ifstream mapfs(wcsfileout.c_str());
      if (!mapfs) {
	cerr << "Could not open map spec file " << wcsfileout << endl;
	exit(1);
      }
      img::Header h; 
      mapfs >> h;
      wcsOut = readTPV(h);
      if (!wcsOut) {
	cerr << "Error reading output WCS map from file " << wcsfileout
	     << endl;
	exit(1);
      }
    }

    // Begin reading data
    string buffer;
    while (stringstuff::getlineNoComment(cin, buffer)) {
      double xpix;
      double ypix;
      istringstream iss(buffer);
      if (!(iss >> xpix >> ypix)) {
	cerr << "Bad input entry <" << buffer
	     << endl;
	exit(1);
      }
      // Input coordinates to world system:
      SphericalICRS icrs;
      if (wcsIn)
	icrs.convertFrom(wcsIn->toSky(xpix, ypix));
      else if (inProjection) {
	inProjection->setLonLat(xpix*DEGREE, ypix*DEGREE);
	icrs.convertFrom(*inProjection);
      } else
	icrs.setLonLat(xpix*DEGREE/15., ypix*DEGREE);

      // Now choose output quantities
      if (wcsOut) {
	xpix = ypix = 0.;
	wcsOut->fromSky(icrs, xpix, ypix);
      } else if (outProjection) {
	outProjection->convertFrom(icrs);
	outProjection->getLonLat(xpix, ypix);
	xpix /= DEGREE; ypix /= DEGREE;
      } else {
	icrs.getLonLat(xpix, ypix);
	xpix /= DEGREE*15.;
	ypix /= DEGREE;
      }
      cout << xpix << " " << ypix << endl;

      string appendix;
      while(iss >> appendix) {
          cout << " " << appendix;
      }
      cout << endl;
    }
    if (wcsIn) delete wcsIn;
    if (wcsOut) delete wcsOut;
    if (inProjection) delete inProjection;
    if (outProjection) delete outProjection;
  } catch (std::runtime_error &m) {
    quit(m,1);
  }
}
