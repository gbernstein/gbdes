// Code to output WCS mapping of a grid of points on every CCD
#include <fstream>
#include <sstream>
#include <map>

#include "Std.h"
#include "StringStuff.h"
#include "Astrometry.h"
#include "PixelMapCollection.h"
#include "TemplateMap.h"
#include "PiecewiseMap.h"
#include "DECamInfo.h"
using namespace std;
using namespace astrometry;

string usage="TabulatePixmap: Map a grid of points on every device through a Pixel Map\n"
  "usage: TabulatePixmap <map name> <YAML file> \n"
  "      <map name>  the instrument or exposure name; device names will be appended\n"
  "      <YAML file> name of file holding serialized maps\n"
  "     stdout: each line is <device> <xpix> <ypix> <xworld> <yworld> for the map";

int
main(int argc, char *argv[])
{
  try {
    // Read parameters
    if (argc!=3) {
      cerr << usage << endl;
      exit(1);
    }

    string basename = argv[1];
    string filename = argv[2];

    // put all new kinds of PixelMaps we might be deserializing here:
    PixelMapCollection::registerMapType<astrometry::TemplateMap>();
    PixelMapCollection::registerMapType<astrometry::PiecewiseMap>();

    ifstream ifs(filename.c_str());
    if (!ifs) {
      cerr << "Could not open serialized WCS file " << filename << endl;
      exit(1);
    }
    PixelMapCollection pmc;
    if (!pmc.read(ifs)) {
      cerr << "File <" << filename << "> is not serialized WCS file" << endl;
      exit(1);
    }

    const double edge = 32.;  // Grid begins this far from device edge
    const double xSize = 2048.;
    const double ySize = 4096.;
    const int nX = 8; // Number of grid points on x side of each device
    const int nY = 2*nX;

    double dx = (xSize - 2*edge) / (nX-1);
    double dy = (ySize - 2*edge) / (nY-1);
    std::map<string, decam::Device> devices = decam::decamInfo();
    for (std::map<string, decam::Device>::const_iterator i = devices.begin();
	 i != devices.end();
	 ++i) {
      string devname = i->first;
      string mapName = basename + "/" + devname;
      if (pmc.mapExists(mapName)) {
	SubMap* sm = pmc.issueMap(mapName);
	for (int ix=0; ix<nX; ix++) {
	  double xpix = edge + ix*dx;
	  for (int iy=0; iy<nY; iy++) {
	    double ypix = edge + iy*dy;
	    double xw, yw;
	    sm->toWorld(xpix, ypix, xw, yw);
	    cout << devname
		 << fixed << setprecision(3)
		 << " " << setw(8) << xpix 
		 << " " << setw(8) << ypix
		 << setprecision(7)
		 << " " << setw(10) << xw 
		 << " " << setw(10) << yw
		 << endl;
	  }
	}
      }
    }
  } catch (std::runtime_error &m) {
    quit(m,1);
  }
}
