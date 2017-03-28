// Test inversion of PixelMaps.
#include "Std.h"
#include "PixelMapCollection.h"
#include "PiecewiseMap.h"
#include "TemplateMap.h"
#include "PolyMap.h"

using namespace std;
using namespace astrometry;

string usage = "testmaps <pmc file> - test inversion of pixel maps. Give wcs file on cmd line,\n"
  "then follow prompts.";

int main(int argc,
	 char *argv[]) {
  if (argc!=2) {
    cerr << usage << endl;
    exit(1);
  }
  string pmcfile = argv[1];
  astrometry::PixelMapCollection::registerMapType<astrometry::TemplateMap>();
  astrometry::PixelMapCollection::registerMapType<astrometry::PiecewiseMap>();
  PixelMapCollection pmc;

  ifstream ifs(pmcfile.c_str());
  pmc.read(ifs);
  ifs.close();

  cerr << "Enter name of pixel map:";
  string mapname;
  while (cin >> mapname) {
    PixelMap* pm = pmc.cloneMap(mapname);
    if (!pm) {
      cerr << "Did not find map with name <" + mapname + ">; try again:" << endl;
      continue;
    }
    cerr << "Enter x, y pixel coordinates (-1 -1 to quit): ";
    double xpix,ypix;
    while (cin >> xpix >> ypix) {
      if (xpix==-1. && ypix==-1.)
	break;
      double xw, yw;
      pm->toWorld(xpix,ypix,xw,yw);
      cerr << "World coords: " << xw << " " << yw << endl;
      cerr << "Enter wc's to invert: ";
      cin >> xw >> yw;
      pm->toPix(xw,yw,xpix,ypix);
      cerr << "Pixel coords returned: " << xpix << " " << ypix << endl;
      double xw2,yw2;
      pm->toWorld(xpix,ypix,xw2,yw2);
      cerr << "which map to " << xw2 << " " << yw2 << endl;
      cerr << "with errors " << xw2-xw << " " << yw2-yw << endl;
      cerr << "Next x, y pixel coordinates (-1 -1 to quit): ";
    }
    delete pm;
    cerr << "Enter name of pixel map:";
  }
  exit(0);
}

