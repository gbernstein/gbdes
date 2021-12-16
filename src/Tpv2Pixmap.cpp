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

string usage = "Tpv2Pixmap: convert ASCII FITS header TPV WCS to YAML serialized PixelMapCollection\n"
               "usage: Tpv2Pixmap <TPV in> <YAML out> <WCS name>\n"
               "   TPV in:   file holding ASCII FITS-header TPV WCS\n"
               "   YAML out: file to hold serialized PixelMapCollection\n"
               "   WCS name: name of the WCS in the output file.";

int main(int argc, char *argv[])
{
  try
  {
    // Read parameters
    if (argc != 4)
    {
      cerr << usage << endl;
      exit(1);
    }

    string tpvIn = argv[1];
    string yamlOut = argv[2];
    string wcsName = argv[3];

    // put all new kinds of PixelMaps we might be deserializing here:
    PixelMapCollection::registerMapType<astrometry::TemplateMap>();
    PixelMapCollection::registerMapType<astrometry::PiecewiseMap>();

    ifstream mapfs(tpvIn.c_str());
    if (!mapfs)
    {
      cerr << "Could not open TPV spec file " << tpvIn << endl;
      exit(1);
    }
    img::Header h;
    mapfs >> h;
    auto wcs = readTPV(h, wcsName);
    if (!wcs)
    {
      cerr << "Error reading input TPV map from file " << tpvIn << endl;
      exit(1);
    }

    PixelMapCollection pmc;
    pmc.learnWcs(*wcs);
    ofstream ofs(yamlOut.c_str());
    if (!ofs)
    {
      cerr << "Could not open output file " << yamlOut << endl;
      exit(1);
    }
    pmc.write(ofs);
    ofs.close();

    exit(0);
  }
  catch (std::runtime_error &m)
  {
    quit(m, 1);
  }
}
