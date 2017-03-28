// Test writing and reading of PixelMaps

#include <sstream>
#include "TPVMap.h"
#include "PixelMap.h"
#include "PixelMapCollection.h"

using namespace astrometry;
using namespace std;

int
main(int argc,
     char *argv[])
{
  try {
    ifstream mapfs(argv[1]);
    if (!mapfs) {
      cerr << "Could not open map spec file " << argv[1] << endl;
      exit(1);
    }
    img::Header h; 
    mapfs >> h;
    Wcs* wcsIn = readTPV(h, "tryMe");
   
    PixelMapCollection pmc;
    pmc.learnWcs(*wcsIn);
    ostringstream oss;
    pmc.write(oss);
    cout << oss.str();
    cout << "---------------------" << endl;
    istringstream iss(oss.str());
    pmc.read(iss, "Dup_");
    Wcs* dup = pmc.issueWcs("Dup_tryMe");
    double xp, yp;
    cerr << "Enter xp, yp: ";
    while (cin >> xp >> yp) {
      cout << wcsIn->toSky(xp, yp) << endl;
      cout << dup->toSky(xp, yp) << endl;
    }

  } catch (std::runtime_error& m) {
    quit(m,1);
  }
}

