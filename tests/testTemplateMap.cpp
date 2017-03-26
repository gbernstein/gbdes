// Test writing and reading of TemplateMaps

#include "TemplateMap.h"
#include <sstream>

using namespace astrometry;
using namespace std;

int
main(int argc,
     char *argv[])
{
  string mapname = argv[1];
  string filename = argv[2];
  TemplateMap tm(1024., mapname+"_low", mapname + "_high", filename, "blah");

  double xp, yp;
  while (cin >> xp >> yp) {
    double xw, yw;
    tm.toWorld(xp,yp,xw,yw);
    cout << xw << " " << yw << endl;
    tm.toPix(xw, yw, xp, yp);
    cout << "... " << xp << " " << yp << endl;
  }

  YAML::Emitter out;
  tm.write(out);
  cout << out.c_str() << endl;
}
