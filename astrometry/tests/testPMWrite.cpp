// Test writing and reading of PixelMaps

#include "PixelMap.h"
#include "yaml-cpp/yaml.h"
#include <sstream>

using namespace astrometry;
using namespace std;

int
main(int argc,
     char *argv[])
{
  SphericalICRS pole1(15.*DEGREE, -10.*DEGREE);
  SphericalICRS pole2(12.*DEGREE, -7.*DEGREE);
  Orientation orient1(pole1, 12.*DEGREE);
  Orientation orient2(pole2);

  ReprojectionMap mapIn(Gnomonic(orient1), 
			Gnomonic(orient2), DEGREE, "name1");

  string serial;
  {
    YAML::Emitter out;
    out << mapIn;
    serial = out.c_str();
    cout << "Serialized Map:" << endl;
    cout << serial << endl;
  }
  istringstream iss(serial);
  YAML::Node root = YAML::Load(iss);
  bool defaulted;
  PixelMap* mapOut = ReprojectionMap::create(root, defaulted, "name2");

  double xp = -0.25;
  double yp = 2.5;
  double xw, yw;
  mapIn.toWorld(xp,yp,xw,yw);
  cout << "From first map: " << xw << "," << yw << endl;
  mapOut->toWorld(xp,yp,xw,yw);
  cout << "From second map: " << xw << "," << yw << endl;
  delete mapOut;

  
}
