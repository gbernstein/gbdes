// Test writing and reading of Spherical coordinates to YAML.
#include "Astrometry.h"
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
  Orientation orient(pole1, 12.*DEGREE);
  Gnomonic gn(pole2, orient);

  cout << "xi,eta: " << gn << endl;
  
  YAML::Emitter out;

  out << YAML::BeginMap
      << YAML::Key << "c1" << YAML::Value << pole1
      << YAML::Key << "c2" << YAML::Value << gn
      << YAML::EndMap;

  cout << "YAML produced------- " << endl;
  cout << out.c_str() << endl;
  cout << "--------" << endl;
  

  // Read back in
  YAML::Node root = YAML::Load(out.c_str());
  SphericalCoords* in1 = SphericalCoords::deserialize(root["c1"]);
  cout << "ICRS coords: " <<  *in1 << endl;
  SphericalCoords* in2 = SphericalCoords::deserialize(root["c2"]);
  cout << "xi, eta: " << *in2 << endl;
  cout << "Back to ICRS: " << SphericalICRS(*in2) << endl;
  delete in1, in2;
  
  
}
