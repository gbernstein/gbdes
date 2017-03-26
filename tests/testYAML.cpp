// Program to test use of YAML library
#include <vector>
#include <fstream>
#include "PhotoMapCollection.h"

using namespace std;
using namespace photometry;

int
main(int argc,
     char *argv[])
{
  PhotoMapCollection pmc;

  ConstantMap pm1(7., "yo");
  pmc.learnMap(pm1);

  PolyMap pm2(3, PolyMap::Device, "mama");
  pmc.learnMap(pm2);

  ofstream ofs("test.yaml");
  pmc.write(ofs, "Creating test serialization");
  ofs.close();

  ifstream ifs("test.yaml");
  PhotoMapCollection pmc2;
  pmc2.read(ifs);

  exit(0);

}
