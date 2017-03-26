// Dump the DECam corners and centers in field coordinates
#include "DECamInfo.h"
using namespace decam;

int main(int argc,
	 char *argv[])
{
  typedef std::map<string,Device> Dmap;
  Dmap d = decamInfo();

  for (Dmap::const_iterator i = d.begin();
       i != d.end();
       ++i) {
    cout << i->first << " " << i->second.b << endl;
  }
  exit(0);
}
