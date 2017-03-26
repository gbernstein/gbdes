// Spit out LDAC header
#include "FitsTable.h"
#include <vector>

using namespace std;

int main(int argc,
	 char *argv[])
{
  try {
    string filename = argv[1];
    int hdu = (argc>2 ? atoi(argv[2]) : 1);
    FITS::FitsTable readit(filename, FITS::ReadOnly, hdu);
    img::FTable ft = readit.extract();
    std::vector<std::string> vs;
    ft.readCell(vs, "Field Header Card",0);
    for (int i=0; i<vs.size(); i++) 
      cout << vs[i] << endl;
  /**
    img::Header h = ReadLdacHeader(filename, hdu);
    cerr << "Read done" << endl;
    cout << h; **/
  } catch (std::runtime_error &m) {
    quit(m,1);
  }
}
