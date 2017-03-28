// test the FTable stuff.
#include <iostream>
#include <vector>
#include "Std.h"

#include "FITSTable.h"

using namespace std;
using namespace FITS;

int
main(int argc,
     char *argv[])
{
  try {
    string filename = argv[1];
    int hdunumber = atoi(argv[2]);
    string rowfilter = argv[3];
    string colfilter = argv[4];
    vector<string> cn;


    FitsTable h(filename, ReadOnly, hdunumber);
    cout << "Extension name " << h.getName() << endl;
    cout << " type " << h.getType() << endl;
    cout << " number " << h.getNumber() << endl;
    cout << " Header dump:" << endl;
    cout << *(h.header()) << endl;

    FitsTable ftnew("junk.fits", ReadWrite + OverwriteFile, "Extract");
    img::FTable te = h.extract(rowfilter, vector<string>(1,colfilter));
    ftnew.adopt(te);
  } catch (std::runtime_error &m) {
    quit(m,1);
  }
}

