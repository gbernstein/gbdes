// test the FTable stuff.
#include <iostream>
#include <vector>
#include "Std.h"

#include "Hdu.h"

using namespace std;
using namespace FITS;

int
main(int argc,
     char *argv[])
{
  try {
    string filename = argv[1];
    int hdunumber = atoi(argv[2]);

    FITS::Hdu h(filename, HDUAny, hdunumber, ReadOnly);
    cout << "Extension name " << h.getName() << endl;
    cout << " type " << h.getType() << endl;
    cout << " number " << h.getNumber() << endl;
    cout << " Header dump:" << endl;
    cout << *(h.header()) << endl;

    FITS::Hdu h2("test3.fits", HDUNull, "NewGuy2", ReadWrite+Create);
    h2.header()->replace("FLOATKEY", 23., "xxxxxx");
    string ss("No more");
    h2.header()->replace("STRING", ss);

  } catch (std::runtime_error &m) {
    quit(m,1);
  }
}

