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

    vector<string> cn;


    FitsTable h(filename, ReadOnly, hdunumber);
    cout << "Extension name " << h.getName() << endl;
    cout << " type " << h.getType() << endl;
    cout << " number " << h.getNumber() << endl;
    cout << " Header dump:" << endl;
    cout << *(h.header()) << endl;

    img::FTable te = h.extract();
    /**/cerr << "*** done first extract ***" << endl;
    /**/cerr << "MAG_APER repeat: " << te.repeat("MAG_APER") << endl;
    set<long> finders;
    finders.insert(16);
    finders.insert(17);
    finders.insert(19);
    img::FTable t2 = te.extractRows("MAG_AUTO < 16 && FLAGS");
    /**/cerr << "*** done second extract ***" << endl;
    /**/cerr << "MAG_APER repeat: " << t2.repeat("MAG_APER") << endl;

    {
      FitsTable ftnew("junk.fits", ReadWrite + OverwriteFile, "Extract");
      ftnew.copy(t2);
      cerr << "*** done copy ***" << endl;
      cerr << "MAG_APER repeat: " << t2.repeat("Mag_APER") << endl;
    }

    cerr << "destroyed Extract" << endl;

    FitsTable ftnew("junk.fits", ReadWrite+Create, "Filter");
    cerr << "copy again:***" << endl;
    ftnew.copy(t2);

    img::FTable t3=ftnew.use();
    t3.filterRows("FLAGS", finders);
    cerr << "t3 filtered row count " << t3.nrows() << endl;
    /**/
    cerr << "Done" << endl;
  } catch (std::runtime_error &m) {
    quit(m,1);
  }
}

