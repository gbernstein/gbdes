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

    cerr << "call use()" << endl;
    img::FTable ft = h.use();

    cerr << "call colnames" << endl;
    cn = ft.listColumns();
    for (int i=0; i<cn.size(); i++)
      cout << cn[i] << endl;

    FitsTable h2("junk.fits", ReadWrite + OverwriteFile , "glob");

    cerr << "...copy" << endl;
    h2.copy(ft);
    cerr << "...use" << endl;
    img::FTable ft2=h2.use();
    cerr << "...clear" << endl;
    ft2.clear();
    vector<double> vv(5,-24.7);
    cerr << "...addColumn" << endl;
    ft2.addColumn(vv, "UND_WOOD", 5);
    cn = ft2.listColumns();
    cout << "ft2 row count " << ft2.nrows() << ", columns: " << endl;
    for (int i=0; i<cn.size(); i++)
      cout << cn[i] << endl;
    
    cn = ft2.listColumns();
    cout << "...should match h2 columns: " << endl;
    for (int i=0; i<cn.size(); i++)
      cout << cn[i] << endl;

    FitsTable h3("junk.fits", ReadWrite + Create, "Adopted");
    img::FTable ft3 = h2.extract(1,2,vector<string>(1,"u.*"));
    cn = ft3.listColumns();
    cout << "ft3 row count " << ft3.nrows() << ", columns: " << endl;
    for (int i=0; i<cn.size(); i++)
      cout << cn[i] << endl;

    {
      vector<LONGLONG> count(20);
      for (int i=0; i<20; i++) count[i] = i;
      ft3.addColumn(count, "COUNT");
      vector<int> odds(30);
      for (int i=0; i<30; i++) odds[i] = 2*i+1;
      ft3.addColumn(odds, "ODDS");

      ft3.eraseRows(7,12);
      ft3.insertRows(7,5);
      ft3.eraseRows(23);
      
      vector<int> even(5);
      for (int i=0; i<even.size(); i++) even[i] = 2*i;
      ft3.writeCells(even, "ODDS", 20);
    }

    cout << "ft3 row count: " << ft3.nrows() << endl;
    cn = ft3.listColumns();
    for (int i=0; i<cn.size(); i++)
      cout << cn[i] << endl;

    cerr << "...adopting" << endl;
    h3.adopt(ft3);
    cerr << "...destroying" << endl;
    /**/
  } catch (std::runtime_error &m) {
    quit(m,1);
  }
}

