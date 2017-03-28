// test the FTable stuff.
#include "FTable.h"
#include <iostream>
#include <vector>
#include "Std.h"

using namespace std;

int
main(int argc,
     char *argv[])
{
  try {
    img::FTable ft;
    vector<double> count(20);
    for (int i=0; i<20; i++) count[i] = i;
    ft.addColumn(count, "COUNT");
    vector<int> odds(30);
    for (int i=0; i<30; i++) odds[i] = 2*i+1;
    ft.addColumn(odds, "ODDS");

    cout << "Table columns: " << ft.ncols() << endl;
    cout << "Table rows: " << ft.nrows() << endl;

    ft.eraseRows(7,12);
    ft.insertRows(7,5);
    ft.eraseRows(23);

    cout << "rows after erase: " << ft.nrows() << endl;

    vector<int> even(5);
    for (int i=0; i<even.size(); i++) even[i] = 2*i;
    ft.writeCells(even, "ODDS", 20);

    img::FTable ft2 = ft.duplicate();
    //img::FTable ft2 = ft;

    double xx;
    int ii;
    vector<string> colnames = ft2.columnNames();
    cout << "Row ";
    for (int i=0; i<colnames.size(); i++) cout << colnames[i] << "  ";
    cout << endl;
    for (int i=0; i<ft2.nrows(); i++) {
      ft2.readCell(xx, "COUNT", i);
      ft2.readCell(ii, "ODDS", i);
      cout << i << " " << xx << " " << ii << endl;
    }

    //    ft2.addColumn(vector<string>(), "STRINGS", 10);
    ft2.addColumn(vector<string>(), "STRINGS");
    /**/cerr << "constructed...";
    string s="A string";
    ft2.writeCell(s, "STRINGS", 7);
    /**/cerr << "wrote..." ;
    vector<string> vv;
    ft.readCells(vv, "STRINGS", 0, 9);
    /**/cerr << "read" << endl;
    for (int i=0; i<vv.size(); i++)
      cout << i << ": <" << vv[i] << ">" << endl;

    s = "Very long string";
    ft2.writeCell(s, "STRINGS", 15);
  } catch (std::runtime_error &m) {
    quit(m,1);
  }
}

