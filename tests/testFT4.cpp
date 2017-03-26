// test the FTable stuff.
#include <iostream>
#include <vector>
#include "Std.h"

#include "FITSTable.h"

using namespace std;
using namespace img;

int
main(int argc,
     char *argv[])
{
  try {
    string filename = argv[1];
    int hdunumber = atoi(argv[2]);

    /**/cerr << "Make table: " << endl;
    FTable ft;
    vector<float> f(3,2.);
    /**/cerr << "Add floats: " << endl;
    ft.addColumn(f, "Floats");
    vector<string> vs(3);
    vs[0] = "Hello";
    vs[1] = "there";
    vs[2] = "guy";
    /**/cerr << "Add strings: " << endl;
    ft.addColumn(vector<vector<string> >(), "Strings", 3, 5);
    ft.writeCell(vs, "Strings", 0);
    ft.writeCell(vs, "Strings", 1);
    ft.writeCell(vs, "Strings", 2);
    /**/cerr << "Add bools" << endl;
    vector<bool> vb(3,true);
    vb[2]=false;
    ft.addColumn(vector<vector<bool> >(), "bools", -1);
    ft.writeCell(vb, "bools", 0);
    vb.push_back(false);
    ft.writeCell(vb, "bools", 2);
    vb.push_back(true);
    ft.writeCell(vb, "bools", 1);

    FITSTable fitstab(filename, FITS::ReadWrite + FITS::Create
		      + FITS::OverwriteHDU, hdunumber);
    fitstab.replaceWith(ft);

  } catch (std::runtime_error &m) {
    quit(m,1);
  }
}

