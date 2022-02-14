// Test Lookup1d functions
#include <iostream>
#include <cstdlib>

#include "Lookup1d.h"
using namespace lookup;
using namespace std;

string usage = 
  "Basic trials of the Lookup1d class.  \n"
  "usage: test_lookup <filename> <mapname> \n"
  "  <filename> is a YAML file holding some templates.  In current directory, or\n"
  "             can specify a search path with CAL_PATH environment variable.\n"
  "  <mapname>  is name of one of the maps in the file that you want to probe.\n\n"
  "If file or map do not exist, should get appropriate message then quit.\n"
  "Should throw an exception trying to find name in the 'wrong' cache.  Then proceed\n"
  "to ask for x,y pairs, which it will use to look up in the table and give value and slope.\n"
  "Also tries the inverse function which should yield equality with the input arg value.";


int main(int argc,
	 char *argv[]) {
  if (argc != 3) {
    cerr << usage << endl;
    exit(1);}
  string filename = argv[1];
  string mapname = argv[2];

  const Lookup1d* lu=0;
  try {
    Lookup1d::ingestFile(filename);
    lu = Lookup1d::find(mapname);
    Lookup1d::find(mapname, "wrong"); // look in non-existent cache
  } catch (LookupError& m) {
    cerr << m.what() << endl;
  }

  if (!lu) exit(1);

  cerr << "Enter x,y: ";
  double x,y;
  while (cin >> x >> y) {
    double v,s;
    lu->valueSlope(x,y,v,s);
    cerr << v << " " << s << " at arg " << lu->getArg(x,y) <<  endl;
    v += lu->getArg(x,y);
    cerr << "invert: " << lu->inverse(v,1.) << endl;
    cerr << "Enter x,y: ";
  }
  exit(0);
}
