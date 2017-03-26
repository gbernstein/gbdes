// Simple program that prints extension, object numbers for clipped measurements after fitting
#include "Std.h"
#include "FitsTable.h"
#include "StringStuff.h"

const string usage = 
  "ListClipped: output extension and object numbers of objects marked as\n"
  "             clipped in PhotoFit or WCSFit output files.\n"
  " usage: ListClipped <FITS file>\n"
  "        <FITS file> is name of FITS binary table file output by fitting program\n"
  "        stdout will be extension and object numbers, one pair per line.";

int
main(int argc,
     char *argv[]) {
  try {
    if (argc!=2) {
      cerr << usage << endl;
      exit(1);
    }

    string filename = argv[1];

    cout << "#" << stringstuff::taggedCommandLine(argc, argv) << endl;

    // Find which extensions are fitting outputs
    FITS::FitsFile ff(filename);
    for (int i=1; i<ff.HDUCount(); i++) {
      string hduName;
      {
	FITS::Hdu h(filename, FITS::HDUAny, i);
	hduName = h.getName();
      }
      if (stringstuff::nocaseEqual(hduName, "WCSOUT") ||
	  stringstuff::nocaseEqual(hduName, "PHOTOOUT")) {
	// Dump clipped id's from this table
	FITS::FitsTable ft(filename, FITS::ReadOnly, i);
	vector<string> desiredColumns;
	desiredColumns.push_back("Extension");
	desiredColumns.push_back("Object");
	img::FTable tab = ft.extract("CLIP", desiredColumns);
	vector<LONGLONG> extn;
	vector<LONGLONG> obj;
	tab.readCells(extn,"Extension");
	tab.readCells(obj,"Object");
	for (int j = 0; j<obj.size(); j++)
	  cout << extn[j] << " " << obj[j] << endl;
      }
    }
  } catch (std::runtime_error& m) {
    quit(m,1);
  }
}
