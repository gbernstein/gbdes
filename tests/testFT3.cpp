// test the FTable stuff.
#include <iostream>
#include <vector>
#include "Std.h"

#include "Hdu.h"

using namespace std;
using namespace img;

int
main(int argc,
     char *argv[])
{
  string filename = argv[1];
  int hdunumber = atoi(argv[2]);
  double fff=atof(argv[3]);
  try {
    /**/cerr << "Opening hdu...";
    FITS::Hdu ext(filename,
		  FITS::HDUAny,
		  hdunumber,
		  FITS::ReadWrite + FITS::Create);
    
    /**/cerr << "setName...";
    ext.setName("Testing");
    /**/cerr << "replace...";
    ext.header()->replace("NEWKEY", fff, "A comment");
  } catch (std::runtime_error &m) {
      cerr << "First block catches: " << m.what() << endl;
  }
  /**/cerr << "done1..." << endl;

  try {
    FITS::Hdu ext(filename,
		  FITS::HDUAny,
		  "Testing");
    cerr << "Re-opened, at HDU number " << ext.getNumber() << endl; 
    cerr << "isWriteable()? " << ext.isWriteable() << endl;
    cerr << "Header isLocked()? " << ext.header()->isLocked() << endl;
   try {
      ext.header()->getValue("NEWKEY", fff);
      cerr << "recovered NEWKEY " << fff << endl;
    } catch (std::runtime_error &m) {
      cerr << "hdr block1 catches: " << m.what() << endl;
    }
  } catch (std::runtime_error &m) {
      cerr << "Second block catches: " << m.what() << endl;
  }
}

