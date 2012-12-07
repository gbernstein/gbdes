// Glue together a list of FITS files into a single one with more extensions.
// Takes the 0 extension of just the first one unless you give a 2nd argument.
// Usage: FitsGlue outfile <use primaries>
// stdin is list of input files, one per line.

#include "FITS.h"
#include <iostream>
#include "StringStuff.h"

using namespace std;
using namespace FITS;

int
main(int argc,
     char *argv[])
{
  string outfile = argv[1];
  bool usePrimaries = (argc > 2);

  int status = 0;
  fitsfile* fin;
  fitsfile* fout;
  fits_create_file(&fout, const_cast<char*> (outfile.c_str()), &status);

  string buffer;
  int iextn = 0;
  while (stringstuff::getlineNoComment(cin, buffer)) {
    fits_open_file(&fin, const_cast<char*> (buffer.c_str()), READONLY, &status);
    if (!(iextn == 0 || usePrimaries))
      fits_movabs_hdu(fin, 2, 0, &status);
    fits_copy_file(fin, fout, 0, 1, 1, &status);
    fits_close_file(fin, &status);
    /**/cerr << "Status " << status << " Gluing " << buffer << endl;
    iextn++;
  }
  fits_close_file(fout, &status);
}
