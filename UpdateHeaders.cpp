// Program to modify headers of multi-extension FITS files 
#include "Hdu.h"
#include <fstream>
#include <sstream>
#include "StringStuff.h"
#include "HeaderFromStream.h"

using namespace FITS;
using namespace stringstuff;

string usage = "UpdateHeaders: replace or add a set of keywords to all extensions in a DECam image\n"
  "Usage: UpdateHeaders <headerFile> <fitsFile> [fitsFile...]\n"
  "headerFile is a text file containing lines of format\n"
  "    <ccdnum>  <filter> <valid FITS header card...>\n"
  "fitsFile is the name of FITS file(s) to be altered\n"
  "  An input card is applied to an extension if its CCDNUM header card value matches the one given\n"
  "on the input line.  Use '.' to apply a change to all FITS extensions.\n"
  "  Also the (start of) FILTER keyword value must match for card to be applied.  Value of '.'\n"
  "will match any filter.  The FILTER keyword in the primary header will be used if there is none\n"
  "in a particular extension.";

int
main(int argc, char *argv[]) {

  if (argc <3) {
    cerr << usage << endl;
    exit(1);
  }

  try {
    string infile = argv[1];

    ifstream ifs(infile.c_str());
    vector<int> ccdnums;
    vector<string> filters;
    vector<string> cards;

    string buffer;
    while (getlineNoComment(ifs, buffer)) {
      istringstream iss(buffer);
      string ccdnumString;
      string filter;
      iss >> ccdnumString >> filter;
      if (ccdnumString==".") {
	ccdnums.push_back(-1);
      } else {
	ccdnums.push_back(atoi(ccdnumString.c_str()));
      }
      filters.push_back(filter);
      string card;
      getline(iss, card);
      cards.push_back(card);
    }

    for (int iarg=2; iarg<argc; iarg++) {
      string filename = argv[iarg];

      // get HDU count
      int hduCount = 0;
      {
	FitsFile ff(filename);
	hduCount = ff.HDUCount();
      }

      // Read the filter from primary header if there is such an entry
      string primaryFilter = "";
      const string filterKW = "FILTER";
      const string ccdnumKW = "CCDNUM";

      // Loop through all extensions, making modifications
      for (int ihdu=0; ihdu<hduCount; ihdu++) {
	Hdu h(filename, FITS::HDUAny, ihdu, FITS::ReadWrite);
	// Get primaryFilter if this is primary extension
	if (ihdu==0)
	  if (!h.header()->getValue(filterKW, primaryFilter))
	    primaryFilter="";

	// Now get CCDNUM from this extension, if any
	int ccdnum;
	if (!h.header()->getValue(ccdnumKW, ccdnum))
	  ccdnum = -1;

	// And get any local filter to supercede global
	string filter;
	if (!h.header()->getValue(filterKW, filter))
	  filter = primaryFilter;

	// Assemble the header cards applicable to this HDU
	string newHeader;
	{
	  ostringstream oss;
	  for (int i=0; i<filters.size(); i++) {
	    if ( (filters[i] == "." || nocaseEqual(filters[i], filter.substr(0,filters[i].size())))
		 && (ccdnums[i]<0 || ccdnum==ccdnums[i]) )
	      oss << cards[i] << endl;
	  }
	  newHeader = oss.str();
	}

	// Build header from the new ones and add it to the current header
	istringstream iss(newHeader);
	(*h.header())+=img::HeaderFromStream(iss);
      } // end HDU loop.
    } // end file loop
  } catch (std::runtime_error& m) {
    quit(m,1);
  }
}

