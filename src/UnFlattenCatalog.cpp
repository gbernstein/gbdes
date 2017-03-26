#include <fstream>
#include <sstream>
#include <map>

#include "StringStuff.h"
#include "FitsTable.h"
#include "FitsImage.h"
#include "ReadLdacHeader.h"
#include "DECamInfo.h"

using namespace img;
using namespace FITS;
using namespace stringstuff;

string usage = "UnFlattenCatalog: Collect individual-CCD LDAC catalogs into a single FITS file,\n"
  "     along the way multiplying all magnitudes by a factor derived from lookup into the\n"
  "     superpixel representation of the sky flat used by DESDM reductions.\n"
  "     Also filtering the catalog down to columns of interest and magnitude range of interest.\n"
  "\n"
  " usage: UnFlattenCatalog <superpix file> <output FITS> [max magerr] [superpix size]\n"
  "     <superpix file> is name of the star flat superpixel file containing sky flat\n"
  "                  corrections that we want to un-correct on each magnitude.\n"
  "                  Enter 'NULL' if no magnitude corrections wanted.\n"
  "     <output FITS> is name of output combined FITS catalog file\n"
  "     [max magerr] is the largest MAGERR_AUTO allowed into the output catalogs (default=0.05)\n"
  "     [superpix size] is number of CCD pixels per pixel on tne input image (default 512)\n"
  "     stdin is a list of individual LDAC FITS files to be combined into this output FITS catalog.\n"
  "  The file detsecs.dat holding device corners must be in the current directory.";

int
main(int argc, char *argv[])
{
  if (argc>5 || argc<3) {
    cerr << usage << endl;
    exit(1);
  }
  string superpixFile = argv[1];
  bool useStarFlat = !stringstuff::nocaseEqual(superpixFile, "NULL");
  string outFits = argv[2];
  const float maxMagErr = argc > 3 ? atof(argv[3]) : 0.05;
  const int superPixelSize = argc > 4 ? atoi(argv[4]) : 512;

  try {
    // Read in the device information
    map<string,decam::Device> devices = decam::decamInfo();

    // Read the starflat image
    Image<> starflat;
    if (useStarFlat) {
      FitsImage<> f(superpixFile,FITS::ReadOnly,0);
      starflat = f.extract();
    }

    // Names of columns that we want to keep from Final Cut catalogs:
    vector<string> keepColumns;
    keepColumns.push_back("BACKGROUND");
    keepColumns.push_back("FLUX_MAX");
    keepColumns.push_back("FLAGS");
    keepColumns.push_back("SPREAD_MODEL");
    keepColumns.push_back("FLUX_RADIUS");
    keepColumns.push_back("MAG_APER");
    keepColumns.push_back("MAGERR_APER");
    keepColumns.push_back("MAG_AUTO");
    keepColumns.push_back("MAGERR_AUTO");
    keepColumns.push_back("MAG_PSF");
    keepColumns.push_back("MAGERR_PSF");
    keepColumns.push_back("XWIN_IMAGE");
    keepColumns.push_back("YWIN_IMAGE");
    keepColumns.push_back("ERRAWIN_IMAGE");
    keepColumns.push_back("ERRBWIN_IMAGE");
    keepColumns.push_back("ERRTHETAWIN_IMAGE");

    // Our selection criterion:
    string selection;
    {
      ostringstream oss;
      oss << "MAGERR_AUTO < " << maxMagErr;
      selection = oss.str();
    }

    // Range of valid positions on the array
    Bounds<int> array(1,2048,1,4096);

    // Open the output catalog file

    string inFits;
    bool firstWriteToFile = true;
    while (getlineNoComment(cin, inFits)) {
      // Loop over input catalogs to go into output
      stripWhite(inFits);
      // Copy over an LDAC header extension
      {
	FitsTable hdrTable(inFits, FITS::ReadOnly, 1);
	if (!stringstuff::nocaseEqual(hdrTable.getName(), "LDAC_IMHEAD")) {
	  cerr << "Extension 1 is not LDAC_HEADER for file " << inFits << endl;
	  exit(1);
	}
	if (firstWriteToFile) {
	  FitsTable hdrOut(outFits, FITS::Create + FITS::OverwriteFile, -1) ;
	  hdrOut.adopt( hdrTable.extract());
	  firstWriteToFile = false;
	  hdrOut.setName("LDAC_IMHEAD");
	} else {
	  FitsTable hdrOut(outFits, FITS::Create, -1) ;
	  hdrOut.adopt( hdrTable.extract());
	  hdrOut.setName("LDAC_IMHEAD");
	}
      }

      // Extract the CCDNAME
      string detPos;
      img::Header ldacHead = ReadLdacHeader(inFits, 1);
      if (!ldacHead.getValue("DETPOS",detPos)) {
	cerr << "Did not find DETPOS in header of " << inFits << endl;
	exit(1);
      }
      stripWhite(detPos);

      // Get super-array origin coordinates
      int detsecX = devices[detPos].detsecX;
      int detsecY = devices[detPos].detsecY;
    
      // Output table where this will go:
      FitsTable ff2(outFits, FITS::Create, -1);
      ff2.setName("LDAC_OBJECTS");

      FTable ft;
      {
	// Extract catalog
	FitsTable ff(inFits, FITS::ReadOnly, 2);
	ft = ff.extract(selection, keepColumns);
      }

      if (useStarFlat) {
	// Adjust magnitudes of each object
	const long nrows = ft.nrows();
	vector<double> x(nrows);
	vector<double> y(nrows);
	vector<float> magauto(nrows);
	vector<float> magpsf(nrows);
	vector<vector<float> > magaper(nrows);
	ft.readCells(x, "XWIN_IMAGE");
	ft.readCells(y, "YWIN_IMAGE");
	ft.readCells(magauto, "MAG_AUTO");
	ft.readCells(magpsf, "MAG_PSF");
	ft.readCells(magaper, "MAG_APER");
	for (long i = 0; i < nrows; i++) {
	  int xPix = static_cast<int> (floor(x[i]));
	  int yPix = static_cast<int> (floor(y[i]));
	  if ( array.includes(xPix, yPix)) {
	    // Find the pixel in Yanny's image that gives superpixel here
	    int xSup = (xPix + detsecX-2)/superPixelSize + 1;
	    int ySup = (yPix + detsecY-2)/superPixelSize + 1;
	    double magAdjust = -starflat(xSup,ySup);
	    magauto[i] += magAdjust;
	    magpsf[i] += magAdjust;
	    for (int j = 0; j<magaper[i].size(); j++)
	      magaper[i][j] += magAdjust;
	  }
	}
	ft.writeCells(magauto,"MAG_AUTO");
	ft.writeCells(magpsf,"MAG_PSF");
	ft.writeCells(magaper,"MAG_APER");
      }
      // New output table adopts the filtered & adjusted data
      ff2.adopt(ft);
    } // End input catalog loop
  } catch (std::runtime_error& m) {
    quit(m,1);
  }
}
