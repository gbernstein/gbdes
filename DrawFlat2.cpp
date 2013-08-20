// Make a FITS image of a DECam flat field, restoring the chip-to-chip normalizations
// This version to condense a single MEF version of the flat.
#include <fstream>
#include <sstream>
#include <map>
#include <iomanip>

#include "StringStuff.h"

#include "FitsImage.h"
#include "PixelMapCollection.h"

#include "FitSubroutines.h"
#include "DECamInfo.h"
#include "Statistics.h"

using namespace img;
using namespace astrometry;
using namespace stringstuff;
using namespace decam;

string usage = 
  "DrawFlat2: map a MEF flat field onto the sky plane.\n"
  "     Flat fields are converted to magnitudes with mean of zero.\n"
  "usage: DrawFlat2 <out pix scale> <flat field in> <output FITS> [fill value] [med box]\n"
  "     <out pix scale> is number of DECam pixels (0.264\"/pix) per output pixel\n"
  "     <flat field in> is name of flat-field file\n"
  "     <output FITS> is name of output FITS file portraying array on sky.\n"
  "     [fill value] is assigned to pixels with no CCD (default=-1)\n"
  "     [med box] is 'radius' of box around each sample point to take median (default=5)";

int
main(int argc, char *argv[])
{
  if (argc<5 || argc>7) {
    cerr << usage << endl;
    exit(1);
  }
  double pixelScale = atof(argv[1]);
  string flatfile = argv[2];
  string outFits = argv[3];
  const float noData = argc>4 ? atof(argv[4]) : -1.;
  int medbox = argc>5 ? atoi(argv[5]) : 5;
  
  // Convert rough input scale in DECam pixels into degrees per output pix:
  pixelScale *= 0.264 / 3600.;

  // Bounds of device pixel coordinates that we will try to map:
  Bounds<int> bDevice(16, 2032, 1, 4096);

  // Read in the device information
  map<string,Device> devices = decamInfo();

  DECamMap wcs;

  // Make bounding rectangle for whole exposure
  Bounds<double> bExpo;
  for (map<string,Device>::iterator i = devices.begin();
       i != devices.end();
       ++i) 
      bExpo += i->second.b;

  // Make the Image
  Bounds<int> bImage( static_cast<int> (floor(bExpo.getXMin()/pixelScale)),
		      static_cast<int> (ceil(bExpo.getXMax()/pixelScale)),
		      static_cast<int> (floor(bExpo.getYMin()/pixelScale)),
		      static_cast<int> (ceil(bExpo.getYMax()/pixelScale)) );
  Image<> bigFlat(bImage, noData);

  double sum = 0.;
  long n = 0;

  // Open the FITS file to count extensions
  int nHDUs;
  {
    FITS::FitsFile ff(flatfile);
    nHDUs = ff.HDUCount();
  }

  for (int iHDU = 1; iHDU < nHDUs; iHDU++) {
    FitsImage<> fi(flatfile, FITS::ReadOnly, iHDU);
    string detpos;
    if (!fi.header()->getValue("DETPOS",detpos)) {
      cerr << "Did not find DETPOS at extension " << iHDU << endl;
      exit(1);
    }
    stripWhite(detpos);

    // Move along if this is an unwanted device
    if (devices.count(detpos)==0) continue;

    // Convert device boundaries into pixel range
    Bounds<double> b = devices[detpos].b;
    Bounds<int> bpix( static_cast<int> (floor(b.getXMin()/pixelScale)),
		      static_cast<int> (ceil(b.getXMax()/pixelScale)),
		      static_cast<int> (floor(b.getYMin()/pixelScale)),
		      static_cast<int> (ceil(b.getYMax()/pixelScale)) );

    // Open the flat-field image for this extension
    Image<> flat = fi.use();

    // Loop over pixels
    for (int iy = bpix.getYMin(); iy <= bpix.getYMax(); iy++)
      for (int ix = bpix.getXMin(); ix <= bpix.getXMax(); ix++) {
	double xField = ix * pixelScale;
	double yField = iy * pixelScale;
	double xPix, yPix;
	wcs.toPix(xField, yField, xPix, yPix);
	Bounds<int> medbounds( static_cast<int> (floor(xPix-medbox)),
			       static_cast<int> (ceil(xPix+medbox)),
			       static_cast<int> (floor(yPix-medbox)),
			       static_cast<int> (ceil(yPix+medbox)));
	// Overlap with measurable region
	medbounds = medbounds & bDevice;
	if (!medbounds) continue;	// Skip if outside device

	// Get median of pixels in the box
	vector<float> values;
	for (int jy = medbounds.getYMin(); jy <= medbounds.getYMax(); jy++)
	  for (int jx = medbounds.getXMin(); jx <= medbounds.getXMax(); jx++)
	    values.push_back(flat(jx,jy));

	// Apply normalization correction
	bigFlat(ix, iy) = -2.5*log10( stats::median(values));
	sum += bigFlat(ix,iy);
	n++;
      }

  } // end Device loop.
  // Set mean of good pixels to zero:
  bigFlat -= sum/n;
  bigFlat.shift(1,1);
  FitsImage<>::writeToFITS(outFits, bigFlat, 0);
}
