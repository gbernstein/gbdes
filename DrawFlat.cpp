// Make a FITS image of a DECam flat field, restoring the chip-to-chip normalizations
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
  "DrawFlat: map a flat field onto the sky plane, restoring chip normalizations\n"
  "     Flat fields are converted to magnitudes with mean of zero.\n"
  "usage: DrawFlat <out pix scale> <flat field base> <output FITS> [fill value] [med box]\n"
  "     <out pix scale> is number of DECam pixels (0.264\"/pix) per output pixel\n"
  "     <flat field base> is name of flat-field files, append _01.fits, etc.\n"
  "        The SCALMEAN keyword will be needed to restore CCD scaling.\n"
  "     <output FITS> is name of output FITS file portraying array on sky.\n"
  "     [fill value] is assigned to pixels with no CCD (default=-1)\n"
  "     [med box] is 'radius' of box around each sample point to take median (default=5)";

int
main(int argc, char *argv[])
{
  if (argc<4 || argc>6) {
    cerr << usage << endl;
    exit(1);
  }
  try {
  double pixelScale = atof(argv[1]);
  string flatbase = argv[2];
  string outFits = argv[3];
  const float noData = argc>4 ? atof(argv[4]) : -1.;
  int medbox = argc>5 ? atoi(argv[5]) : 5;
  
  // Convert rough input scale in DECam pixels into degrees per output pix:
  pixelScale *= 0.264 / 3600.;

  // Bounds of device pixel coordinates that we will try to map:
  Bounds<int> bDevice(16, 2033, 16, 4081);

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

  for (map<string,Device>::iterator i = devices.begin();
       i != devices.end();
       ++i) {
    wcs.setDevice(i->first);
    // Convert device boundaries into pixel range
    Bounds<double> b = i->second.b;
    Bounds<int> bpix( static_cast<int> (floor(b.getXMin()/pixelScale)),
		      static_cast<int> (ceil(b.getXMax()/pixelScale)),
		      static_cast<int> (floor(b.getYMin()/pixelScale)),
		      static_cast<int> (ceil(b.getYMax()/pixelScale)) );

    // Open the flat-field image for this extension
    ostringstream oss;
    oss << flatbase << "_" << setfill('0') << setw(2)  << i->second.ccdnum << ".fits" ;
    string flatfile = oss.str();
    Image<> flat;
    {
      FitsImage<> flatfits(flatfile,FITS::ReadOnly,0);
      flat = flatfits.extract();
    }
    double norm;
    if (!flat.getHdrValue("SCALMEAN",norm)) {
      cerr << "Flat field " << flatfile << " does not have SCALMEAN keyword" << endl;
      exit(1);
    }
    /**/cerr << "File " << flatfile << " SCALMEAN " << norm << endl;
    double normMag = -2.5*log10(norm);
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
	    if (flat(jx, jy)>0) 
	      values.push_back(flat(jx,jy));

	// Apply normalization correction
	double median = stats::median(values);
	bigFlat(ix, iy) = -2.5*log10(median) + normMag;
	if (median <=0 || median > 100) 
	  cerr << "median " << median << " for range " << medbounds << endl;
	sum += bigFlat(ix,iy);
	n++;
      }

  } // end Device loop.
  // Set mean of good pixels to zero:
  bigFlat -= sum/n;
  bigFlat.shift(1,1);
  FitsImage<>::writeToFITS(outFits, bigFlat, 0);
  } catch (std::runtime_error& e) {
    quit(e,1);
  }
}
