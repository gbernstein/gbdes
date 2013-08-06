// Make a FITS image of a DECam photometric solution, with flat normalization factors taken out.
#include <fstream>
#include <sstream>
#include <map>

#include "StringStuff.h"

#include "FitsImage.h"
#include "PhotoMapCollection.h"
#include "PixelMapCollection.h"

#include "FitSubroutines.h"
#include "DECamInfo.h"

using namespace img;
using namespace astrometry;
using namespace stringstuff;
using namespace decam;

using photometry::PhotoMapCollection;
using photometry::PhotoMap;

string usage = 
  "DrawPhoto: map a photometric solution onto the sky plane, with possible correction for\n"
  "     the normalizations of the underlying dome flats.\n"
  "usage: DrawPhoto <out pix scale> <photo file> <solution name> <output FITS> "
  " [norm file] [fill value]\n"
  "     <out pix scale> is number of DECam pixels (0.264\"/pix) per output pixel\n"
  "     <photo file> is file holding the serialized solution from PhotoFit\n"
  "     <solution name> is the prefix of the name of the photometric solution to use from the file\n"
  "     <output FITS> is name of output FITS file portraying array on sky.\n"
  "     [norm file] is optional file containing dome flat norm factors to adjust for.\n"
  "           This file should have <device name> <norm factor> pair on each line.\n"
  "     [fill value] is placed in pixels not falling on the array.  (default=-1)";

int
main(int argc, char *argv[])
{
  if (argc<5 || argc>7) {
    cerr << usage << endl;
    exit(1);
  }
  double pixelScale = atof(argv[1]);
  string photoFile = argv[2];
  string photoPrefix = argv[3];
  photoPrefix += "/";
  string outFits = argv[4];
  string normFile = argc>5 ? argv[5] : "";
  const float noData = argc>6 ? atof(argv[6]) : -1.;
  
  loadPhotoMapParser();

  // Convert rough input scale in DECam pixels into degrees per output pix:
  pixelScale *= 0.264 / 3600.;

  // Bounds of device pixel coordinates that we will try to map:
  Bounds<double> bDevice(11., 2038., 1., 4096.);

  // Read in the device information
  map<string,Device> devices = decamInfo();
  // And norms if wanted
  if (!normFile.empty())
    getDeviceNorms(normFile, devices);

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
  Image<> starflat(bImage, noData);
  Image<> colorterm(bImage, noData);

  double sum = 0.;
  long n = 0;

  PhotoMapCollection photomaps;
  {
    string filename = photoFile;
    ifstream ifs(filename.c_str());
    photomaps.read(ifs);
  }

  for (map<string,Device>::iterator i = devices.begin();
       i != devices.end();
       ++i) {
    // Convert device boundaries into pixel range
    Bounds<double> b = i->second.b;
    Bounds<int> bpix( static_cast<int> (floor(b.getXMin()/pixelScale)),
		      static_cast<int> (ceil(b.getXMax()/pixelScale)),
		      static_cast<int> (floor(b.getYMin()/pixelScale)),
		      static_cast<int> (ceil(b.getYMax()/pixelScale)) );

    // Get Wcs and photometric solution
    wcs.setDevice(i->first);
    photometry::SubMap* devPhoto = photomaps.issueMap(photoPrefix + i->first);
    double normMag = -2.5*log10(i->second.norm);
    // Loop over pixels
    for (int iy = bpix.getYMin(); iy <= bpix.getYMax(); iy++)
      for (int ix = bpix.getXMin(); ix <= bpix.getXMax(); ix++) {
	photometry::PhotoArguments args;
	args.xExposure = ix * pixelScale;
	args.yExposure = iy * pixelScale;
	wcs.toPix(args.xExposure, args.yExposure, args.xDevice, args.yDevice);
	if (!bDevice.includes(args.xDevice,args.yDevice)) continue;	// Skip if outside device

	// Calculate photometric correction
	args.color = 0.;
	double photo = devPhoto->forward(0., args);
	args.color = 1;
	colorterm(ix,iy) = devPhoto->forward(0.,args) - photo;
	// Apply normalization correction
	starflat(ix, iy) = -(photo + normMag);
	sum += starflat(ix,iy);
	n++;
      }

  } // end Device loop.
  // Set mean of good pixels to zero:
  starflat -= sum/n;
  starflat.shift(1,1);
  FitsImage<>::writeToFITS(outFits, starflat, 0);
  colorterm.shift(1,1);
  FitsImage<>::writeToFITS(outFits, colorterm, 1);
}
