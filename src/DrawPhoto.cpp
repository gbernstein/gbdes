// Make a FITS image of a DECam photometric solution
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
  "DrawPhoto: make a FITS image a photometric solution, mapping all CCDs into sky plane.\n"
  "usage: DrawPhoto <out pix scale> <photo file> <solution name> <output FITS> \n"
  "                 [ref color] [fill value]\n"
  "     <out pix scale> is number of DECam pixels (0.264\"/pix) per output pixel\n"
  "     <photo file> is file holding the serialized solution from PhotoFit\n"
  "     <solution name> is the prefix of the name of the photometric solution to use from the file\n"
  "     <output FITS> is name of output FITS file portraying array on sky.\n"
  "     [ref color] is color used to draw the star flat (default 0.61)\n"
  "     [fill value] is placed in pixels not falling on the array.  (default=-1)";

int
main(int argc, char *argv[])
{
  if (argc<5 || argc>6) {
    cerr << usage << endl;
    exit(1);
  }
  double pixelScale = atof(argv[1]);
  string photoFile = argv[2];
  string photoPrefix = argv[3];
  photoPrefix += "/";
  string outFits = argv[4];
  const double refColor = argc>5 ? atof(argv[5]) : 0.61;
  const float noData = argc>6 ? atof(argv[6]) : -1.;
  
  loadPhotoMapParser();

  // Convert rough input scale in DECam pixels into degrees per output pix:
  pixelScale *= 0.264 / 3600.;

  // Bounds of device pixel coordinates that we will try to map:
  Bounds<double> bDevice(11., 2038., 1., 4096.);

  // Read in the device information
  map<string,Device> devices = decamInfo();

  DECamMap wcs;

  // Make bounding rectangle for whole exposure
  Bounds<double> bExpo;
  for (map<string,Device>::iterator i = devices.begin();
       i != devices.end();
       ++i) {
    bExpo += i->second.b;
  }
  // Make the Image
  Bounds<int> bImage( static_cast<int> (floor(bExpo.getXMin()/pixelScale)),
		      static_cast<int> (ceil(bExpo.getXMax()/pixelScale)),
		      static_cast<int> (floor(bExpo.getYMin()/pixelScale)),
		      static_cast<int> (ceil(bExpo.getYMax()/pixelScale)) );
  Image<> starflat(bImage, noData);
  Image<> colorterm(bImage, noData);

  PhotoMapCollection photomaps;
  {
    string filename = photoFile;
    ifstream ifs(filename.c_str());
    photomaps.read(ifs);
  }

  for (auto& i : devices) {
    // Convert device boundaries into pixel range
    Bounds<double> b = i.second.b;
    Bounds<int> bpix( static_cast<int> (floor(b.getXMin()/pixelScale)),
		      static_cast<int> (ceil(b.getXMax()/pixelScale)),
		      static_cast<int> (floor(b.getYMin()/pixelScale)),
		      static_cast<int> (ceil(b.getYMax()/pixelScale)) );

    // Get Wcs and photometric solution
    wcs.setDevice(i.first);
    photometry::SubMap* devPhoto = 0;
    try {
      devPhoto = photomaps.issueMap(photoPrefix + i.first);
    } catch (photometry::PhotometryError& pm) {
      // Might not have data for some extensions.  Just skip them.
      /**/cerr << "Skipping " << photoPrefix + i.first << endl;
      devPhoto = 0;
    }
    if (devPhoto) {
      // Loop over pixels
      for (int iy = bpix.getYMin(); iy <= bpix.getYMax(); iy++)
	for (int ix = bpix.getXMin(); ix <= bpix.getXMax(); ix++) {
	  photometry::PhotoArguments args;
	  args.xExposure = ix * pixelScale;
	  args.yExposure = iy * pixelScale;
	  wcs.toPix(args.xExposure, args.yExposure, args.xDevice, args.yDevice);
	  if (!bDevice.includes(args.xDevice,args.yDevice)) continue;	// Skip if outside device

	  // Calculate photometric correction
	  args.color = refColor;
	  double photo = devPhoto->forward(0., args);
	  args.color = refColor+1.;
	  colorterm(ix,iy) = devPhoto->forward(0.,args) - photo;
	  starflat(ix, iy) = -photo;
	}
    }
  } // end Device loop.
  starflat.shift(1,1);
  FitsImage<>::writeToFITS(outFits, starflat, 0);
  colorterm.shift(1,1);
  FitsImage<>::writeToFITS(outFits, colorterm, 1);
}
