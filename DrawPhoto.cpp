// Make a FITS image of a DECam photometric solution, with flat normalization factors taken out.
#include <fstream>
#include <sstream>
#include "StringStuff.h"

#include "FitsImage.h"
#include "PhotoMapCollection.h"
#include "PixelMapCollection.h"
using namespace img;
using namespace astrometry;
using namespace stringstuff;
using photometry::PhotoMapCollection;
using photometry::PhotoMap;

class Device {
public:
  string name;
  Bounds<double> b;
  double norm;
};

int
main(int argc, char *argv[])
{
  double pixelScale = atof(argv[1]);
  // Convert rough input scale in DECam pixels into degrees per output pix:
  pixelScale *= 0.264 / 3600.;

  // Bounds of device pixel coordinates that we will try to map:
  Bounds<double> bDevice(11., 2038., 1., 4096.);

  // Read in the device information
  list<Device> devices;
  Bounds<double> bExpo;
  {
    ifstream ifs("ccdInfo_i.txt");
    string buffer;
    while (getlineNoComment(ifs, buffer)) {
      istringstream iss(buffer);
      string name;
      double xmin, xmax, ymin, ymax, norm;
      if (!(iss >> name >> xmin >> xmax >> ymin >> ymax >> norm)) {
	cerr << "Bad Device line: " << buffer << endl;
	exit(1);
      }
      Device d;
      d.name = name;
      d.norm = norm;
      d.b = Bounds<double>(xmin, xmax, ymin, ymax);
      bExpo += d.b;
      devices.push_back(d);
    }
  }

  // Make the Image
  Bounds<int> bImage( static_cast<int> (floor(bExpo.getXMin()/pixelScale)),
		      static_cast<int> (ceil(bExpo.getXMax()/pixelScale)),
		      static_cast<int> (floor(bExpo.getYMin()/pixelScale)),
		      static_cast<int> (ceil(bExpo.getYMax()/pixelScale)) );
  Image<> starflat(bImage);

  /**/cerr << "degrees: " << bExpo 
	   << " Image: " << bImage
	   << endl;

  // Read in astrometric and photometric solutions
  PixelMapCollection astromaps;
  {
    ifstream ifs("jul4.wcs");
    astromaps.read(ifs);
  }

  PhotoMapCollection photomaps;
  {
    ifstream ifs("jul4.i.photo");
    photomaps.read(ifs);
  }

  //  string exposurePrefix = "DECam_00163584/";
  string exposurePrefix = "DECam_00163609/";

  for (list<Device>::iterator i = devices.begin();
       i != devices.end();
       ++i) {
    // Convert device boundaries into pixel range
    Bounds<int> bpix( static_cast<int> (floor(i->b.getXMin()/pixelScale)),
		      static_cast<int> (ceil(i->b.getXMax()/pixelScale)),
		      static_cast<int> (floor(i->b.getYMin()/pixelScale)),
		      static_cast<int> (ceil(i->b.getYMax()/pixelScale)) );

    // Get Wcs and photometric solution
    Wcs* devWcs = astromaps.issueWcs(exposurePrefix + i->name);
    photometry::SubMap* devPhoto = photomaps.issueMap(exposurePrefix + i->name);

    // Loop over pixels
    for (int iy = bpix.getYMin(); iy <= bpix.getYMax(); iy++)
      for (int ix = bpix.getXMin(); ix <= bpix.getXMax(); ix++) {
	photometry::PhotoArguments args;
	args.xExposure = ix * pixelScale;
	args.yExposure = iy * pixelScale;

	// Map pixel coordinates to exposure and device coordinates
	devWcs->toPix(args.xExposure, args.yExposure, args.xDevice, args.yDevice);
	if (!bDevice.includes(args.xDevice,args.yDevice)) continue;	// Skip if outside device

	// Calculate photometric correction
	args.color = 0.;
	double photo = devPhoto->forward(0., args);

	// Apply normalization correction
	starflat(ix, iy) = pow(10., -0.4*photo) * i->norm; // ??? divide or multiply?
      }

  } // end Device loop.
  starflat.shift(1,1);
  FitsImage<>::writeToFITS("test.fits", starflat, 0);
}

