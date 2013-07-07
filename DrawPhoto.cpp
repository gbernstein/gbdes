// Make a FITS image of a DECam photometric solution, with flat normalization factors taken out.
#include <fstream>
#include <sstream>
#include <map>

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
  map<string,Device> devices;
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
      d.norm = norm;
      d.b = Bounds<double>(xmin, xmax, ymin, ymax);
      bExpo += d.b;
      devices[name] = d;
    }
  }

  // Make the Image
  Bounds<int> bImage( static_cast<int> (floor(bExpo.getXMin()/pixelScale)),
		      static_cast<int> (ceil(bExpo.getXMax()/pixelScale)),
		      static_cast<int> (floor(bExpo.getYMin()/pixelScale)),
		      static_cast<int> (ceil(bExpo.getYMax()/pixelScale)) );
  Image<> starflat(bImage);

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

  // Adjust coordinates such that mean of N4 and S4 centers is at origin.
  Position<double> center;
  {
    Bounds<double> b;
    double x,y;
    astromaps.issueWcs(exposurePrefix + "N4")->toWorld(1024.,2048.,x,y);
    b += Position<double>(x,y);
    astromaps.issueWcs(exposurePrefix + "S4")->toWorld(1024.,2048.,x,y);
    b += Position<double>(x,y);
    center = b.center();
    cout << "Center: " << center*3600. << " arcsec" << endl;
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
    Wcs* devWcs = astromaps.issueWcs(exposurePrefix + i->first);
    photometry::SubMap* devPhoto = photomaps.issueMap(exposurePrefix + i->first);

    // Loop over pixels
    for (int iy = bpix.getYMin(); iy <= bpix.getYMax(); iy++)
      for (int ix = bpix.getXMin(); ix <= bpix.getXMax(); ix++) {
	photometry::PhotoArguments args;
	args.xExposure = ix * pixelScale + center.x;
	args.yExposure = iy * pixelScale + center.y;

	// Map pixel coordinates to exposure and device coordinates
	devWcs->toPix(args.xExposure, args.yExposure, args.xDevice, args.yDevice);
	if (!bDevice.includes(args.xDevice,args.yDevice)) continue;	// Skip if outside device

	// Calculate photometric correction
	args.color = 0.;
	double photo = devPhoto->forward(0., args);

	// Apply normalization correction
	starflat(ix, iy) = pow(10., -0.4*photo) * i->second.norm;
      }

  } // end Device loop.
  starflat.shift(1,1);
  FitsImage<>::writeToFITS("test.fits", starflat, 0);
}
