// Make a FITS image of a DECam photometric solution, with flat normalization factors taken out.
#include <fstream>
#include <sstream>
#include <map>
#include <ctime>
#include <iomanip>

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
  "Photo2DESDM: convert photometric instrumental solution to FITS files suitable\n"
  "     for use by DESDM\n"
  "usage: Photo2DESDM <photo file> <instrument name> <output FITS base> \n"
  "     <photo file> is file holding the serialized solution from PhotoFit\n"
  "     <instrument name> is the filter name, should be instrument name in photo file too\n"
  "     <output FITS base> is root of the files that will be drawn, one per CCD";

int
main(int argc, char *argv[])
{
  if (argc!=4) {
    cerr << usage << endl;
    exit(1);
  }
  string photoFile = argv[1];
  string filter = argv[2];
  string photoPrefix = filter + "/";
  string outFits = argv[3];
  
  try {
  loadPhotoMapParser();

  // Size of full devices:
  Bounds<int> bpix(1, 2048, 1, 4096);
  // Bounds of device pixel coordinates that we will use photo solutions.
  // They will just be unity outside these bounds.
  Bounds<int> bDevice(11, 2038, 1, 4096);

  // Read in the device information
  map<string,Device> devices = decamInfo();

  PhotoMapCollection photomaps;
  {
    string filename = photoFile;
    ifstream ifs(filename.c_str());
    photomaps.read(ifs);
  }

  for (map<string,Device>::iterator i = devices.begin();
       i != devices.end();
       ++i) {
    cerr << "Working on " << i->first << endl;
    photometry::SubMap* devPhoto = photomaps.issueMap(photoPrefix + i->first);

    int ccdNum = i->second.ccdnum;
    Image<> starflat(bpix, 1.);

    starflat.header()->append("CCDNUM", ccdNum);
    starflat.header()->append("FILTER", filter);
    time_t now;
    time(&now);
    string date = ctime(&now);
    starflat.header()->append("DESMKICR", date);
    string s = "illumcor";
    starflat.header()->append("OBSTYPE", s);
    s = "IMAGE ";
    starflat.header()->append("DES_EXT", s);
    starflat.header()->addHistory("Created from starflat solution " + photoFile);

    // Loop over pixels
    photometry::PhotoArguments args;
    args.xExposure = args.yExposure = args.color = 0.;

    for (int iy = bDevice.getYMin(); iy <= bDevice.getYMax(); iy++)
      for (int ix = bDevice.getXMin(); ix <= bDevice.getXMax(); ix++) {
	args.xDevice = ix;
	args.yDevice = iy;
	double photo = devPhoto->forward(0., args);
	starflat(ix, iy) = pow(10., -0.4*photo);
      }

    ostringstream oss;
    oss << outFits << "_" << std::setfill('0') << setw(2) << ccdNum << ".fits";
    FitsImage<>::writeToFITS(oss.str(), starflat, 0);
  } // end Device loop.
  } catch (std::runtime_error& e) {
    quit(e,1);
  }
}
