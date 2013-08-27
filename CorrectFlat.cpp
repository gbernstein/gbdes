// Apply photometric correction model to an existing flat-field image, taking log
#include <fstream>
#include <sstream>
#include <map>
#include <iomanip>

#include "StringStuff.h"

#include "FitsImage.h"
#include "PhotoMapCollection.h"

#include "FitSubroutines.h"
#include "DECamInfo.h"
#include "Statistics.h"

using namespace img;
using namespace photometry;
using namespace stringstuff;
using namespace decam;

string usage = 
  "CorrectFlat: Apply derived star flat correction to a dome flat image\n"
  "usage: CorrectFlat <flat in> <photo file> <instrument name> <output FITS>\n"
  "     <flat in> is name of input flat-field FITS file\n"
  "     <photo file> holds serialized star flat solutions from PhotoFit\n"
  "     <instrument name> identifies the instrument/device name for photo map to use\n"
  "     <output FITS> is name of output FITS file with product of dome and star flats.";

int
main(int argc, char *argv[])
{
  if (argc!=5) {
    cerr << usage << endl;
    exit(1);
  }
  try {
  string infile = argv[1];
  string photfile = argv[2];
  string mapname = argv[3];
  string outFits = argv[4];
  
  // Bounds of device pixel coordinates that we will try to map:
  Bounds<int> bDevice(16, 2033, 16, 4081);

  // Read in the device information
  map<string,Device> devices = decamInfo();

  Image<> flat;
  {
    FitsImage<> fi(infile, FITS::ReadOnly, 0);
    flat = fi.extract();
  }

  loadPhotoMapParser();
  PhotoMapCollection photomaps;
  {
    ifstream ifs(photfile.c_str());
    photomaps.read(ifs);
  }
  SubMap* devPhoto = photomaps.issueMap(mapname);

  photometry::PhotoArguments args;
  args.xExposure = args.yExposure = args.color = 0.;

  for (int iy = bDevice.getYMin(); iy <= bDevice.getYMax(); iy++) {
    for (int ix = bDevice.getXMin(); ix <= bDevice.getXMax(); ix++) {
      args.xDevice = ix;
      args.yDevice = iy;
      double photo = devPhoto->forward(0., args);
      flat(ix,iy) *= pow(10., 0.4*photo);
    }
  }
  FitsImage<>::writeToFITS(outFits, flat, 0);
  } catch (std::runtime_error& e) {
    quit(e,1);
  }
}
