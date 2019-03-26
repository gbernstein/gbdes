// Make a FITS image of a DECam photometric solution, with flat normalization factors taken out.
#include <fstream>
#include <sstream>
#include <map>
#include <ctime>
#include <iomanip>

#include "StringStuff.h"
#include "Statistics.h"

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
  "usage: Photo2DESDM <photo file> <instrument name> <output FITS base> [flatnorms] \n"
  "     <photo file> is file holding the serialized solution from PhotoFit\n"
  "     <instrument name> is the filter name, should be instrument name in photo file too\n"
  "     <output FITS base> is root of the files that will be drawn, one per CCD\n"
  "     [flatnorms] is an optional file with CCDPOS in 1st column and 2nd column containing\n"
  "              a normalization factor to divide into each CCD's star flat by.";

int
main(int argc, char *argv[])
{
  if (argc<4 || argc>5) {
    cerr << usage << endl;
    exit(1);
  }
  string photoFile = argv[1];
  string filter = argv[2];
  string photoPrefix = filter + "/";
  string outFits = argv[3];
  bool useNorms = argc>4;
  string normFile = useNorms ? argv[4] : "";

  bool oneFile = false;	// set true to output all CCDs into one multi-extension FITS
  try {
  loadPhotoMapParser();

  // Size of full devices:
  Bounds<int> bpix(1, 2048, 1, 4096);
  // Bounds of device pixel coordinates that we will use photo solutions.
  // They will just be unity outside these bounds.
  Bounds<int> bDevice(11, 2038, 1, 4096);

  // Read in the device information
  map<string,Device> devices = decamInfo();

  if (useNorms) {
    // Read in normalization factors
    cerr << "Reading device norms from " << normFile << endl;
    getDeviceNorms(normFile, devices);
  }

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
    if (useNorms) cerr << "Using norm " << i->second.norm << endl;

    photometry::SubMap* devPhoto=0;
    try {
      devPhoto = photomaps.issueMap(photoPrefix + i->first);
    } catch (photometry::PhotometryError& m) {
      cerr << "Missing, skipping" << endl;
      devPhoto = 0;
    }
    if (!devPhoto) continue;

    int ccdNum = i->second.ccdnum;
    Image<> starflat(bpix, 1.);

    starflat.header()->append("CCDNUM", ccdNum);
    starflat.header()->append("DETPOS", i->first);
    starflat.header()->append("FILTER", filter);
    starflat.header()->append("BAND", filter);
    time_t now;
    time(&now);
    string date = ctime(&now);
    stripWhite(date);  // Remove trailing LF/CR's.
    starflat.header()->append("DESMKICR", date);
    string s = "illumcor";
    starflat.header()->append("OBSTYPE", s);
    s = "IMAGE ";
    starflat.header()->append("DES_EXT", s);
    s = "SCI";
    starflat.header()->append("EXTNAME", s);
    starflat.header()->addHistory("Created from starflat solution " + photoFile);

    // Loop over pixels
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for (int iy = bDevice.getYMin(); iy <= bDevice.getYMax(); iy++) {
      photometry::PhotoArguments args;
      args.xExposure = args.yExposure = args.color = 0.;
      for (int ix = bDevice.getXMin(); ix <= bDevice.getXMax(); ix++) {
	args.xDevice = ix;
	args.yDevice = iy;
	double photo = devPhoto->forward(0., args);
	starflat(ix, iy) = pow(10., 0.4*photo);
      }
    }
    if (useNorms) starflat *= i->second.norm;

    // Calculate the median value of each amplifier and store in header
    // Subsample the pixels in the useful image region at 4x4.
    std::vector<float> v;
    Bounds<int> b = decam::datasec(i->first,"A") & bDevice;
    for (int iy = b.getYMin(); iy <= b.getYMax(); iy+=4)
      for (int ix = b.getXMin(); ix <= b.getXMax(); ix+=4)
	v.push_back(starflat(ix,iy));
    starflat.header()->append("FLATMEDA", stats::median(v));
    v.clear();
    b = decam::datasec(i->first,"B") & bDevice;
    for (int iy = b.getYMin(); iy <= b.getYMax(); iy+=4)
      for (int ix = b.getXMin(); ix <= b.getXMax(); ix+=4)
	v.push_back(starflat(ix,iy));
    starflat.header()->append("FLATMEDB", stats::median(v));

    // Add fake WCS to enable array display with DS9 (as per Yanny message 24 Jan 2014)
    starflat.header()->append("CRVAL1",0.);
    starflat.header()->append("CRVAL2",0.);
    starflat.header()->append("CTYPE1",string("RA---TAN"));
    starflat.header()->append("CTYPE2",string("DEC--TAN"));
    starflat.header()->append("CD1_1",0.);
    starflat.header()->append("CD1_2",7.286e-5);
    starflat.header()->append("CD2_2",0.);
    starflat.header()->append("CD2_1",-7.286e-5);
    int ny = (i->second.detsecY-1)    / 2048;
    int nx = (i->second.detsecX-2049) / 2048;
    starflat.header()->append("CRPIX1",13423.2 - nx*2254.4);
    starflat.header()->append("CRPIX2",14826. - ny*2129.666667);

    if (oneFile) {
      FitsImage<>::writeToFITS(outFits + ".fits", starflat, i->first);
    } else {
      ostringstream oss;
      oss << outFits << "_" << setfill('0') << setw(2) << ccdNum << ".fits";
      string fname = oss.str();
      FitsImage<>::writeToFITS(fname, starflat, 0);
    }
  } // end Device loop.

  {
    // And write a dummy image for dead CCD 61:
    int ccdNum = 61;
    string detpos = "N30";
    Image<> starflat(bpix, 1.);

    starflat.header()->append("CCDNUM", ccdNum);
    starflat.header()->append("DETPOS", detpos);
    starflat.header()->append("FILTER", filter);
    time_t now;
    time(&now);
    string date = ctime(&now);
    starflat.header()->append("DESMKICR", date);
    string s = "illumcor";
    starflat.header()->append("OBSTYPE", s);
    s = "IMAGE ";
    starflat.header()->append("DES_EXT", s);
    starflat.header()->addHistory("Dummy image for dead CCD ");

    // Add fake WCS to enable array display with DS9 (as per Yanny message 24 Jan 2014)
    starflat.header()->append("CRVAL1",0.);
    starflat.header()->append("CRVAL2",0.);
    starflat.header()->append("CTYPE1",string("RA---TAN"));
    starflat.header()->append("CTYPE2",string("DEC--TAN"));
    starflat.header()->append("CD1_1",0.);
    starflat.header()->append("CD1_2",7.286e-5);
    starflat.header()->append("CD2_2",0.);
    starflat.header()->append("CD2_1",-7.286e-5);
    int detsecX = 24577;
    int detsecY = 12289;
    int ny = (detsecY-1)    / 2048;
    int nx = (detsecX-2049) / 2048;
    starflat.header()->append("CRPIX1",13423.2 - nx*2254.4);
    starflat.header()->append("CRPIX2",14826. - ny*2129.666667);

    if (oneFile) {
      FitsImage<>::writeToFITS(outFits + ".fits", starflat, detpos);
    } else {
      ostringstream oss;
      oss << outFits << "_" << setfill('0') << setw(2) << ccdNum << ".fits";
      string fname = oss.str();
      FitsImage<>::writeToFITS(fname, starflat, 0);
    }
  }

  } catch (std::runtime_error& e) {
    quit(e,1);
  }
}
