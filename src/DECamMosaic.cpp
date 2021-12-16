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
    "DECamMosaic: Take as input a multi-extension FITS file with one DECam CCD\n"
    "     image per extension, and produce a single output image with these values\n"
    "     projected onto the sky for the full field of view.\n"
    "usage: DECamMosaic <input FITS> <output FITS> [out pix scale] [fill value] [med box]\n"
    "     <input FITS> is name of multi-extension input FITS file\n"
    "     <output FITS> is name of output FITS file portraying array on sky.\n"
    "     [out pix scale] is number of DECam pixels (0.264\") per output pixel (default=32)\n"
    "     [fill value] is assigned to pixels with no CCD (default=-1)\n"
    "     [med box] is 'radius' of box around each sample point to take median (default=scale)";

// ??? Make an output WCS that works??
int main(int argc, char *argv[])
{
  if (argc < 3 || argc > 6)
  {
    cerr << usage << endl;
    exit(1);
  }
  try
  {
    string inFits = argv[1];
    string outFits = argv[2];
    double pixelScale = argc > 3 ? atof(argv[3]) : 32;
    const float noData = argc > 4 ? atof(argv[4]) : -1.;
    int medbox = argc > 5 ? atoi(argv[5]) : static_cast<int>(floor(pixelScale));

    // Convert rough input scale in DECam pixels into degrees per output pix:
    pixelScale *= 0.264 / 3600.;

    // Bounds of device pixel coordinates that we will try to map:
    Bounds<int> bDevice(16, 2033, 16, 4081);

    // Read in the device information
    map<string, Device> devices = decamInfo();

    DECamMap wcs;

    // Make bounding rectangle for whole exposure
    Bounds<double> bExpo;
    for (auto &i : devices)
      bExpo += i.second.b;

    // Make the Image
    Bounds<int> bImage(static_cast<int>(floor(bExpo.getXMin() / pixelScale)),
                       static_cast<int>(ceil(bExpo.getXMax() / pixelScale)),
                       static_cast<int>(floor(bExpo.getYMin() / pixelScale)),
                       static_cast<int>(ceil(bExpo.getYMax() / pixelScale)));
    Image<> big(bImage, noData);

    // Open the FITS file to count extensions
    int nHDUs;
    {
      FITS::FitsFile ff(inFits);
      nHDUs = ff.HDUCount();
    }

    for (int iHDU = 1; iHDU < nHDUs; iHDU++)
    {
      FitsImage<> fi(inFits, FITS::ReadOnly, iHDU);
      string detpos;
      if (!fi.header()->getValue("DETPOS", detpos))
      {
        cerr << "Did not find DETPOS at extension " << iHDU << endl;
        exit(1);
      }
      stripWhite(detpos);

      /**/ cerr << "HDU " << iHDU << " DETPOS " << detpos << endl;

      // Move along if this is an unwanted device
      if (devices.count(detpos) == 0)
        continue;

      wcs.setDevice(detpos);

      // Convert device boundaries into pixel range
      Bounds<double> b = devices[detpos].b;
      Bounds<int> bpix(static_cast<int>(floor(b.getXMin() / pixelScale)),
                       static_cast<int>(ceil(b.getXMax() / pixelScale)),
                       static_cast<int>(floor(b.getYMin() / pixelScale)),
                       static_cast<int>(ceil(b.getYMax() / pixelScale)));

      // Open the flat-field image for this extension
      Image<> flat = fi.extract();

      // Loop over pixels
      for (int iy = bpix.getYMin(); iy <= bpix.getYMax(); iy++)
        for (int ix = bpix.getXMin(); ix <= bpix.getXMax(); ix++)
        {
          double xField = ix * pixelScale;
          double yField = iy * pixelScale;
          double xPix, yPix;
          wcs.toPix(xField, yField, xPix, yPix);
          Bounds<int> medbounds(static_cast<int>(floor(xPix - medbox)),
                                static_cast<int>(ceil(xPix + medbox)),
                                static_cast<int>(floor(yPix - medbox)),
                                static_cast<int>(ceil(yPix + medbox)));
          // Overlap with measurable region
          medbounds = medbounds & bDevice;
          if (!medbounds)
            continue; // Skip if outside device

          // Get median of pixels in the box
          vector<float> values;
          for (int jy = medbounds.getYMin(); jy <= medbounds.getYMax(); jy++)
            for (int jx = medbounds.getXMin(); jx <= medbounds.getXMax(); jx++)
              if (flat(jx, jy) > 0)
                values.push_back(flat(jx, jy));

          // Apply normalization correction
          double median = stats::median(values);
          if (median <= 0)
            cerr << "median " << median << " for range " << medbounds << endl;
          big(ix, iy) = -2.5 * log10(median);
        }
    } // end Device loop.
    big.shift(1, 1);
    FitsImage<>::writeToFITS(outFits, big, 0);
  }
  catch (std::runtime_error &e)
  {
    quit(e, 1);
  }
}
