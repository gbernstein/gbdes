// Make FITS images of an astrometric solution
#include <fstream>
#include <sstream>
#include <map>
#include "Std.h"
#include "StringStuff.h"
#include "Pset.h"

#include "FitsImage.h"
#include "FitSubroutines.h"
#include "TPVMap.h"

using namespace img;
using namespace astrometry;
using namespace stringstuff;

string usage = 
  "DrawAstro: Make FITS images of astrometric maps for all devices\n"
  "of a selected instrument.  Also optionally fit TPV map to them\n"
  "and output differential to that instead.\n"
  "\n"
  "usage: DrawAstro <YAML file> <instrument name> [parameter file] [parameter file...]\n"
  "   [-parameter[=]value...]\n"
  "      <YAML file>:  file containing serialized WCS/Pixel maps\n"
  "      <instrument name>: base name of maps to tabulate, find all devices for it.\n"
  "      Program parameters specified as command-line options or read from\n"
  "          parameter file(s) specified on cmd line"
  "Outputs: 3 FITS files will be created, each having one extension per device for the\n"
  "   instrument.  <outputBase>.[xy].fits will hold xy pixel locations (or differentials wrt TPV)\n"
  "   while <outputBase>.a.fits holds relative pixel areas.\n"
  "   If makeTPV=true, another file <outputBase>.tpv will hold ASCII version of all TPV's, which\n"
  "   will also be placed into each image's header.";

// ??? When decimate != 1, the TPV headers will not describe the images they're attached to.
// Could change the CRPIX and CD values to adjust this for proper displaying.

int
main(int argc, char *argv[])
{
  string astroFile;
  string instrumentName;
  string outputBase;
  bool   makeTPV;
  int    tpvOrder;
  double color;
  int    decimate;
  double nominalPixel;
  
  Pset parameters;
  {
    const int def=PsetMember::hasDefault;
    const int low=PsetMember::hasLowerBound;
    const int up=PsetMember::hasUpperBound;
    const int lowopen = low | PsetMember::openLowerBound;
    const int upopen = up | PsetMember::openUpperBound;

    parameters.addMember("outputBase",&outputBase, def,
			 "base of filenames for output FITS files","astrom");
    parameters.addMember("makeTPV",&makeTPV, def,
			 "Set to fit TPV file and return pixel displacements relative to it",false);
    parameters.addMember("decimate",&decimate, def+low,
			 "Decimation factor for output images",1,1);
    parameters.addMember("tpvOrder",&tpvOrder, def+low+up,
			 "Order of TPV polynomial",3,3,4);
    parameters.addMember("color",&color, def,
			 "Color assumed for source", 0.61);
    parameters.addMember("nominalPixel",&nominalPixel, def + lowopen,
			 "Reference pixel size (arcsec)", 0.264,0.);
  }

  try {
    
    // Read all the command-line and parameter-file program parameters
    processParameters(parameters, usage, 2, argc, argv);

    astroFile = argv[1];
    instrumentName = argv[2];

    // Names we'll use for output images of x,y coordinates and pixel areas
    string xFITS = outputBase + ".x.fits";
    string yFITS = outputBase + ".y.fits";
    string aFITS = outputBase + ".a.fits";
    // And if we want to output text TPV
    string tpvFile = outputBase + ".tpv";

    // tolerance on SCAMP-format fit to real solution:
    const double SCAMPTolerance=0.0001 / 3600.;

    // Create empty images in primary extensions of output files
    {
      FitsImage<> fx(xFITS,FITS::OverwriteFile + FITS::Create, 0);
      FitsImage<> fy(yFITS,FITS::OverwriteFile + FITS::Create, 0);
      FitsImage<> fa(aFITS,FITS::OverwriteFile + FITS::Create, 0);
    }
    // Open the TPV file if we're making one
    ofstream* ofs=0;
    if (makeTPV)
      ofs = new ofstream(tpvFile.c_str());
    
    // (inclusive) bounds of pixel numbers
    Bounds<int> ccdBounds(1,2048,1,4096);
    // Size of images to draw
    Bounds<int> drawBounds(1, ccdBounds.getXMax()/decimate,
			   1, ccdBounds.getYMax()/decimate);
    
    loadPixelMapParser();

    // Read the pixel maps
    PixelMapCollection pmc;
    {
      ifstream ifs(astroFile.c_str());
      pmc.read(ifs);
    }

    // Find all maps that look like devices of the instrument
    auto mapNames = findMatches("^" + instrumentName + "/[NS][[:digit:]]*$",
				pmc.allMapNames());
    
    set<string> deviceNames;
    for (auto mapName : mapNames) {
      deviceNames.insert(split(mapName,'/').back());
    }

    for (auto device : deviceNames) {
      // Loop through devices
      cerr << "Working on " << device << endl;
      // Acquire the desired map
      auto pm = pmc.cloneMap(instrumentName + "/" + device);
      
      // Make images of x, y world coordinates

      Image<double> xout(drawBounds);
      Image<double> yout(drawBounds);
      xout.header()->replace("DETPOS",device);
      yout.header()->replace("DETPOS",device);
      double xw, yw;
      for (int iy = drawBounds.getYMin(); iy<=drawBounds.getYMax(); iy++) {
	double yp = (iy-1)*decimate + 1.;
	for (int ix = drawBounds.getXMin(); ix<=drawBounds.getXMax(); ix++) {
	  double xp = (ix-1)*decimate + 1.;
	  pm->toWorld(xp,yp,xw,yw,color);
	  xout(ix,iy) = xw;
	  yout(ix,iy) = yw;
	}
      }
      // Calculate pixel area map - edges will default to unity.
      int nx = drawBounds.getXMax();
      int ny = drawBounds.getYMax();
      Bounds<int> mid(2,nx-1,2,ny-1);
      
      Image<> area(drawBounds,1.);
      area.header()->replace("DETPOS",device);
      // Conversion from differentials to area relative to reference pixel
      double rescale = pow( 2 * decimate * nominalPixel / 3600., -2.);
      // Store image area rounded to accuracy 1 part in 10^5
      double round = 1e5;
      for (int iy = mid.getYMin(); iy<=mid.getYMax(); iy++)
	for (int ix = mid.getXMin(); ix<=mid.getXMax(); ix++) {
	  double a = (xout(ix+1,iy)-xout(ix-1,iy))*(yout(ix,iy+1)-yout(ix,iy-1)) -
	    (xout(ix,iy+1)-xout(ix,iy-1))*(yout(ix+1,iy)-yout(ix-1,iy));
	  area(ix,iy) = floor(a*rescale*round + 0.5) / round;
	}
      
      if (makeTPV) {
	// Fit TPV solution to this (add projection to 00+00)
	Orientation orient(SphericalICRS(0.,0.));
	Gnomonic projection(0., 0., orient);
	Wcs devWcs(pm, projection);
	double buffer=30.;  // Number of pixels to ignore at edges when fitting TPV
	Wcs* tpv = fitTPV(Bounds<double>(ccdBounds.getXMin()+buffer,
					 ccdBounds.getXMax()-buffer,
					 ccdBounds.getYMin()+buffer,
					 ccdBounds.getYMax()-buffer),
			  devWcs, projection,
			  "NoName", color, SCAMPTolerance, tpvOrder);
	auto hdr = writeTPV(*tpv);
	hdr.replace("DETPOS",device);

	// Write TPV header to text file
	*ofs << hdr;
	  
	// Place TPV into FITS Image headers
	area.header()->copyFrom(hdr);
	xout.header()->copyFrom(hdr);
	yout.header()->copyFrom(hdr);

	// Get difference btwn exact and TPV
	Image<double> xtpv(drawBounds);
	Image<double> ytpv(drawBounds);
	xout.header()->replace("DETPOS",device);
	yout.header()->replace("DETPOS",device);
	double xw, yw;
	for (int iy = drawBounds.getYMin(); iy<=drawBounds.getYMax(); iy++) {
	  double yp = (iy-1)*decimate + 1.;
	  for (int ix = drawBounds.getXMin(); ix<=drawBounds.getXMax(); ix++) {
	    double xp = (ix-1)*decimate + 1.;
	    tpv->getMap()->toWorld(xp,yp,xw,yw,color);
	    xtpv(ix,iy) = xw;
	    ytpv(ix,iy) = yw;
	  }
	}
	xout -= xtpv;
	yout -= ytpv;
	   
	// Convert back to pixel displacement, round to 0.0001 pixel, save as floats
	Image<> dx(drawBounds,0.);
	Image<> dy(drawBounds,0.);
	double roundp = 1e4;

	for (int iy = mid.getYMin(); iy<=mid.getYMax(); iy++)
	  for (int ix = mid.getXMin(); ix<=mid.getXMax(); ix++) {
	    double dxdx = (xtpv(ix+1,iy)-xtpv(ix-1,iy)) / (2.*decimate);
	    double dydx = (ytpv(ix+1,iy)-ytpv(ix-1,iy)) / (2.*decimate);
	    double dxdy = (xtpv(ix,iy+1)-xtpv(ix,iy-1)) / (2.*decimate);
	    double dydy = (ytpv(ix,iy+1)-ytpv(ix,iy-1)) / (2.*decimate);
	    double det = dxdx*dydy - dxdy*dydx;
	    double d = (xout(ix,iy) * dydy - yout(ix,iy)*dxdy) / det;
	    dx(ix,iy) = floor(d*roundp+0.5) / roundp;
	    d = (-xout(ix,iy) * dydx + yout(ix,iy)*dxdx) / det;
	    dy(ix,iy) = floor(d*roundp+0.5) / roundp;
	  }
      
	dx.header()->copyFrom(hdr);
	dy.header()->copyFrom(hdr);
	FitsImage<>::writeToFITS(xFITS, dx, device);
	FitsImage<>::writeToFITS(yFITS, dy, device);

      } else {
	// Save double-precision x and y maps
	FitsImage<double>::writeToFITS(xFITS, xout, device);
	FitsImage<double>::writeToFITS(yFITS, yout, device);
      }

      // Write x, y, area maps
      FitsImage<>::writeToFITS(aFITS, area, device);
      delete pm;
    }

    // close output text stream if any
    if (ofs)
      ofs->close();

    exit(0);
  } catch (std::runtime_error &m) {
    quit(m,1);
  }
}
      
  /**


    // Make target images
    int nx = (xEnd-xStart)/xyStep;
    int ny = (yEnd-yStart)/xyStep;
    Bounds<int> bFITS(1,nx,1,ny);
    Image<double> xw(bFITS,0.);
    Image<double> yw(bFITS,0.);
    Image<> pixArea(bFITS,0.);
    Image<> shear1(bFITS,0.);
    Image<> shear2(bFITS,0.);
    Image<> rotation(bFITS,0.);

    // Start filling, converting to degrees / arcsec as needed
    double tx,ty;
    Matrix22 dwdp;
    for (int iy=1; iy<=ny; iy++) {
      double yp = yStart + (iy-1)*xyStep;
      for (int ix=1; ix<=nx; ix++) {
	double xp = xStart + (ix-1)*xyStep;
	pm->toWorld(xp,yp,tx,ty,color);
	dwdp = pm->dWorlddPix(xp,yp,color);
	xw(ix,iy) = tx;
	yw(ix,iy) = ty;
	double det = dwdp(0,0)*dwdp(1,1)-dwdp(1,0)*dwdp(0,1);
	if (det<0.) {
	  // Flip parity to get shear/rotation right
	  dwdp(0,0) *= -1.;
	  dwdp(0,1) *= -1.;
	}
	pixArea(ix,iy) = det * (3600.*3600.); // Give in arcsec^2
	shear1(ix,iy) = (dwdp(0,0)-dwdp(1,1)) / sqrt(det);
	shear2(ix,iy) = (dwdp(1,0)+dwdp(0,1)) / sqrt(det);
	rotation(ix,iy) = (dwdp(1,0)-dwdp(0,1)) / sqrt(det);
      }
    }

    // Write to output
    FitsImage<double>::writeToFITS(fitsName, xw, "X");
    FitsImage<double>::writeToFITS(fitsName, yw, "Y");
    FitsImage<>::writeToFITS(fitsName, pixArea, "AREA");
    FitsImage<>::writeToFITS(fitsName, shear1, "SHEAR1");
    FitsImage<>::writeToFITS(fitsName, shear2, "SHEAR2");
    FitsImage<>::writeToFITS(fitsName, rotation, "ROTATION");

}
  ***/
