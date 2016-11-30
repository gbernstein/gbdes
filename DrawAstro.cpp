// Make FITS images of an astrometric solution
#include <fstream>
#include <sstream>
#include <map>
#include "Std.h"
#include "StringStuff.h"

#include "FitsImage.h"
#include "FitSubroutines.h"

using namespace img;
using namespace astrometry;
using namespace stringstuff;

string usage = 
  "DrawAstro: Make a FITS image sampling an astrometric map.\n"
  "usage: DrawAstro <YAML file> <map name> <output file>\n"
  "...add options for bounds/step of pixel values; using FITS/TMV";

int
main(int argc, char *argv[])
{
  if (argc!=4) {
    cerr << usage << endl;
    exit(1);
  }
  string astroFile = argv[1];
  string mapName = argv[2];
  string fitsName = argv[3];
  
  int xStart = 1;
  int xEnd = 2049;
  int yStart = 1;
  int yEnd = 4097;
  int xyStep = 8;
  int dxy = 5;  // Step in pixels used for finite-difference derivs

  double color;  // Stellar color to use for maps
    
  try {
    loadPixelMapParser();

    // Get the desired pixel map
    PixelMapCollection pmc;
    {
      ifstream ifs(astroFile.c_str());
      pmc.read(ifs);
    }
    auto pm = pmc.issueMap(mapName);

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

    exit(0);
  } catch (std::runtime_error &m) {
    quit(m,1);
  }
}
