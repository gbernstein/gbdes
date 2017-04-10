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
#include "DECamInfo.h"  // For obtaining CCDNUM

#ifdef _OPENMP
#include <omp.h>
#endif

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
  "          parameter file(s) specified on cmd line\n"
  "Outputs: 3 FITS files will be created, each having one extension per device for the\n"
  "   instrument.  <outputBase>.[xy].fits will hold xy pixel locations (or differentials\n"
  "   wrt TPV) while <outputBase>.a.fits holds relative pixel areas.\n"
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
  bool   makeFITS;
  int    tpvOrder;
  double color;
  int    decimate;
  double nominalPixel;
  double roundto;
  
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
    parameters.addMember("makeFITS",&makeFITS, def,
			 "Draw the FITS images?",true);
    parameters.addMember("decimate",&decimate, def+low,
			 "Decimation factor for output images",1,1);
    parameters.addMember("tpvOrder",&tpvOrder, def+low+up,
			 "Order of TPV polynomial",3,3,4);
    parameters.addMember("color",&color, def,
			 "Color assumed for source", 0.61);
    parameters.addMember("nominalPixel",&nominalPixel, def + lowopen,
			 "Reference pixel size (arcsec)", 0.264,0.);
    parameters.addMember("roundto",&roundto, def + low,
			 "Round values to this fraction of pixel (0=>exact)", 1e-4,0.);
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
    const double SCAMPTolerance=0.5 / 3600.;

    // Create empty images in primary extensions of output files
    if (makeFITS) {
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

    auto decammap = decam::decamInfo();
    
    for (auto device : deviceNames) {
      // Loop through devices
      cerr << "Working on " << device << endl;
      // Get the CCDNUM if the DETPOS matches a DECam detector
      auto diter = decammap.find(device);
      // Acquire the desired map
      auto pm = pmc.cloneMap(instrumentName + "/" + device);
      
      // Make images of x, y world coordinates

      Image<double> xout(drawBounds);
      Image<double> yout(drawBounds);
      Image<> area(drawBounds,1.);

      Bounds<int> mid(2,drawBounds.getXMax()-1,
		      2,drawBounds.getYMax()-1);

      if (makeFITS) {
	xout.header()->replace("DETPOS",device);
	yout.header()->replace("DETPOS",device);
	if (diter != decammap.end()) {
	  xout.header()->replace("CCDNUM",diter->second.ccdnum);
	  yout.header()->replace("CCDNUM",diter->second.ccdnum);
	}
      
#ifdef _OPENMP
#pragma omp parallel for
#endif
	for (int iy = drawBounds.getYMin(); iy<=drawBounds.getYMax(); iy++) {
	  double xw, yw;
	  double yp = (iy-1)*decimate + 1.;
	  for (int ix = drawBounds.getXMin(); ix<=drawBounds.getXMax(); ix++) {
	    double xp = (ix-1)*decimate + 1.;
	    pm->toWorld(xp,yp,xw,yw,color);
	    xout(ix,iy) = xw;
	    yout(ix,iy) = yw;
	  }
	}
      
	// Now draw the pixel area
	area.header()->replace("DETPOS",device);
	if (diter != decammap.end()) {
	  area.header()->replace("CCDNUM",diter->second.ccdnum);
	}
	// Conversion from differentials to area relative to reference pixel
	double rescale = pow( 2 * decimate * nominalPixel / 3600., -2.);
	// Store image area rounded,
#ifdef _OPENMP
#pragma omp parallel for
#endif
	for (int iy = mid.getYMin(); iy<=mid.getYMax(); iy++)
	  for (int ix = mid.getXMin(); ix<=mid.getXMax(); ix++) {
	    double a = (xout(ix+1,iy)-xout(ix-1,iy))*(yout(ix,iy+1)-yout(ix,iy-1)) -
	      (xout(ix,iy+1)-xout(ix,iy-1))*(yout(ix+1,iy)-yout(ix-1,iy));
	    if (roundto>0)
	      a = floor(a*rescale/roundto + 0.5) * roundto;
	    area(ix,iy) = a;
	  }
      } // End if (makeFITS)
      
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
	if (diter != decammap.end())
	  hdr.replace("CCDNUM",diter->second.ccdnum);

	// Write TPV header to text file
	*ofs << hdr;
	  
	if (makeFITS) {
	  // Place TPV into FITS Image headers
	  area.header()->copyFrom(hdr);
	  xout.header()->copyFrom(hdr);
	  yout.header()->copyFrom(hdr);

	  // Get difference btwn exact and TPV
	  Image<double> xtpv(drawBounds);
	  Image<double> ytpv(drawBounds);
	  xout.header()->replace("DETPOS",device);
	  yout.header()->replace("DETPOS",device);
#ifdef _OPENMP
#pragma omp parallel for
#endif
	  for (int iy = drawBounds.getYMin(); iy<=drawBounds.getYMax(); iy++) {
	    double xw, yw;
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

	   
	  // Convert back to pixel displacement, round, save as floats
	  Image<> dx(drawBounds,0.);
	  Image<> dy(drawBounds,0.);

	  for (int iy = mid.getYMin(); iy<=mid.getYMax(); iy++)
	    for (int ix = mid.getXMin(); ix<=mid.getXMax(); ix++) {
	      double dxdx = (xtpv(ix+1,iy)-xtpv(ix-1,iy)) / (2.*decimate);
	      double dydx = (ytpv(ix+1,iy)-ytpv(ix-1,iy)) / (2.*decimate);
	      double dxdy = (xtpv(ix,iy+1)-xtpv(ix,iy-1)) / (2.*decimate);
	      double dydy = (ytpv(ix,iy+1)-ytpv(ix,iy-1)) / (2.*decimate);
	      double det = dxdx*dydy - dxdy*dydx;
	      double d = (xout(ix,iy) * dydy - yout(ix,iy)*dxdy) / det;
	      if (roundto>0)
		d = floor(d/roundto+0.5) * roundto;
	      dx(ix,iy) = d;
	      d = (-xout(ix,iy) * dydx + yout(ix,iy)*dxdx) / det;
	      if (roundto>0)
		d = floor(d/roundto+0.5) * roundto;
	      dy(ix,iy) = d;
	    }
      
	  dx.header()->copyFrom(hdr);
	  dy.header()->copyFrom(hdr);
	  FitsImage<>::writeToFITS(xFITS, dx, device);
	  FitsImage<>::writeToFITS(yFITS, dy, device);
	}
      } else if (makeFITS) {
	// Save double-precision x and y maps
	FitsImage<double>::writeToFITS(xFITS, xout, device);
	FitsImage<double>::writeToFITS(yFITS, yout, device);
      }

      if (makeFITS) {
	// Write area map
	FitsImage<>::writeToFITS(aFITS, area, device);
      }

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
