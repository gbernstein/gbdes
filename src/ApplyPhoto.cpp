#include <fstream>
#include <sstream>
#include <map>

#include "Std.h"
#include "StringStuff.h"
#include "Astrometry.h"
#include "Header.h"
#include "TPVMap.h"

#include "FitSubroutines.h"

using namespace std;

string usage = "ApplyPhoto: transform magnitudes using solutions from PhotoFit\n"
               "usage: ApplyPhoto <PhotoIn> [WCSin RA Dec]\n"
               "      <PhotoIn> is photometric solution to use, in form of\n"
               "               <mapname>@<photo file>, <photo file> holds serialized PhotoFit solution\n"
               "      [WCSin RA Dec] These three arguments are needed only if the photometric map is\n"
               "         dependent upon the position of object in the exposure FOV.  In this case:\n"
               "      <WCSin> is astrometric solution to use, in form of either\n"
               "              name of file holding FITS-style WCS spec, or\n"
               "              <wcsname>@<filename>, reads WCS with given name serialized in the file\n"
               "      <RA, Dec> give the axis about which field coords were calculated\n"
               "             for this exposure by PhotoFit\n"
               "      stdin: each line is <mag In> <xpix> <ypix> [color], which are \n"
               "             input magnitude, x and y pixel coordinates, and optionally a color\n"
               "             to use in the transformation.  Otherwise color is set to 0.\n"
               "      stdout: each line is <mag out>";

int main(int argc, char *argv[])
{
  try
  {
    // Read parameters
    if (!(argc == 2 || argc == 5))
    {
      cerr << usage << endl;
      exit(1);
    }

    loadPixelMapParser();
    loadPhotoMapParser();

    // Deserialize the photometric transformation from a file
    photometry::PhotoMap *photo = 0;
    {
      string photofile = argv[1];
      size_t at = photofile.find('@');
      string mapname = photofile.substr(0, at);
      at++;
      if (at >= photofile.size())
      {
        cerr << "Missing filename in mapname@photofile argument: " << photofile << endl;
        exit(1);
      }
      string filename = photofile.substr(at);
      ifstream ifs(filename.c_str());
      if (!ifs)
      {
        cerr << "Could not open serialized Photo file " << filename << endl;
        exit(1);
      }
      photometry::PhotoMapCollection pmc;
      if (!pmc.read(ifs))
      {
        cerr << "File <" << filename << "> is not serialized Photo file" << endl;
        exit(1);
      }
      photo = pmc.cloneMap(mapname);
    }

    // Get the WCS
    astrometry::Wcs *wcs = nullptr;

    if (argc > 2)
    {
      string wcsfile = argv[2];

      if (wcsfile.find('@') != string::npos)
      {
        // Deserialize WCS from a file
        size_t at = wcsfile.find('@');
        string wcsname = wcsfile.substr(0, at);
        at++;
        if (at >= wcsfile.size())
        {
          cerr << "Missing filename in wcsname@filename argument: " << wcsfile << endl;
          exit(1);
        }
        string filename = wcsfile.substr(at);
        ifstream ifs(filename.c_str());
        if (!ifs)
        {
          cerr << "Could not open serialized WCS file " << filename << endl;
          exit(1);
        }
        astrometry::PixelMapCollection pmc;
        if (!pmc.read(ifs))
        {
          cerr << "File <" << filename << "> is not serialized WCS file" << endl;
          exit(1);
        }
        wcs = pmc.cloneWcs(wcsname);
      }
      else
      {
        // Get WCS as TPV in FITS-header-style file
        ifstream mapfs(wcsfile.c_str());
        if (!mapfs)
        {
          cerr << "Could not open TPV spec file " << wcsfile << endl;
          exit(1);
        }
        img::Header h;
        mapfs >> h;
        wcs = astrometry::readTPV(h);
        if (!wcs)
        {
          cerr << "Error reading input TPV map from file " << wcsfile << endl;
          exit(1);
        }
      }

      double ra = astrometry::hmsdeg(argv[2]) * DEGREE;
      double dec = astrometry::dmsdeg(argv[2]) * DEGREE;
      astrometry::Gnomonic exposure(astrometry::Orientation(astrometry::SphericalICRS(ra, dec)));
      wcs->reprojectTo(exposure);
    }

    // Read data
    string buffer;
    while (stringstuff::getlineNoComment(cin, buffer))
    {
      double magIn;
      photometry::PhotoArguments args;
      istringstream iss(buffer);
      if (!(iss >> magIn >> args.xDevice >> args.yDevice))
      {
        cerr << "Bad input entry <" << buffer << ">"
             << endl;
        exit(1);
      }
      if (!(iss >> args.color))
        args.color = 0.;
      // ??? could distinguish failure from EOF in the above...

      if (wcs)
        wcs->toWorld(args.xDevice, args.yDevice, args.xExposure, args.yExposure);
      else
        args.xExposure = args.yExposure = 0.;

      cout << photo->forward(magIn, args) << endl;
    }

    // Clean up
    if (wcs)
      delete wcs;
    delete photo;
  }
  catch (std::runtime_error &m)
  {
    quit(m, 1);
  }
}
