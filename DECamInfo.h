// Some routines used frequently in working with DECam images
#ifndef DECAMINFO_H
#define DECAMINFO_H

#include <map>
#include "Bounds.h"
#include "Wcs.h"
#include "PixelMapCollection.h"

namespace decam {
  class Device {
  public:
    int ccdnum;
    Bounds<double> b;	// Bounds in the sky-field units
    int detsecX;		// Position of LL pixel in a close-packed image of array
    int detsecY;
    double norm;		// Normalization applied to flat field
  };

  // Return a DETPOS-indexed container of Devices - unit norms
  std::map<string,Device> decamInfo();

  // Load a set of normalizations from a file into the Device map 
  void getDeviceNorms(string filename, std::map<string,Device>& devices);

  // Class that will map between field coordinates (in degrees) and pixel coords,
  // using a particular fixed astrometric solution.
  class DECamMap {
  public:
    DECamMap();
    // Choose device to which the pixel coords apply
    void setDevice(string detpos);
    void toPix(double xField, double yField, double& xPix, double& yPix) const;
    void toField(double xPix, double yPix, double& xField, double& yField) const;
  private:
    DECamMap(const DECamMap& rhs);	// hide
    void operator=(const DECamMap& rhs); // hide
    const string prefix;
    astrometry::PixelMapCollection pmc;
    double centerX;
    double centerY;
    astrometry::Wcs* wcs;
    static string referenceMap;
  };

} // end namespace decam
#endif
