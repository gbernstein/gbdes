// Code to read PhotoPriors from files.
// Used in PhotoFit.cpp.

// Program to fit photometric solutions to detection catalogs already matched by WCSFoF.
#include <fstream>
#include <sstream>
#include <map>
#include <iomanip>

#include "Std.h"
#include "Astrometry.h"
#include "StringStuff.h"
#include "PixelMapCollection.h"
#include "PhotoMatch.h"
#include "PhotoMapCollection.h"
#include "PhotoInstrument.h"

using namespace std;
using namespace stringstuff;
using namespace photometry;

void formatError(string filename, string line) {
  cerr << "Format error reading prior file " << filename << " on input line:" << endl;
  cerr << line << endl;
  exit(1);
}

// Read double value and then possible "fixed" or "free" from stream.
// Missing fixed/free defaults to fixed.
// Return true on success.
bool valueAndFree(istream& is, double& value, bool& free) {
  if (!(is >> value)) return false;
  string freeString;
  if (!(is >> freeString)) {
    free = false;
    return true;
  } else if (nocaseEqual(freeString, "fixed")) {
    free = false;
    return true;
  } else if (nocaseEqual(freeString, "free")) {
    free = true;
    return true;
  } else {
    // Unknown word:
    return false;
  }
}

list<PhotoPrior*>
readPriors(string filename, 
	   const vector<Instrument*>& instruments, 
	   const vector<Exposure*>& exposures, 
	   const vector<Extension*>& extensions, 
	   const vector<long>& detectionsPerExposure) {

  Assert(detectionsPerExposure.size() == exposures.size());

  ifstream ifs(filename.c_str());
  if (!ifs) {
    cerr << "Could not open photometric priors file " << filename << endl;
    exit(1);
  }

  list<PhotoPrior*> out;

  string buffer;
  while (getlineNoComment(ifs, buffer)) {
    string name;
    double sigma;
    {
      istringstream iss(buffer);
      if (!(iss >> name >> sigma) || sigma<=0.) {
	cerr << "Error reading prior name & sigma at line: " << buffer << endl;
	exit(1);
      }
    }

    list<PhotoPriorReferencePoint> refs;
    double zeropoint=0.;
    double airmassCoeff=0.;
    double colorCoeff=0.;
    bool zpIsFree = false;
    bool airmassIsFree = false;
    bool colorIsFree = false;

    while (getlineNoComment(ifs,buffer)) {
      string keyword;
      double value;
      string fixedString;
      istringstream iss(buffer);
      if (!(iss >> keyword)) 
	formatError(filename, buffer);

      // "END" marks end of a prior
      if (nocaseEqual(keyword, "END")) 
	break;
      
      if (nocaseEqual(keyword, "zeropoint")) {
	if (!valueAndFree(iss, zeropoint, zpIsFree))
	  formatError(filename, buffer);
      } else if (nocaseEqual(keyword, "airmass")) {
	if (!valueAndFree(iss, airmassCoeff, airmassIsFree))
	  formatError(filename, buffer);
      } else if (nocaseEqual(keyword, "color")) {
	if (!valueAndFree(iss, colorCoeff, colorIsFree))
	  formatError(filename, buffer);
      } else {
	// The keyword should be the name of an exposure.
	string exposureRegex = keyword;
	// We will also get name of reference device,
	string deviceName;
	// and pixel coords of reference point - device center is default
	bool useCenterOfDevice = true;
	double xpix, ypix;
	// and reference instrumental magnitude
	double instrumentalMag = 0.;

	if (!(iss >> deviceName))
	  formatError(filename, buffer);
	if (iss >> xpix >> ypix) {
	  useCenterOfDevice = false;
	  iss >> instrumentalMag;
	}

	// Now look for all exposures that match the request here.
	vector<int> matchingExposures;
	for (int iExpo = 0; iExpo < exposures.size(); iExpo++) {
	  if (!exposures[iExpo]) continue;
	  if (detectionsPerExposure[iExpo] <=0) continue; // Exclude no-data exposures
	  if (regexMatch(exposureRegex, exposures[iExpo]->name))
	    matchingExposures.push_back(iExpo);
	}

	// For each matching exposure: look for instrument, device name, 
	// and get center coords if needed
	vector<int> matchingDevices;
	vector<double> matchingXPix;
	vector<double> matchingYPix;
	for (vector<int>::iterator i = matchingExposures.begin();
	     i != matchingExposures.end(); ) {
	  int iExpo = *i;
	  Assert(exposures[iExpo]);
	  int iInst = exposures[iExpo]->instrument;
	  if (iInst<0 || !instruments[iInst]) {
	    // We will ignore failed exposure specifications and only warn when
	    // nothing turns up valid at the end.
	    i = matchingExposures.erase(i);
	    continue;
	  }

	  Instrument& inst = *instruments[iInst];
	  int iDevice = inst.deviceNames.indexOf(deviceName);
	  if (iDevice < 0) {
	    // This instrument has no such device.
	    i = matchingExposures.erase(i);
	    continue;
	  }
	  // Save away the device's number and the pixel coordinates we will use.
	  matchingDevices.push_back(iDevice);
	  Position<double> center = inst.domains[iDevice].center();
	  matchingXPix.push_back( useCenterOfDevice ? center.x : xpix);
	  matchingYPix.push_back( useCenterOfDevice ? center.y : ypix);

	  ++i;
	} // end of loop finding matching exposures
	    
	// Troll extensions looking for exposure + device match
	bool foundOne = false;
	for (int iExtn = 0; iExtn < extensions.size(); iExtn++) {
	  if (!extensions[iExtn]) continue;
	  Extension& extn = *extensions[iExtn];
	  bool extnMatch = false;
	  for (int iMatch = 0; !extnMatch && iMatch < matchingExposures.size(); iMatch++) {
	    if (extn.exposure==matchingExposures[iMatch]
		&& extn.device==matchingDevices[iMatch]) {
	      extnMatch = true; // Stop looping - a given extension will match only once
	      PhotoPriorReferencePoint point;
	      point.exposureName = exposures[matchingExposures[iMatch]]->name;
	      point.deviceName = deviceName;
	      point.magIn = instrumentalMag;
	      // Check for airmass if needed
	      point.airmass = extn.airmass;
	      if (airmassIsFree && point.airmass < 1.) {
		cerr << "Prior " << name 
		     << " uses airmass but airmass<1 at exposure " 
		     << point.exposureName
		     << " device " << point.deviceName
		     << endl;
		exit(1);
	      }

	      // Map the coordinates to get [xy]Exposure
	      Assert(extn.startWcs);
	      point.args.xDevice = matchingXPix[iMatch];
	      point.args.yDevice = matchingYPix[iMatch];
	      extn.startWcs->toWorld(point.args.xDevice, point.args.yDevice,
				     point.args.xExposure, point.args.yExposure);
	      
	      Assert(extn.map);
	      point.map = extn.map;
	     
	      point.args.color = 0.;
	      // Initialize magOut using the map
	      point.magOut = point.map->forward(point.magIn, point.args);
	      refs.push_back(point);

	      if (colorIsFree) {
		// Need a second reference point with different color:
		point.args.color = 1.;
		point.magOut = point.map->forward(point.magIn, point.args);
		refs.push_back(point);
	      }

	      // We've found a relevant catalog.
	      foundOne = true;
	    } // End section for relevant extension found.

	  } // End loop over exposures matching this input line

	} // End loop over all extensions

	// Warn if exposure not found
	if (!foundOne) {
	  cerr << "WARNING: No fitted catalogs found for prior requesting exposure/device "
	       << exposureRegex << " / " << deviceName
	       << endl;
	}
      } // End processing of an input line specifying an exposure in the prior

    } // Done reading this prior.
    out.push_back(new PhotoPrior(refs,
				 sigma, name,
				 zeropoint, airmassCoeff, colorCoeff,
				 zpIsFree, airmassIsFree, colorIsFree));
  } // No more input from the file.

  return out;
}

