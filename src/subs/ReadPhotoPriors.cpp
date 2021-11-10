// Code to read PhotoPriors from files.
// Used in PhotoFit.cpp.
// Now reading the priors in YAML format.

/***
Recall that the prior constraining residuals to the model
(m_out - m_in) = zpt + a*(airmass-1) + b*color + c*apcorr
where the LHS is the shift being applied to cataloged mags
for some location on some device in some exposure.

In the prior file, YAML root node will be a dict with the
keys being names of prior constraint on a group of exposures, and the values
being dictionaries with these keys:

exposures:  one or a list of regexes specifying which exposures to constrain
sigma:      the uncertainty assigned to (mout-min) in forming chisq (default: 0.001)
device:     which device on the exposure to use (default: "S2")
xpixel:     x pixel location to use (default or negative: use ctr of device)
ypixel:     y pixel location to use (default or negative: use ctr of device)

zeropoint:  initial value for zpt coefficient  
airmass:    initial value for airmass coefficient
color:      initial value for color coefficient
apcorr:     initial value for apcorr coefficient
 (all default to zero)

free:       list of which of (zeropoint,airmass,color,apcorr)
coefficients should be allowed to vary during fit.  (default: none)

So the default will be an "absolute" constraint forcing the mag shift to
zero for the specified exposures.
***/
#include <fstream>
#include <sstream>

#include "Std.h"
#include "StringStuff.h"
#include "PixelMapCollection.h"
#include "PhotoMatch.h"
#include "PhotoMapCollection.h"
#include "FitSubroutines.h"
#include "yaml-cpp/yaml.h"

using namespace std;
using namespace stringstuff;
using namespace photometry;

void formatError(string filename, string line) {
  cerr << "Format error reading prior file " << filename << " on input line:" << endl;
  cerr << line << endl;
  exit(1);
}

list<PhotoPrior*>
readPriors(string filename, 
	   const vector<unique_ptr<Instrument>>& instruments, 
	   const vector<unique_ptr<Exposure>>& exposures, 
	   const vector<Photo::Extension*>& extensions ) {

  ifstream ifs(filename.c_str());
  if (!ifs) {
    cerr << "Could not open photometric priors file " << filename << endl;
    exit(1);
  }

  list<PhotoPrior*> out;

  try {
    YAML::Node root = YAML::Load(ifs);
    if (!root.IsMap()) {
      // A valid YAML document that does not have the structure 
      // that our files should have.
      throw PhotometryError("Invalid structure for photo prior file " + filename);
    }
    
    // Every entry in the root node map is a prior specification
    for (YAML::const_iterator i = root.begin();
	 i != root.end();
	 ++i) {
      string priorName = i->first.as<string>();
      const YAML::Node& node=i->second;
      if (!node.IsMap() || !node["Exposures"]) 
	throw PhotometryError("Non-map or missing exposures key in prior YAML entry " +
			      priorName);

      // Read the exposures
      list<string> exposureRegexes;
      if (node["Exposures"].IsSequence()) {
	exposureRegexes = node["Exposures"].as<list<string>>();
      } else {
	exposureRegexes.push_back(node["Exposures"].as<string>());
      }

      // Width of magnitude prior:
      double sigma = 0.001;
      if (node["Sigma"])
	sigma = node["Sigma"].as<double>();

      // Location of reference point:
      string deviceName = "S2";
      double xpixel = -1.;
      double ypixel = -1.;
      if (node["Device"])
	deviceName = node["Device"].as<string>();
      if (node["Xpixel"])
	xpixel = node["Xpixel"].as<double>();
      if (node["Ypixel"])
	ypixel = node["Ypixel"].as<double>();

      // Coefficients
      double zeropoint=0.;
      double airmassCoeff=0.;
      double colorCoeff=0.;
      double apcorrCoeff=0.;
      if (node["Zeropoint"]) 
	zeropoint = node["Zeropoint"].as<double>();
      if (node["Airmass"])
	airmassCoeff = node["Airmass"].as<double>();
      if (node["Color"])
	colorCoeff = node["Color"].as<double>();
      if (node["Apcorr"])
	apcorrCoeff = node["Apcorr"].as<double>();


      // Fixed or free parameters?
      bool zpIsFree = false;
      bool airmassIsFree = false;
      bool colorIsFree = false;
      bool apcorrIsFree = false;
      list<string> freeList;
      if (node["Free"]) {
	if (node["Free"].IsSequence()) {
	  freeList = node["Free"].as<list<string>>();
	} else {
	  freeList.push_back(node["Free"].as<string>());
	}
	for (auto s : freeList) {
	  if (nocaseEqual(s, "zeropoint"))
	    zpIsFree = true;
	  else if (nocaseEqual(s, "airmass"))
	    airmassIsFree = true;
	  else if (nocaseEqual(s, "color"))
	    colorIsFree = true;
	  else if (nocaseEqual(s, "apcorr"))
	    apcorrIsFree = true;
	  else 
	    throw PhotometryError("Unknown PhotoPrior coefficient <" + s + ">");
	}
      }

      // The input magnitude for each reference point:
      double instrumentalMag = 0.;

      // Now construct reference points from matching exposures
      list<PhotoPriorReferencePoint> refs;

      // Now look for all exposures that match the request here.
      vector<int> matchingExposures;
      vector<int> matchingDevices;
      vector<double> matchingXPix;
      vector<double> matchingYPix;
      for (int iExpo = 0; iExpo < exposures.size(); iExpo++) {
	if (!exposures[iExpo]) continue;
	if (regexMatchAny(exposureRegexes, exposures[iExpo]->name)) {
	  // For each matching exposure: look for instrument, device name, 
	  // and get center coords if needed
	  int iInst = exposures[iExpo]->instrument;
	  if (iInst<0 || !instruments[iInst])
	    continue;

	  Instrument& inst = *instruments[iInst];
	  int iDevice = inst.deviceNames.indexOf(deviceName);
	  if (iDevice < 0) {
	    // This instrument has no such device.
	    continue;
	  }
	  // Save away the device's number and the pixel coordinates we will use.
	  matchingExposures.push_back(iExpo);
	  matchingDevices.push_back(iDevice);
	  Position<double> center = inst.domains[iDevice].center();
	  matchingXPix.push_back( xpixel<0. ? center.x : xpixel);
	  matchingYPix.push_back( ypixel<0. ? center.y : ypixel);
	}
      } // end of loop finding matching exposures
	    
      // Troll extensions looking for exposure + device match
      for (int iExtn = 0; iExtn < extensions.size(); iExtn++) {
	if (!extensions[iExtn]) continue;
	Photo::Extension& extn = *extensions[iExtn];
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
	    if ( (airmassIsFree || airmassCoeff!=0.) && point.airmass < 1.) {
	      cerr << "Prior " << priorName 
		   << " uses airmass but airmass<1 at exposure " 
		   << point.exposureName
		   << " device " << point.deviceName
		   << endl;
	      exit(1);
	    }
	    point.apcorr = extn.apcorr;

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

	    // Need a second reference point with different color:
	    // to force any color dependence to agree.
	    point.args.color = 1.;
	    point.magOut = point.map->forward(point.magIn, point.args);
	    refs.push_back(point);

	  } // End section for relevant extension found.

	} // End loop through exposure/device pairs
      } // End loop over extensions
      
      // Warn if no exposures found
      if (refs.empty()) {
	cerr << "WARNING: No fitted catalogs found for prior "
	     << priorName
	     << endl;
      } else {
	out.push_back(new PhotoPrior(refs,
				     sigma, priorName,
				     zeropoint, airmassCoeff, colorCoeff, apcorrCoeff,
				     zpIsFree, airmassIsFree, colorIsFree, apcorrIsFree));
      }
    } // End loop of YAML prior nodes

  } catch (YAML::Exception& e) {
    cerr << "YAML error during ReadPhotoPrior" << endl;
    throw;
  }
  return out;
}

