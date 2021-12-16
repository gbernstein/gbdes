#ifndef WCSFIT_FIT_H
#define WCSFIT_FIT_H

#include <map>
#include <algorithm>

#include "Std.h"
#include "Astrometry.h"
#include "FitsTable.h"
#include "StringStuff.h"
#include "Pset.h"
#include "PixelMapCollection.h"
#include "YAMLCollector.h"
#include "Match.h"
#include "Instrument.h"

#include "FitSubroutines.h"
#include "WcsSubs.h"
#include "MapDegeneracies.h"

#ifdef _OPENMP
#include <omp.h>
#endif

using namespace std;
using namespace astrometry;
using namespace stringstuff;
using img::FTable;

class WCSFit {
  public:
    WCSFit(int minMatches=2,
           int verbose=0);
    WCSFit(Fields & fields_,
           vector<shared_ptr<Instrument>> instruments_,
           ExposuresHelper exposures_,
           vector<int> extensionExposureNumbers,
           vector<int> extensionDevices,
           YAMLCollector inputYAML,
           vector<shared_ptr<astrometry::Wcs>> wcss,
           vector<int> sequence,
           vector<LONGLONG> extns,
           vector<LONGLONG> objects,
           double sysErr=2.0,
           double refSysErr=2.0,
           int minMatches=2,
           string skipObjectsFile="",
           string fixMaps="",
           bool usePM=true,
           int verbose=0
           );
    
    int minMatches;
    int verbose;

    // Field information
    Fields fields;
    
    // Instrument and device (i.e. detector/ccd) information
    vector<unique_ptr<Instrument>> instruments;
    
    // The table of exposures
    vector<unique_ptr<Exposure>> exposures;

    // Extension tables:
    vector<unique_ptr<ColorExtension>> colorExtensions;
    vector<unique_ptr<Extension>> extensions;

    PixelMapCollection mapCollection;

    // List of all Matches - they will hold pointers to all Detections too.
    MCat matches;

    void setExposures(vector<unique_ptr<Exposure>> expos, double sysErr, double refSysErr);

    void setRefWCSNames();

    void addMap(YAMLCollector& inputYAML, string mapName, vector<string> mapParams);

    void setupMaps(YAMLCollector& inputYAML, string fixMaps="");
    
    void setMatches(vector<int> sequence, vector<LONGLONG> extensions, vector<LONGLONG> objects,
                    ExtensionObjectSet skipSet, bool usePM=true);

    void setObjects(int i, map<string, vector<double>> tableMap,
                    string xKey, string yKey, 
                    vector<string> xyErrKeys, string idKey="", string pmCovKey="", string magKey="",
                    int magKeyElement=0, string magErrKey="",
                    int magErrKeyElement=0, string pmRaKey="", string pmDecKey="", string parallaxKey="");

    void reprojectWCSs();

    void fit(double maxError=100., int minFitExposures=200, double reserveFraction=0.2,
             int randomNumberSeed=1234, double minimumImprovement=0.02, double clipThresh=5.0,
             double chisqTolerance=0.001, bool clipEntireMatch=false, bool divideInPlace=false,
             bool purgeOutput=false, double minColor=-10.0, double maxColor=10.0);

    void saveResults(string outWcs, string outCatalog, string starCatalog);
      
};

#endif
