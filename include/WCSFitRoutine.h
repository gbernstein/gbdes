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

class WCSFit {
public:
    WCSFit(int minMatches = 2, int verbose = 0);
    WCSFit(Fields &fields_, std::vector<shared_ptr<Instrument>> instruments_, ExposuresHelper exposures_,
           std::vector<int> extensionExposureNumbers, std::vector<int> extensionDevices,
           astrometry::YAMLCollector inputYAML, std::vector<shared_ptr<astrometry::Wcs>> wcss,
           const std::vector<int> &sequence, const std::vector<LONGLONG> &extns, const std::vector<LONGLONG> &objects,
           const std::vector<int> &exposureColorPriorities = std::vector<int>(), double sysErr = 2.0, double refSysErr = 2.0,
           int minMatches = 2, std::string skipObjectsFile = "", std::string fixMaps = "", bool usePM = true,
           double pmPrior=100.0, double parallaxPrior=10.0, int verbose = 0);

    int minMatches;
    int verbose;

    // Field information
    Fields fields;

    // Instrument and device (i.e. detector/ccd) information
    std::vector<unique_ptr<Instrument>> instruments;

    // The table of exposures
    std::vector<unique_ptr<Exposure>> exposures;

    // Extension tables:
    std::vector<unique_ptr<ColorExtension>> colorExtensions;
    std::vector<unique_ptr<Extension>> extensions;

    astrometry::PixelMapCollection mapCollection;

    // List of all Matches - they will hold pointers to all Detections too.
    astrometry::MCat matches;

    void setExposures(std::vector<unique_ptr<Exposure>> expos, double sysErr, double refSysErr);

    void setRefWCSNames();

    void addMap(astrometry::YAMLCollector &inputYAML, std::string mapName,
                std::vector<std::string> mapParams);

    void setupMaps(astrometry::YAMLCollector &inputYAML, std::string fixMaps = "");

    void setMatches(const std::vector<int> &sequence, const std::vector<LONGLONG> &extensions,
                    const std::vector<LONGLONG> &objects, ExtensionObjectSet skipSet, bool usePM = true);

    void setObjects(int i, const std::map<std::string, std::vector<double>> &tableMap,
                    const std::string &xKey, const std::string &yKey,
                    const std::vector<std::string> &xyErrKeys, const std::string &idKey = "",
                    const std::string &pmCovKey = "", const std::string &magKey = "",
                    const int &magKeyElement = 0, const std::string &magErrKey = "",
                    const int &magErrKeyElement = 0, const std::string &pmRaKey = "",
                    const std::string &pmDecKey = "", const std::string &parallaxKey = "",
                    const std::vector<std::vector<double>> &fullCov = std::vector<std::vector<double>>(0.0));

    void reprojectWCSs();

    void fit(double maxError = 100., int minFitExposures = 200, double reserveFraction = 0.2,
             int randomNumberSeed = 1234, double minimumImprovement = 0.02, double clipThresh = 5.0,
             double chisqTolerance = 0.001, bool clipEntireMatch = false, bool divideInPlace = false,
             bool purgeOutput = false, double minColor = -10.0, double maxColor = 10.0);

    std::map<std::string, vector<float>> getOutputCatalog();
    std::map<std::string, vector<float>> getPMCatalog(vector<vector<float>> &PMMean,
                                                      vector<vector<float>> &PMInvCov);
    std::map<std::string, vector<float>> getStarCatalog(vector<vector<float>> &starInvCov);

    void saveResults(std::string outWcs, std::string outCatalog, std::string starCatalog);
};

#endif
