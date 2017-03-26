// Match objects from multiple catalogs together
// And then adjust parameters of their PhotoMaps to minimize
// differences from mean world coordinates.

#ifndef PHOTOMATCH_H
#define PHOTOMATCH_H

#include <list>
#include <set> 
using std::list;
#include <string>
#include "Std.h"
#include "UseTMV.h"
#include "TMV_SymBand.h"
#include "Bounds.h"
#include "PhotoMapCollection.h"
#include "AlphaUpdater.h"

#ifdef _OPENMP
#include <omp.h>
#endif

namespace photometry {

  class Match;  // Forward declaration

  class Detection {
  public:
    long catalogNumber;
    long objectNumber;
    double magIn;
    double sigma;	// Uncertainty, in mag
    photometry::PhotoArguments args;
    double magOut;
    // weight, nominally inverse sigma squared of output mag
    double wt;
    // Inverse square of the sigma used for sig-clipping:
    double clipsq;
    bool isClipped;
    const Match* itsMatch;
    const SubMap* map;
    Detection(): itsMatch(0), map(0), isClipped(false) {}
  };
  
  class Match {
  private:
    list<Detection*> elist;
    int nFit;	// Number of un-clipped points with non-zero weight in fit
    bool isReserved;	// Do not contribute to re-fitting if true
    static bool isFit(const Detection* e);
  public:
    Match(Detection* e);
    // Add and remove automatically update itsMatch of the Detection
    void add(Detection* e);
    void remove(Detection* e);
    // Mark all members of match as unmatched, or optionally delete them,
    // then empty the list:
    list<Detection*>::iterator erase(list<Detection*>::iterator i,
				     bool deleteDetection=false);
    // Mark all members of match as unmatched, or optionally delete them,
    // then empty the list:
    void clear(bool deleteDetections=false);

    int size() const {return elist.size();}
    // Number of points that would have nonzero weight in next fit
    int fitSize() const {return nFit;}
    // Update and return count of above:
    void countFit(); 

    void remap();  // Remap each point, i.e. make new magOut
    // Mean of un-clipped output mags, optionally with total weight - no remapping done
    void getMean(double& mag) const;
    void getMean(double& mag, double& wt) const;

    // Is this object to be reserved from re-fitting?
    bool getReserved() const {return isReserved;}
    void setReserved(bool b) {isReserved = b;}

    // Returned integer is the DOF count
    int accumulateChisq(double& chisq,
			DVector& beta,
			astrometry::AlphaUpdater& updater);
    // sigmaClip returns true if clipped, and deletes the clipped guy
    // if 2nd arg is true.
    bool sigmaClip(double sigThresh,
		   bool deleteDetection=false); 
    void clipAll(); // Mark all detections as clipped

    // Chisq for this match, and largest-sigma-squared deviation
    // 2 arguments are updated with info from this match.
    double chisq(int& dof, double& maxDeviateSq) const;

    typedef list<Detection*>::iterator iterator;
    typedef list<Detection*>::const_iterator const_iterator;
    iterator begin() {return elist.begin();}
    iterator end() {return elist.end();}
    const_iterator begin() const {return elist.begin();}
    const_iterator end() const {return elist.end();}
  };

  // A prior on zeropoints of exposures will be forcing reference points into agreement.
  // This class is for the reference points
  class PhotoPriorReferencePoint {
  public:
    PhotoPriorReferencePoint(): airmass(1.), map(0), isClipped(false) {}
    string exposureName;
    string deviceName;
    double magIn;  // Magnitude assigned to 1 count per second.
    double magOut;	// Mag after mapping.
    PhotoArguments args;
    double airmass;
    const SubMap* map;	// Photometric map applying to this reference point
    bool isClipped;
  };

  // Class that manifests some prior constraint on the agreement of zeropoints of different
  // exposures; typically one PhotoPrior for a photometric night's data

  class PhotoPrior {
  public:
    PhotoPrior(list<PhotoPriorReferencePoint> points_,
	       double sigma_,
	       string name_="",
	       double zeropoint = 0.,
	       double airmassCoefficient = 0.,
	       double colorCoefficient = 0.,
	       bool freeZeropoint = false,
	       bool freeAirmass = false,
	       bool freeColor = false);
    int nParams() const {return nFree;}
    DVector getParams() const;
    void setParams(const DVector& p);
    // If this call is true, prior cannot be used for fitting:
    bool isDegenerate() const {return nFit < nFree;}

    string getName() const {return name;}
    double getSigma() const {return sigma;}
    double getZeropoint() const {return m;}
    double getAirmass() const {return a;}
    double getColor() const {return b;}
    bool zeropointIsFree() const {return mIsFree;}
    bool airmassIsFree() const {return aIsFree;}
    bool colorIsFree() const {return bIsFree;}

    // Locations of parameters in global vector
    int startIndex() const {return globalStartIndex;}
    int mapNumber() const {return globalMapNumber;}

    // Recalculate *all* the reference points with current parameters
    void remap();
    // Recalculate the fittable points and increment chisq and fitting vector/matrix
    int accumulateChisq(double& chisq,
			 DVector& beta,
			//**tmv::SymMatrix<double>& alpha);
			astrometry::AlphaUpdater& updater);
    // sigmaClip returns true if clipped one, and only will clip worst one - no recalculation
    bool sigmaClip(double sigThresh);
    // Mark all of this prior as clipped (will no longer be fit)
    void clipAll();
    // Chisq for this match, also updates dof for free parameters - no recalculation
    double chisq(int& dof) const;
    
    // For a printed report to the stated stream:
    static void reportHeader(ostream& os);
    void report(ostream& os) const;

  private:
    double sigma;	// strength of prior (mags) at each reference point
    list<PhotoPriorReferencePoint> points;
    string name;
    int nFree;	// Count of number of free parameters among m, a, b.
    friend class PhotoAlign;	// PhotoAlign can change the global indices
    int globalStartIndex;	// Index of m/a/b in PhotoMapCollection param vector
    int globalMapNumber;	// map number for resource locking
    bool mIsFree;
    bool aIsFree;
    bool bIsFree;
    double m;	// Mag zeropoint
    double a;	// X-1 airmass coefficient
    double b;	// color coefficient

    int nFit;	// Number of reference points being fit
    void countFit();	// Update the above number
  };

  // Class that fits to make magnitudes agree
  class PhotoAlign {
  private:
    list<Match*>& mlist;
    PhotoMapCollection& pmc;
    list<PhotoPrior*>& priors;
    double relativeTolerance;
    set<int> frozenParameters;  // Keep track of degenerate parameters
    map<string, set<int>> frozenMaps; // Which atoms have which params frozen
    int nPriorParams;
    int maxMapNumber;
    void countPriorParams();  // Update parameter counts, indices, map numbers for priors
  public:
    PhotoAlign(PhotoMapCollection& pmc_,
	       list<Match*>& mlist_,
	       list<PhotoPrior*>& priors_): mlist(mlist_),
					    pmc(pmc_), 
					    priors(priors_),
					    relativeTolerance(0.001)  {countPriorParams();}

    // Conduct one round of sigma-clipping.  If doReserved=true, 
    // then only clip reserved Matches.  If =false, then
    // only clip non-reserved Matches.
    int sigmaClip(double sigThresh, bool doReserved=false,
		  bool clipEntireMatch=false);
    // Clip the priors to eliminate non-photometric exposures.  Set the flag to
    // eliminate all use of a prior that has any outliers (one bad exposure means full night
    // is assumed non-photometric).  Returns number of priors with a clip.
    int sigmaClipPrior(double sigThresh, bool clipEntirePrior=false);

    // Calculate total chisq.  doReserved same meaning as above.
    void setParams(const DVector& p);
    DVector getParams() const;
    int nParams() const {return pmc.nParams() + nPriorParams;}

    void remap();	// Re-map all Detections and Priors using current params
    double fitOnce(bool reportToCerr=true);	// Returns chisq of previous fit, updates params.
    double chisqDOF(int& dof, double& maxDeviate, bool doReserved=false) const;
    // This is the () operator that is called by Marquardt to try a fit with new params
    void operator()(const DVector& params, double& chisq,
		    DVector& beta, tmv::SymMatrix<double>& alpha);

    void setRelTolerance(double tol) {relativeTolerance=tol;}
    // Return count of useful (un-clipped) Matches & Detections.
    // Count either reserved or non-reserved objects, and require minMatches useful
    // Detections for a valid match:
    void count(long int& mcount, long int& dcount, 
	       bool doReserved=false, int minMatches=2) const;
    // This gives the count for just a single selected input catalog:
    void count(long int& mcount, long int& dcount, 
	       bool doReserved, int minMatches, long catalog) const;
  };

} // namespace astrometry

#endif
