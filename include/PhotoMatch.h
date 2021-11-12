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
#include "LinearAlgebra.h"
#include "Bounds.h"
#include "PhotoMapCollection.h"
#include "SymmetricUpdater.h"

#ifdef _OPENMP
#include <omp.h>
#endif

namespace photometry {

  using linalg::SymmetricUpdater;
  class Match;  // Forward declaration

  class Detection {
  public:
    long catalogNumber;
    long objectNumber;
    double magIn;
    photometry::PhotoArguments args;
    double magOut;
    // Inverse covariance for fitting.
    // This will include systematic errors.
    double invVar;
    // This is the weighting factor applied here.  We will apply it to
    // invCov when solving for centroid
    double fitWeight;  
    // This is the expectation of (dx * invCov * dx), taking
    // into account the weight of this measurement in xmean for dx=x-xmean
    // Tends to 1 for many observations
    double expectedTrueChisq;
    bool isClipped;
    const Match* itsMatch;
    const SubMap* map;
    Detection(): isClipped(false),
		 fitWeight(1.), expectedTrueChisq(0.),
		 itsMatch(nullptr), map(nullptr) {}
    // Get total error on magnitude:
    double getSigma() const {return 1./sqrt(invVar);}
    // Get chisq using total error but no weighting factor
    double trueChisq() const;

    double residMag() const;  // Get mag error after fit
    void clip() {isClipped=true; fitWeight=0.;}
  };

  /*** Old quantities ???
    double sigma;	// Uncertainty, in mag
    // weight, nominally inverse sigma squared of output mag
    double wt;
    // Inverse square of the sigma used for sig-clipping:
    double clipsq;
  ***/
  class Match {
  private:
    list<unique_ptr<Detection>> elist;
    bool isReserved;	// Do not contribute to re-fitting if true

    /* State maintenance */
    // First flag is set if the Fisher "matrix", dof, and expected chisq's are
    // current, and the solution matrices for mean.  These only change if
    // a Detection is added or removed from the set being fit.
    mutable bool isPrepared;  
    // These flags are set if the magOuts of Detections are valid
    // They need to be reset whenever the maps' parameters change
    mutable bool isMappedFit;  // Detections involved in Fit are current
    mutable bool isMappedAll;  // All detections (incl. clipped) are current
    // This flag is set if the best-fit centroid is current.
    // isPrepared && isMappedFit are prereqs for being solved.
    mutable bool isSolved;
    // This flag is set if all weights are 1 or 0
    mutable bool trivialWeights;

    // Here are cached quantities for the overall solution:
    mutable int nFit;	// Number of un-clipped points with non-zero weight in fit
    mutable double meanMag;   // Best-fit magnitude (valid if isSolved).
    mutable double meanF;  // pseudo-inverse-var for centroid (2nd deriv of weighted chisq)
    // The expected variance of the mean (accounting for weighting):
    mutable double trueMeanVar; 
    mutable int dof;

    // Call this to redo the centroid variance if needed (=>set isPrepare)
    virtual void prepare() const;
    // Other state-changing calls are public.
    
  public:
    explicit Match(unique_ptr<Detection> e);
    // Add and remove automatically update itsMatch of the Detection
    void add(unique_ptr<Detection> e);
    void remove(const Detection & e);
    // Remove a Detection from the match given an iterator to it.
    list<unique_ptr<Detection>>::iterator erase(list<unique_ptr<Detection>>::iterator i);
    // Delete all members of match.
    void clear();

    // True if a Detection will contribute to chisq:
    static bool isFit(const Detection & e);

    // Is this object to be reserved from re-fitting?
    bool getReserved() const {return isReserved;}
    void setReserved(bool b) {isReserved = b;}

    // Number of points that would have nonzero weight in next fit
    int fitSize() const {prepare(); return nFit;}
    // Return DOF in the fit
    int getDOF() const {prepare(); return dof;}
    
    // State-changing routines
    void mapsHaveChanged() const {
      // magOut of Detections no longer valid
      isMappedFit = false;
      isMappedAll = false;
      isSolved = false;
    }
    void forceRecalculation() const {
      // Mark everything as invalid
      mapsHaveChanged();
      isPrepared = false;
    }

    // Remap Detections to magOut current map, either all or just
    // the ones involved in fitting
    virtual void remap(bool doAll = true) const;  
    // Recalculate best-fit centroid
    virtual void solve() const;     

    // Mean mag and its inverse variance under weighted fit
    double getMean() const {solve(); return meanMag;}
    double getMeanFisher() {prepare(); return meanF;}

    // Returned integer is the DOF count
    int accumulateChisq(double& chisq,
			DVector& beta,
			SymmetricUpdater& updater,
			bool reuseAlpha=false);
    // sigmaClip returns true if clipped, and deletes the clipped guy
    // if 2nd arg is true.
    bool sigmaClip(double sigThresh,
		   bool deleteDetection=false); 
    void clipAll(); // Mark all detections as clipped

    // Chisq for this match, and largest-sigma-squared deviation
    // 2 arguments are updated with info from this match.
    double chisq(int& dofAccum, double& maxDeviateSq) const;

    typedef list<unique_ptr<Detection>>::iterator iterator;
    typedef list<unique_ptr<Detection>>::const_iterator const_iterator;
    iterator begin() {return elist.begin();}
    iterator end() {return elist.end();}
    const_iterator begin() const {return elist.begin();}
    const_iterator end() const {return elist.end();}
    int size() const {return elist.size();}
  };

  // A prior on zeropoints of exposures will be forcing reference points into agreement.
  // This class is for the reference points
  class PhotoPriorReferencePoint {
  public:
    PhotoPriorReferencePoint(): airmass(1.), map(nullptr), isClipped(false) {}
    string exposureName;
    string deviceName;
    double magIn;  // Magnitude assigned to 1 count per second.
    mutable double magOut;	// Mag after mapping.
    PhotoArguments args;
    double airmass;
    double apcorr;
    const SubMap* map;	// Photometric map applying to this reference point
    bool isClipped;
    void clip() {isClipped = true;}
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
	       double apcorrCoefficient = 0.,
	       bool freeZeropoint = false,
	       bool freeAirmass = false,
	       bool freeColor = false,
	       bool freeApcorr = false);
    int nParams() const {return nFree;}
    DVector getParams() const;
    void setParams(const DVector& p);
    // If this call is true, prior cannot be used for fitting:
    bool isDegenerate() const {prepare(); return dof<0;}  // ??? get rid of this?
    int getDOF() const {prepare(); return dof;}
    

    string getName() const {return name;}
    double getSigma() const {return sigma;}
    double getZeropoint() const {return m;}
    double getAirmass() const {return a;}
    double getColor() const {return b;}
    double getApcorr() const {return c;}
    bool zeropointIsFree() const {return mIsFree;}
    bool airmassIsFree() const {return aIsFree;}
    bool colorIsFree() const {return bIsFree;}
    bool apcorrIsFree() const {return cIsFree;}

    // Locations of parameters in global vector
    int startIndex() const {return globalStartIndex;}
    int mapNumber() const {return globalMapNumber;}

    // Set this when underlying data change (parameters of exposure's maps)
    void mapsHaveChanged() const {isMapped = false;}

    // Recalculate the fittable points and increment chisq and fitting vector/matrix
    int accumulateChisq(double& chisq,
			DVector& beta,
			SymmetricUpdater& updater,
			bool reuseAlpha=false);
    // sigmaClip returns true if clipped one, and only will clip worst one - no recalculation
    bool sigmaClip(double sigThresh);
    // Mark all of this prior as clipped (will no longer be fit)
    void clipAll();
    // Chisq for this match, also updates dof for free parameters - no recalculation
    double chisq(int& dofAccum) const;
    
    // For a printed report to the stated stream:
    static void reportHeader(ostream& os);
    void report(ostream& os) const;

  private:
    double sigma;	// strength of prior (mags) at each reference point
    list<PhotoPriorReferencePoint> points;
    string name;
    friend class PhotoAlign;	// PhotoAlign can change the global indices
    int globalStartIndex;	// Index of m/a/b in PhotoMapCollection param vector
    int globalMapNumber;	// map number for resource locking
    double m;	// Mag zeropoint
    double a;	// X-1 airmass coefficient
    double b;	// color coefficient
    double c;	// aperture correction coefficient
    bool mIsFree;
    bool aIsFree;
    bool bIsFree;
    bool cIsFree;

    // State maintenance / cached quantities
    mutable bool isPrepared;  // Calculated DOF and such?
    mutable bool isMapped;    // dm's for reference points are current
    void prepare() const;           // Recalc DOF and such
    void remap() const;       // Recalc dm's
    mutable int nFit;	// Number of reference points being fit
    int nFree;	        // Count of number of free parameters among m, a, b.
    mutable int dof;
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

    // Re-map all Detections and Priors using current params, either all or just
    // those being fit.
    void remap(bool doAll=true) const;	// ??? Need this ???

    // Fit parameters. Returns chisq of previous fit, updates params.
    double fitOnce(bool reportToCerr=true,
		   bool inPlace=false);	 // Set inPlace to save space, less debug

    // Conduct one round of sigma-clipping.  If doReserved=true, 
    // then only clip reserved Matches.  If =false, then
    // only clip non-reserved Matches.
    int sigmaClip(double sigThresh, bool doReserved=false,
		  bool clipEntireMatch=false,
		  bool logging=true);

    // Clip the priors to eliminate non-photometric exposures.  Set the flag to
    // eliminate all use of a prior that has any outliers (one bad exposure means full night
    // is assumed non-photometric).  Returns number of priors with a clip.
    int sigmaClipPrior(double sigThresh, bool clipEntirePrior=false);

    // Calculate total chisq.  doReserved same meaning as above.
    double chisqDOF(int& dof, double& maxDeviate, bool doReserved=false) const;

    // Get/set parameters
    void setParams(const DVector& p);
    DVector getParams() const;
    int nParams() const {return pmc.nParams() + nPriorParams;}

    // This is the calculation of normal equation components used in Marquardt fitting.
    // reuseAlpha=true will leave alpha unchanged, way faster.
    void operator()(const DVector& params, double& chisq,
		    DVector& beta, DMatrix& alpha,
		    bool reuseAlpha=false);

    // Set tolerance for fit convergence
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
