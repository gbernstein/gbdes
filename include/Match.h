// Match objects from multiple catalogs together
// And then adjust parameters of their PixelMaps to minimize
// differences from mean world coordinates.

#ifndef MATCH_H
#define MATCH_H

#include <list>
#include <set> 
using std::list;
#include <string>
#include "Std.h"
#include "LinearAlgebra.h"
#include "Bounds.h"
#include "Astrometry.h"
#include "PixelMap.h"
#include "PixelMapCollection.h"
#include "SymmetricUpdater.h"
#include "Units.h"

namespace astrometry {

  using linalg::SymmetricUpdater;
  class Match;  // Forward declaration

  // Some typedefs for proper motion
  typedef linalg::SVector<double, 5> PMSolution;
  typedef linalg::SMatrix<double,5,5> PMCovariance;
  // See Units.h for conventions of units.

  // Projection matrix from 5d PM to 2d position:
  typedef linalg::SMatrix<double,2,5> PMProjector;  

  class Detection {
  public:
    EIGEN_NEW
    long catalogNumber;
    long objectNumber;
    double xpix;
    double ypix;
    double color;
    double xw;
    double yw;
    // For proper motion and parallax:
    PMProjector* pmProj;
    const PMProjector& getProjector() const {return *pmProj;}
    
    // Inverse covariance, in world coords, for fitting.
    // This will include systematic errors and have been
    // increased by the weight factor:
    Matrix22 invCovFit;
    // This is the initial pure measurement error (before any sys errors or weights)
    Matrix22 invCovMeas;
    // This is the weighting factor applied here.  We will un-apply it when
    // calculating sigma of deviations for sigma clipping.
    double fitWeight;  
    // This is the expectation of (dx * invCovFit * dx), taking
    // into account the weight of this measurement in xmean for dx=x-xmean
    // Tends to 2 for many observations (5 for PMDetection)
    double expectedChisq;
    bool isClipped;
    Match* itsMatch;
    const SubMap* map;

    Detection(): color(astrometry::NODATA), pmProj(nullptr),expectedChisq(0.),
		 isClipped(false), itsMatch(nullptr), map(nullptr)  {}
    virtual ~Detection() {if (pmProj) delete pmProj;}

    // Build the projection matrix for this detection.
    void buildProjector(double pmTDB,	// Time in years from PM reference epoch
			const Vector3& xObs, // Observatory position, barycentric ICRS, in AU
			SphericalCoords* fieldProjection); // Projection used for world coords

    // Calculate chisq relative to best prediction of itsMatch,
    // using full fitting covariance but *before* any weighting factor is applied.
    // This and the resid*() routines do *not* automate remapping to world coords
    // for clipped detections.
    virtual double trueChisq() const;
    // Report 2d residual of this Detection to itsMatch prediction in world coordinates
    Vector2 residWorld() const;
    // Or the residual mapped back to pixel coordinates
    Vector2 residPix() const;

    // Get a rough estimate of a 1-dimensional world-coordinate error
    double getSigma() const {
      double sigsq = 0.5*(invCovMeas(0,0)+invCovMeas(1,1));
      return sigsq > 0. ? sqrt(fitWeight/sigsq) : 0.;
    }
  };

  class PMDetection: public Detection {
    // Derive a class which is for Gaia or other measurements that come with full 5d
    // constraints
  public:
    EIGEN_NEW
    PMSolution pmMean;
    PMCovariance invCovPM;
    // This function sets up the values and covariances to refer to
    // a reference PM epoch shifted by the specified number of years.
    // Also will fill in xw,yw, and covariances in the base class
    // members to refer to the (possibly shifted) pmTDB.  Then the
    // class can be sliced to a Detection and work for fitting without
    // free PM/parallax.
    void shiftReferenceEpoch(double tdbShift);

    // Calculate full *5d* chisq relative to best prediction of itsMatch,
    // using full fitting covariance but *before* any weighting factor is applied.
    virtual double trueChisq() const;
  };
    
  class Match {
    // Class for a set of matched Detections, with no PM or parallax freedom
  protected:
    list<Detection*> elist;
    mutable int nFit;	// Number of un-clipped points with non-zero weight in fit
    bool isReserved;	// Do not contribute to re-fitting if true

    // These flags will let a Match keep track of its internal state and
    // avoid unnecessary recalculations.  Internal states can change
    // even for logically const instance.
    //
    // First flag is set if the Fisher matrices, dof, and expected chisq's are
    // current, and the solution matrices for PM's.  These only change if
    // a Detection is added or removed from the set being fit.
    mutable bool isPrepared;  
    // These flags are set if the world coordinates of detection are valid.
    // They need to be reset whenever the coordinate maps' parameters change
    mutable bool isMappedFit;  // Detections involved in Fit are current
    mutable bool isMappedAll;  // All detections (incl. clipped) are current
    // This flag is set if the best-fit parameters of the Match position/PM
    // are current.  isPrepared && isMappedFit are prereqs for being solved.
    mutable bool isSolved;

    // Here are cached quantities for the overall solution:
    mutable Vector2 xyMean;   // Best-fit centroid (valid if isSolved).
    mutable Matrix22 centroidF;  // Fisher matrix for centroid (=inverse cov)
    mutable int dof;

    // Call this to redo the Fisher matrix if needed (=>set isPrepare)
    virtual void prepare() const;
    // Other state-changing calls are public.
    
  public:
    Match(Detection* e);

    // Add and remove automatically update itsMatch of the Detection
    void add(Detection* e);
    void remove(Detection* e);
    // Remove a Detection from the match given an iterator to it,
    // optionally deleting the Detection:
    list<Detection*>::iterator erase(list<Detection*>::iterator i,
				     bool deleteDetection=false);
    // Mark all members of match as unmatched, or optionally delete them,
    // then empty the list:
    void clear(bool deleteDetections=false);
    int size() const {return elist.size();}
    // True if a Detection will contribute to chisq:
    static bool isFit(const Detection* e);
    // Number of points that would have nonzero weight in next fit
    int fitSize() const {return nFit;}
    // Recount the number of objects contributing to fit (e.g. after
    // meddling with object weights)
    void countFit();
    int getDOF() const {prepare(); return dof;}

    // Is this object to be reserved from re-fitting?
    bool getReserved() const {return isReserved;}
    void setReserved(bool b) {isReserved = b;}

    // State-changing routines
    void mapsHaveChanged() const {
      // World coordinates of Detections no longer valid
      isMappedFit = false;
      isMappedAll = false;
      isSolved = false;
    }
    void forceRecalculation() const {
      // Mark everything as invalid
      mapsHaveChanged();
      isPrepared = false;
    }

    // Remap Detections to world coords with current map, either all or just
    // the ones involved in fitting
    virtual void remap(bool doAll = true) const;  
    // Recalculate best-fit centroid (or PM)
    virtual void solve() const;     

    // Get predicted position (and inverse covariance) for a Detection.
    // The argument is needed only if there is full PM solution.
    virtual Vector2 predict(const Detection* d = nullptr) const {
      solve();
      return xyMean;
    }
    virtual Matrix22 predictFisher(const Detection* d = nullptr)  const {
      prepare();
      return centroidF;
    }

    // Increment chisq, beta, and alpha for this match.
    // Returned integer is the DOF count.  This *does* remap points
    // but only those used for parameter fitting.
    // reuseAlpha=true will skip the incrementing of alpha.
    virtual int accumulateChisq(double& chisq,
				DVector& beta,
				SymmetricUpdater& updater,
				bool reuseAlpha=false);
   
    // sigmaClip returns true if clipped, 
    // and deletes the clipped guy if 2nd arg is true.  
    virtual bool sigmaClip(double sigThresh,
			   bool deleteDetection=false); 
    void clipAll(); // Mark all detections as clipped

    // Chisq for this match, and largest-sigma-squared deviation
    // 2 arguments are updated with info from this match.
    virtual double chisq(int& dofAccum, double& maxDeviateSq) const;

    typedef list<Detection*>::iterator iterator;
    typedef list<Detection*>::const_iterator const_iterator;
    iterator begin() {return elist.begin();}
    iterator end() {return elist.end();}
    const_iterator begin() const {return elist.begin();}
    const_iterator end() const {return elist.end();}
  };

  class PMMatch: public Match {
    // Class for a set of matched Detections, with free PM and parallax
  public:
    PMMatch(Detection* e);

    void setParallaxPrior(double p) {parallaxPrior = p;}

    virtual void solve() const;  // Recalculate best-fit PM 

    // Increment chisq, beta, and alpha for this match.
    // Returned integer is the DOF count.  
    // reuseAlpha=true will skip the incrementing of alpha.
    virtual int accumulateChisq(double& chisq,
				DVector& beta,
				SymmetricUpdater& updater,
				bool reuseAlpha=false);
   
    // sigmaClip returns true if clipped, 
    // and deletes the clipped guy if 2nd arg is true.  
    virtual bool sigmaClip(double sigThresh,
			   bool deleteDetection=false); 

    // Chisq for this match, and largest-sigma-squared deviation
    // 2 arguments are updated with info from this match.
    virtual double chisq(int& dofAccum, double& maxDeviateSq) const;

    // 
    const PMSolution& getPM() const {return pm;}
    const PMCovariance& getInvCovPM() const {prepare(); return pmFisher;}
    
    // Get predicted position (and inverse covariance) for a Detection.
    // The argument is needed only if there is full PM solution.
    virtual Vector2 predict(const Detection* d = nullptr) const;
    virtual Matrix22 predictFisher(const Detection* d = nullptr)  const;

  protected:
  private:
    double parallaxPrior;   // sigma for prior on parallax (radians)

    mutable PMSolution pm;          // Chisq-minimizing position/parallax/pm

    virtual void prepare() const;  // Recalculate and cache, if needed:
    mutable PMCovariance pmFisher;  // Fisher matrix for PM
    mutable PMCovariance pmCov;     // Inverse of Fisher matrix
    mutable double priorChisq;      // Chisq term quadratic in PMDetection means
    mutable PMSolution priorMean;   // PM contribution from PMDetection means
  };

  typedef list<Match*> MCat;

  // Class that aligns all coordinates
  class CoordAlign {
  private:
    list<Match*>& mlist;
    PixelMapCollection& pmc;
    double relativeTolerance;
    set<int> frozenParameters;  // Keep track of degenerate parameters
    map<string, set<int>> frozenMaps; // Which atoms have which params frozen
  public:
    CoordAlign(PixelMapCollection& pmc_,
	       list<Match*>& mlist_): mlist(mlist_),
				      pmc(pmc_), 
				      relativeTolerance(0.001)  {}

    // Re-map Detections using current params - either all or just those being fit.
    void remap(bool doAll=true); 
    // Fitting routine: returns chisq of current fit, updates params. and remaps.
    double fitOnce(bool reportToCerr=true,
		   bool inPlace=false);	  // Set inPlace to save space, but can't debug singularities
    // Conduct one round of sigma-clipping.  If doReserved=true, 
    // then only clip reserved Matches.  If =false, then
    // only clip non-reserved Matches.
    int sigmaClip(double sigThresh, bool doReserved=false,
		  bool clipEntireMatch=false);
    // Calculate total chisq.  doReserved same meaning as above.
    double chisqDOF(int& dof, double& maxDeviate, bool doReserved=false) const;
    void setParams(const DVector& p);
    DVector getParams() const {return pmc.getParams();}
    int nParams() const {return pmc.nParams();}
    // This is the calculation of normal equation components used in fitting.
    // reuseAlpha=true will leave alpha unchanged, way faster.
    void operator()(const DVector& params, double& chisq,
		    DVector& beta, DMatrix& alpha,
		    bool reuseAlpha=false);
    void setRelTolerance(double tol) {relativeTolerance=tol;}
    // Return count of useful (un-clipped) Matches & Detections.
    // Count either reserved or non-reserved objects, and require minMatches useful
    // Detections for a valid match:
    void count(long int& mcount, long int& dcount, bool doReserved=false, int minMatches=2) const;
    void count(long int& mcount, long int& dcount, bool doReserved, int minMatches, long catalog) const;
  };

  // Function that will execute a Marquardt-Levenberg fit to optimize the
  // parameters of some model to a simple list of detections.
  // Assigns equal uncertainty to all test points, given by sigma
  void mapFit(list<Detection*> testPoints, PixelMap* pm, double sigma=1.);

} // namespace astrometry

#endif
