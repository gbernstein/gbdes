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
  // See Units.h for conventions of units and order of PM elements.

  // Projection matrix from 5d PM to 2d position:
  typedef linalg::SMatrix<double,2,5> PMProjector;  

  template<class M>
  double traceATB(const M& A, const M& B) {
    // function to return Trace(A^T B)
#ifdef USE_EIGEN
    return (A.array()*B.array()).sum();
#elif  USE_TMV
    return tmv::ElemProd(A,B).sumElements();
#endif
  }
  

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
    unique_ptr<PMProjector> pmProj;
    const PMProjector& getProjector() const {return *pmProj;}
    
    // Inverse covariance, in world coords, for fitting.
    // This will include systematic errors.
    Matrix22 invCov;
    // This is the weighting factor applied here.  We will apply it to
    // invCov when solving for centroid / PM.
    double fitWeight;  
    // This is the expectation of (dx * invCov * dx), taking
    // into account the weight of this measurement in xmean for dx=x-xmean
    // Tends to 2 for many observations (5 for PMDetection)
    double expectedTrueChisq;
    bool isClipped;
    Match* itsMatch;
    const SubMap* map;

    Detection(): color(astrometry::NODATA), pmProj(nullptr),
		 fitWeight(1.), expectedTrueChisq(0.),
		 isClipped(false), itsMatch(nullptr), map(nullptr)  {}

    Detection(Detection &&) = default;

    virtual ~Detection() = default;

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
    // using the RESIDUAL_UNIT from Units.h
    Vector2 residWorld() const;
    // Or the residual mapped back to pixel coordinates
    Vector2 residPix() const;

    // Get a rough estimate of a 1-dimensional world-coordinate error
    double getSigma() const {
      double sigsq = 0.5*(invCov(0,0)+invCov(1,1));
      return sigsq > 0. ? 1./sqrt(sigsq) : 0.;
    }
    void clip() {isClipped=true; fitWeight=0.;}
  };

  class PMDetection: public Detection {
    // Derive a class which is for Gaia or other measurements that come with full 5d
    // constraints.  
    // When filling in this object, one should also fill in the base class Detection
    // fields as if this were a single detection at the field epoch, then this
    // class can be sliced to a Detection and used as such.
  public:
    EIGEN_NEW
    PMSolution pmMean;
    PMCovariance pmInvCov;

    // Calculate full *5d* chisq relative to best prediction of itsMatch,
    // using full fitting covariance but *before* any weighting factor is applied.
    virtual double trueChisq() const;
  };

  class Match {
    // Class for a set of matched Detections, with no PM or parallax freedom
  protected:
    list<unique_ptr<Detection>> elist;
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
    // This flag is set if all weights are 1 or 0
    mutable bool trivialWeights;

    // Here are cached quantities for the overall solution:
    mutable int nFit;	// Number of un-clipped points with non-zero weight in fit
    mutable Vector2 xyMean;   // Best-fit centroid (valid if isSolved).
    mutable Matrix22 centroidF;  // pseudo-Fisher matrix for centroid (2nd deriv for fitting)
    mutable Matrix22 trueCentroidCov; // Expected covariance of centroid after weighted fit.
    mutable int dof;

    // Call this to redo the Fisher matrix if needed (=>set isPrepare)
    virtual void prepare() const;
    // Other state-changing calls are public.

  public:
    EIGEN_NEW
    Match(unique_ptr<Detection> e);

    // Add and remove automatically update itsMatch of the Detection
    void add(unique_ptr<Detection> e);
    void remove(Detection const & e);
    // Remove a Detection from the match given an iterator to it.
    list<unique_ptr<Detection>>::iterator erase(list<unique_ptr<Detection>>::iterator i);
    // Delete all members of match as unmatched.
    void clear();
    // True if a Detection will contribute to chisq:
    static bool isFit(const Detection & e);
    // Number of points that would have nonzero weight in next fit
    int fitSize() const {prepare(); return nFit;}
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

    // Get the expected covariance matrix of centroid after weighted fit
    const Matrix22& getCentroidCov() const {prepare(); return trueCentroidCov;}

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

    // Fitting chisq (not true one) for this match, and largest ratio
    // of true chisq to expected for any detection being fit.
    // Both arguments are updated with info from this match.
    virtual double chisq(int& dofAccum, double& maxDeviateSq,
			 bool dump=false) const;

    typedef list<unique_ptr<Detection>>::iterator iterator;
    typedef list<unique_ptr<Detection>>::const_iterator const_iterator;
    iterator begin() {return elist.begin();}
    iterator end() {return elist.end();}
    const_iterator begin() const {return elist.begin();}
    const_iterator end() const {return elist.end();}
    int size() const {return elist.size();}
  };

  typedef list<unique_ptr<Match>> MCat;
  
  class PMMatch: public Match {
    // Class for a set of matched Detections, with free PM and parallax
  public:
    EIGEN_NEW
    PMMatch(unique_ptr<Detection> e);

    // Set the prior applied to all PMMatches - given in the I/O units
    static void setPrior(double pmPrior, double parallaxPrior);

    virtual void solve() const;  // Recalculate best-fit PM 

    // Increment fitting chisq, beta, and alpha for this match.
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

    // Fitting chisq (not true one) for this match, and largest ratio
    // of true chisq to expected for any detection being fit.
    // Both arguments are updated with info from this match.
    virtual double chisq(int& dofAccum, double& maxDeviateSq,
			 bool dump=false) const;

    // 
    const PMSolution& getPM() const {return pm;}
    const PMCovariance& getInvCovPM() const {prepare(); return pmFisher;} //??
    // Return expected covariance matrix of PM, after weight solution.
    const PMCovariance& getTrueCovPM() const {prepare(); return pmTrueCov;}
    
    // Get predicted position (and inverse covariance) for a Detection.
    // The argument is needed only if there is full PM solution.
    virtual Vector2 predict(const Detection * d = nullptr) const;
    virtual Matrix22 predictFisher(const Detection* d = nullptr)  const;

  protected:
  private:
    static PMSolution priorFisher;  // A diagonal prior on the PM solution for every PMMatch
    
    double parallaxPrior;   // sigma for prior on parallax (radians)

    mutable PMSolution pm;          // Chisq-minimizing position/parallax/pm

    virtual void prepare() const;  // Recalculate and cache, if needed:
    mutable PMCovariance pmFisher;  // Fisher matrix for PM that includes weights
    mutable PMCovariance pmInvFisher;     // Inverse of pmFisher
    mutable PMCovariance pmTrueCov; // Expected cov matrix of PM.
    mutable PMSolution priorMean;   // PM contribution from PMDetection means
  };

  // Class that aligns all coordinates
  class CoordAlign {
  private:
    MCat& mlist;
    PixelMapCollection& pmc;
    double relativeTolerance;
    set<int> frozenParameters;  // Keep track of degenerate parameters
    map<string, set<int>> frozenMaps; // Which atoms have which params frozen
  public:
    CoordAlign(PixelMapCollection& pmc_,
	       MCat& mlist_): mlist(mlist_),
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
		  bool clipEntireMatch=false,
		  bool logging=true);
    // Calculate total fitting chisq.  doReserved same meaning as above.
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
