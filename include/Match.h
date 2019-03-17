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

namespace astrometry {

  using linalg::SymmetricUpdater;
  class Match;  // Forward declaration

  // Some typedefs for proper motion
  typedef linalg::SVector<double, 5> PMSolution;
  typedef linalg::SMatrix<double,5,5> PMCovariance;
  enum PM {X0, Y0, VX, VY, PAR};  // Standard order of 5d params
  // Note we will assume that VX, VY are in mas/yr, oriented
  // toward ICRS E and N, respectively.  PAR will be in
  // mas.  Both as per Gaia catalog convention.
  // X0, Y0 are in units of the field's WCS.

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
    double fitWeight;  // For clipping, we remove weight from invCovFit
    // Inverse covariance of original measurement - no weights
    Matrix22 invCovMeas;
    bool isClipped;
    const Match* itsMatch;
    const SubMap* map;

    Detection(): color(astrometry::NODATA), pmProj(nullptr),
		 isClipped(false), itsMatch(nullptr), map(nullptr)  {}
    virtual ~Detection() {if (pmProj) delete pmProj;}

    // Build the projection matrix for this detection.
    void buildProjector(double pmTDB,	// Time in years from PM reference epoch
			const Vector3& xObs, // Observatory position, barycentric ICRS, in AU
			SphericalCoords* fieldProjection, // Projection used for world coords
			double wScale=DEGREE); // units of world coordinates (relative to radians)

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
  };
    
  class Match {
    // Class for a set of matched Detections, with no PM or parallax freedom
  protected:
    list<Detection*> elist;
    int nFit;	// Number of un-clipped points with non-zero weight in fit
    bool isReserved;	// Do not contribute to re-fitting if true
    // True if a Detection will contribute to chisq:
    static bool isFit(const Detection* e);

    virtual void prepare() const;   // Calculate the Fisher matrix and DOF
    mutable bool isPrepared;  // Flag to set if Fisher matrix is current
    mutable Matrix22 centroidF;  // Fisher matrix for centroid (=inverse cov)
    mutable int dof;

    virtual void solve() const;     // Calculate best-fit centroid
    mutable bool isSolved;    // Flag set if centroid/PM is for current coords
    mutable Vector2 xyMean;   // Best-fit centroid (valid if isSolved).
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
    // Number of points that would have nonzero weight in next fit
    int fitSize() const {return nFit;}
    // Recount the number of objects contributing to fit (e.g. after
    // meddling with object weights)
    void countFit();

    // Is this object to be reserved from re-fitting?
    bool getReserved() const {return isReserved;}
    void setReserved(bool b) {isReserved = b;}

    virtual void remap();  // Remap *all* points to world coords with current map
    void forceRecalculation() const {isPrepared = false;}  // Invalidate caches

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
    // Does *not* remap the points.
    virtual bool sigmaClip(double sigThresh,
			   bool deleteDetection=false); 
    void clipAll(); // Mark all detections as clipped

    // Chisq for this match, and largest-sigma-squared deviation
    // 2 arguments are updated with info from this match.
    // Does *not* remap the points.
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

    virtual void remap();  // Remap *all* points to world coords with current map

    // Increment chisq, beta, and alpha for this match.
    // Returned integer is the DOF count.  This *does* remap points being fitted.
    // reuseAlpha=true will skip the incrementing of alpha.
    virtual int accumulateChisq(double& chisq,
				DVector& beta,
				SymmetricUpdater& updater,
				bool reuseAlpha=false);
   
    // sigmaClip returns true if clipped, 
    // and deletes the clipped guy if 2nd arg is true.  
    // Does *not* remap the points
    virtual bool sigmaClip(double sigThresh,
			   bool deleteDetection=false); 

    // Chisq for this match, and largest-sigma-squared deviation
    // 2 arguments are updated with info from this match.
    // Does *not* remap the points.
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

    virtual void solve() const;     // Recalculate and cache, if needed:
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

    void remap();	// Re-map all Detections using current params
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
    void setParams(const DVector& p) {pmc.setParams(p);}
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
