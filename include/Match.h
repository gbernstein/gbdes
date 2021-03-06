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

  class Detection {
  public:
    long catalogNumber;
    long objectNumber;
    double xpix;
    double ypix;
    double color;
    double sigma;	// Uncertainty, in pixel units
    double xw;
    double yw;
    // weight, nominally inverse sigma squared in each world coord dimen:
    double wtx;
    double wty;
    // Inverse square of the sigma used for sig-clipping:
    double clipsqx;
    double clipsqy;
    bool isClipped;
    const Match* itsMatch;
    const SubMap* map;
  Detection(): itsMatch(nullptr), map(nullptr), isClipped(false), color(astrometry::NODATA) {}
  };
  
  class Match {
  private:
    list<Detection*> elist;
    int nFit;	// Number of un-clipped points with non-zero weight in fit
    bool isReserved;	// Do not contribute to re-fitting if true
    // True if a Detection will contribute to chisq:
    static bool isFit(const Detection* e);
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

    void remap();  // Remap *all* points to world coords with current map
    // Get centroids - these do *not* recalculate xw,yw 
    void centroid(double& x, double& y) const;
    void centroid(double& x, double& y, 
		  double& wtx, double &wty) const;

    // Increment chisq, beta, and alpha for this match.
    // Returned integer is the DOF count.  This *does* remap points being fitted.
    // reuseAlpha=true will skip the incrementing of alpha.
    int accumulateChisq(double& chisq,
			DVector& beta,
			SymmetricUpdater& updater,
			bool reuseAlpha=false);
   
    // sigmaClip returns true if clipped, 
    // and deletes the clipped guy if 2nd arg is true.  
    // Does *not* remap the points, but does call centroid()
    bool sigmaClip(double sigThresh,
		   bool deleteDetection=false); 
    void clipAll(); // Mark all detections as clipped

    // Chisq for this match, and largest-sigma-squared deviation
    // 2 arguments are updated with info from this match.
    // Does *not* remap the points.
    double chisq(int& dof, double& maxDeviateSq) const;

    typedef list<Detection*>::iterator iterator;
    typedef list<Detection*>::const_iterator const_iterator;
    iterator begin() {return elist.begin();}
    iterator end() {return elist.end();}
    const_iterator begin() const {return elist.begin();}
    const_iterator end() const {return elist.end();}
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
    // Fitting routine: returns chisq of previous fit, updates params.
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
