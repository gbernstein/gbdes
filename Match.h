// Match objects from multiple catalogs together
// And then adjust parameters of their PixelMaps to minimize
// differences from mean world coordinates.
//  Derived from the acs/Match code.
#ifndef MATCH_H
#define MATCH_H

#include <list>
#include <set> 
using std::list;
#include <string>
#include "Std.h"
#include "UseTMV.h"
#include "TMV_SymBand.h"
#include "Bounds.h"
#include "Astrometry.h"
#include "PixelMap.h"
#include "PixelMapCollection.h"

namespace astrometry {

  class AlphaUpdater; // Used only in Match.cpp ???

  class Match;  // Forward declaration

  class Detection {
  public:
    long catalogNumber;
    long objectNumber;
    double xpix;
    double ypix;
    double sigmaPix;
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
    Detection(): itsMatch(0), map(0), isClipped(false) {}
  };
  
  class Match {
  private:
    list<Detection*> elist;
    int nFit;	// Number of un-clipped points with non-zero weight in fit
    bool isReserved;	// Do not contribute to re-fitting if true
  public:
    Match(Detection* e);
    // Add and remove automatically update itsMatch of the Detection
    void add(Detection* e);
    void remove(Detection* e);
    // Mark all members of match as unmatched, or optionally delete them,
    // then empty the list:
    list<Detection*>::iterator erase(list<Detection*>::iterator i,
				     bool deleteDetection=false) {
      if (deleteDetection) delete *i;
      return elist.erase(i);
    }
    // Mark all members of match as unmatched, or optionally delete them,
    // then empty the list:
    void clear(bool deleteDetections=false);
    void remap();  // Remap each point to world coords with current map
    void centroid(double& x, double& y) const;
    void centroid(double& x, double& y, 
		  double& wtx, double &wty) const;
    int size() const {return elist.size();}
    // Number of points that would have nonzero weight in next fit
    int fitSize() const {return nFit;}
    // Update and return count of above:
    int countFit(); 

    // Is this object to be reserved from re-fitting?
    bool getReserved() const {return isReserved;}
    void setReserved(bool b) {isReserved = b;}

    // Returned integer is the DOF count
    int accumulateChisq(double& chisq,
			 DVector& beta,
			//**tmv::SymMatrix<double>& alpha);
			AlphaUpdater& updater);
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

  typedef list<Match*> MCat;

  // Class that aligns all coordinates
  class CoordAlign {
  private:
    list<Match*>& mlist;
    PixelMapCollection& pmc;
    double relativeTolerance;
  public:
    CoordAlign(PixelMapCollection& pmc_,
	       list<Match*>& mlist_): mlist(mlist_),
				      pmc(pmc_), 
				      relativeTolerance(0.001)  {}

    void remap();	// Re-map all Detections using current params
    double fitOnce(bool reportToCerr=true);	// Returns chisq of previous fit, updates params.
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
    void operator()(const DVector& params, double& chisq,
		    DVector& beta, tmv::SymMatrix<double>& alpha);
    void setRelTolerance(double tol) {relativeTolerance=tol;}
    // Return count of useful (un-clipped) Matches & Detections.
    // Count either reserved or non-reserved objects, and require minMatches useful
    // Detections for a valid match:
    void count(long int& mcount, long int& dcount, bool doReserved=false, int minMatches=2) const;
    void count(long int& mcount, long int& dcount, bool doReserved, int minMatches, long catalog) const;
  };

} // namespace astrometry

#endif
