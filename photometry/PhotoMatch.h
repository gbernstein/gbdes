// Match objects from multiple catalogs together
// And then adjust parameters of their PhotoMaps to minimize
// differences from mean world coordinates.

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
#include "PhotoMapCollection.h"

namespace photometry {

  class AlphaUpdater; // Used only in Match.cpp ???

  class Match;  // Forward declaration

  class Detection {
  public:
    long catalogNumber;
    long objectNumber;
    double magIn;
    double sigmaMag;
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
    void remap();  // Remap each point, i.e. make new magOut
    // Mean of un-clipped output mags, optionally with total weight
    void getMean(double& mag) const;
    void getMean(double& mag, double& wt) const;
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

  // A prior on zeropoints of exposures will be forcing reference points into agreement.
  // This class is for the reference points
  class PhotoPriorReferencePoint {
  public:
    double magIn;  // Magnitude assigned to 1 count per second.
    PhotoArguments args;
    double airmass;
    const SubMap* map;	// Photometric map applying to this reference point
  };

  // Class that manifests some prior constraint on the agreement of zeropoints of different
  // exposures; typically one PhotoPrior for a photometric night's data

  class PhotoPrior {
  public:
    PhotoPrior(list<PhotoPriorReferencePoint> _points,
	       double sigma,
	       double zeropoint = 0.,
	       double airmassCoefficient = 0.,
	       double colorCoefficient = 0.,
	       bool freeZeropoint = false,
	       bool freeAirmass = false,
	       bool freeColor = false);
    int nParams() const {return nFree;}
    DVector getParams() const;
    void setParams(const DVector& p);

    // Locations of parameters in global vector
    int startIndex() const {return globalStartIndex;}
    int mapNumber(int iMap) const {return globalMapNumber;}

    void operator()(const DVector& params, double& chisq,
		    DVector& beta, tmv::SymMatrix<double>& alpha);

    int accumulateChisq(double& chisq,
			 DVector& beta,
			//**tmv::SymMatrix<double>& alpha);
			AlphaUpdater& updater);
    // sigmaClip returns true if clipped, and deletes the clipped guy
    bool sigmaClip(double sigThresh);
    // Chisq for this match, also updates dof for free parameters
    double chisq(int& dof, double& maxDeviateSq) const;
    
  private:
    double sigma;	// strength of prior (mags) at each reference point
    list<PhotoPriorReferencePoint> points;
    int nFree;	// Count of number of free parameters among m, a, b.
    int globalStartIndex;	// Index of m/a/b in PhotoMapCollection param vector
    int globalMapNumber;	// map number for resource locking
    bool mIsFree;
    bool aIsFree;
    bool bIsFree;
    double m;	// Mag zeropoint
    double a;	// X-1 airmass coefficient
    double b;	// color coefficient
  }

  // Class that fits to make magnitudes agree
  class PhotoAlign {
  private:
    list<Match*>& mlist;
    PhotoMapCollection& pmc;
    list<PhotoPrior*>& priors;
    double relativeTolerance;
    int nPriorParams;
    int priorMapNumber;	// Map number assigned to all params in the priors
    void countPriorParams();
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
    // is assumed non-photometric).
    int sigmaClipPrior(double sigThresh, bool clipEntirePrior=false);

    // Calculate total chisq.  doReserved same meaning as above.
    void setParams(const DVector& p);
    DVector getParams() const;
    int nParams() const {return pmc.nParams() + nPriorParams;}

    void remap();	// Re-map all Detections using current params
    double fitOnce(bool reportToCerr=true);	// Returns chisq of previous fit, updates params.
    double chisqDOF(int& dof, double& maxDeviate, bool doReserved=false) const;
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
