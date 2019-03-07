// Astrometric matching and fitting classes.

#include "Match.h"
#include <list>
using std::list;
#include <set>
using std::set;
#include "StringStuff.h"
#include <fstream>
#include "AstronomicalConstants.h"
#include <algorithm>
#include "Stopwatch.h"

// Need Eigenvalue routines
#ifdef USE_EIGEN
#include "Eigen/Eigenvalues"
#endif

#ifdef _OPENMP
#include <omp.h>
#endif

// #define DEBUG
#include "Marquardt.h"

using namespace astrometry;



bool
Match::isFit(const Detection* e) {
  // A Detection will contribute to fit if it is not clipped,
  // and if it's either a full PMDetection or has finite position cov
  return !(e->isClipped) && (dynamic_cast<const PMDetection*> (e) || e->invCov(0,0)>0.);
}

Match::Match(Detection *e): elist(), nFit(0),dof(0),
			    isReserved(false), isPrepared(false), isSolved(false) {
  add(e);
}

void
Match::add(Detection *e) {
  isPrepared = false;
  elist.push_back(e);
  e->itsMatch = this;
  e->isClipped = false;
  if ( isFit(e)) nFit++;
}

void
Match::remove(Detection* e) {
  isPrepared = false;
  elist.remove(e);  
  e->itsMatch = nullptr;
  if ( isFit(e) ) nFit--;
}

list<Detection*>::iterator
Match::erase(list<Detection*>::iterator i,
	     bool deleteDetection) {
  isPrepared = false;
  if (isFit(*i)) nFit--;
  if (deleteDetection) delete *i;
  return elist.erase(i);
}


void
Match::clear(bool deleteDetections) {
  isPrepared = false;
  for (auto i : elist) {
    if (deleteDetections)
      delete i;
    else
      i->itsMatch = nullptr;
  }
  elist.clear();
  nFit = 0;
}

void
Match::countFit() {
  nFit = 0;
  for (auto i : elist)
    if (isFit(i)) nFit++;
}

void
Match::clipAll() {
  isPrepared = false;
  for (auto i : elist)
    if (i)
      i->isClipped = true;
  nFit = 0;
}

void
Match::prepare() const {
  // Recalculate, if needed, quantities that do not depend on the
  // detection position values.  In this case, just the Fisher
  // matrix of the centroid.
  if (isPrepared) return;
  isSolved = false;   // Change here means centroid solution is invalid
  centroidF.setZero();
  dof = 0;
  for (auto i : elist) {
    if (dynamic_cast<const PMDetection*> (i)) {
      throw AstrometryError("PMDetection included in a non-PM Match");
      // ??? Could treat PMDetection as a single epoch
    }
    if (!isFit(i)) continue;
    centroidF += i->invCov;
    dof += 2;
  }
  dof -= 2;  // Remove 2 parameters in this fit
  isPrepared = true;
}

void
Match::solve() const {
  // Recalculate the centroid, if needed.
  prepare();
  if (isSolved) return;
  isSolved = true;
  if (dof<0) {
    // Do not execute degenerate fit
    xyMean.setZero();
    return;
  }
  Vector2 sumwx(0.);  
  for (auto i : elist) {
    if (!isFit(i)) continue;
    sumwx[0] += i->invCov(0,0) * i->xw + i->invCov(0,1) * i->yw;
    sumwx[1] += i->invCov(1,0) * i->xw + i->invCov(1,1) * i->yw;
  }
  xyMean = centroidF.inverse() * sumwx;
}

void
Match::remap() {
  for (auto i : elist) 
    i->map->toWorld(i->xpix, i->ypix,
		    i->xw, i->yw,
		    i->color);
  isSolved = false;  // Solution has changed
}

///////////////////////////////////////////////////////////
// Parameter-fitting routines
///////////////////////////////////////////////////////////

// A structure that describes a range of parameters
struct iRange {
  iRange(int i=0, int n=0): startIndex(i), nParams(n) {}
  int startIndex;
  int nParams;
}; 

int
Match::accumulateChisq(double& chisq,
		       DVector& beta,
		       SymmetricUpdater& updater,
		       bool reuseAlpha) {
  int nP = beta.size();

  // No contributions to fit for degenerate fit
  if (dof<0) return 0;

  // Update mapping and save derivatives for each detection:
  // Note this invalidates the centroid solution
  isSolved = false;
  vector<DMatrix*> dXYdP(elist.size(), nullptr);
  int ipt=0;
  for (auto i = elist.begin(); i!=elist.end(); ++i, ++ipt) {
    if (!isFit(*i)) continue;
    int npi = (*i)->map->nParams();
    double xw, yw;
    if (npi>0) {
      dXYdP[ipt] = new DMatrix(2,npi);
      (*i)->map->toWorldDerivs((*i)->xpix, (*i)->ypix,
			       xw, yw,
			       *dXYdP[ipt],
			       (*i)->color);
    } else {
      (*i)->map->toWorld((*i)->xpix, (*i)->ypix, 
			 xw, yw,
			 (*i)->color);
    }
    (*i)->xw = xw;
    (*i)->yw = yw;
  }

  // Get new centroid
  solve();
  
  DMatrix dxyMean(2,nP, 0.);
  Matrix22 invCov;
  Vector2 err;
  map<int, iRange> mapsTouched;
  ipt = 0;
  for (auto i = elist.begin(); i!=elist.end(); ++i, ++ipt) {
    if (!isFit(*i)) continue;
    invCov = (*i)->invCov;
    err[0] = (*i)->xw;
    err[1] = (*i)->yw;
    err -= xyMean;

    chisq += err.transpose() * invCov * err; 

    // Accumulate derivatives:
    int istart=0;
    for (int iMap=0; iMap<(*i)->map->nMaps(); iMap++) {
      int ip=(*i)->map->startIndex(iMap);
      int np=(*i)->map->nSubParams(iMap);
      if (np==0) continue;
      int mapNumber = (*i)->map->mapNumber(iMap);
      // Keep track of parameter ranges we've messed with:
      mapsTouched[mapNumber] = iRange(ip,np);

#ifdef USE_TMV
      tmv::ConstVectorView<double> dxy1=dXYdP[ipt]->subMatrix(0,2,istart,istart+np);
#elif defined USE_EIGEN
      DMatrix dxy1 = dXYdP[ipt]->block(0,istart,2,np);
#endif
      beta.subVector(ip, ip+np) -= dxy1.transpose() * invCov * err;

      // Contribution to derivatives of the mean position:
      // (actually \sum_i C_i^{-1} dx_i/dp)
      dxyMean.subMatrix(0,2, ip, ip+np) += invCov * dxy1;

      if (!reuseAlpha) {
	// Increment the alpha matrix
	// First the blocks astride the diagonal:
	updater.update(mapNumber, ip, dxy1, invCov);

	int istart2 = istart+np;
	for (int iMap2=iMap+1; iMap2<(*i)->map->nMaps(); iMap2++) {
	  int ip2=(*i)->map->startIndex(iMap2);
	  int np2=(*i)->map->nSubParams(iMap2);
	  int mapNumber2 = (*i)->map->mapNumber(iMap2);
	  if (np2==0) continue;
#ifdef USE_TMV
	  tmv::ConstVectorView<double> dxy2=dXYdP[ipt]->subMatrix(0,2,istart2,istart2+np2);
#elif defined USE_EIGEN
	  DMatrix dxy2 = dXYdP[ipt]->block(0,istart2,2,np2);
#endif
	  // Now update below diagonal
	  updater.update(mapNumber2, ip2, dxy2,
			 invCov,
			 mapNumber,  ip,  dxy1);
	  istart2+=np2;
	}
      }
      istart+=np;
    } // outer parameter segment loop
    if (dXYdP[ipt]) {delete dXYdP[ipt]; dXYdP[ipt]=nullptr;}
  } // object loop

  if (!reuseAlpha) {
    // Subtract effects of derivatives on mean
    /*  We want to do this, but without touching the entire alpha matrix every time:
	alpha -=  dxyMean^T * centroidF^-1 * dxmean);
    */

    // Do updates parameter block by parameter block
    Matrix22 negCovCentroid = -centroidF.inverse();
    for (auto m1 = mapsTouched.begin();
	 m1!=mapsTouched.end();
	 ++m1) {
      int map1 = m1->first;
      int i1 = m1->second.startIndex;
      int n1 = m1->second.nParams;
      if (n1<=0) continue;
      DMatrix dxy1 = dxyMean.subMatrix(0,2, i1, i1+n1);
      updater.update(map1, i1, dxy1, negCovCentroid);

      auto m2 = m1;
      ++m2;
      for ( ; m2 != mapsTouched.end(); ++m2) {
	int map2 = m2->first;
	int i2 = m2->second.startIndex;
	int n2 = m2->second.nParams;
	if (n2<=0) continue;
	DMatrix dxy2 = dxyMean.subMatrix(0,2,i2, i2+n2);
	updater.update(map2, i2, dxy2,
		       negCovCentroid,
		       map1, i1, dxy1);
      }
    }
  } // Finished putting terms from mean into alpha

  return dof;
}

bool
Match::sigmaClip(double sigThresh,
		 bool deleteDetection) {
  // Only clip the worst outlier at most
  if (dof<0) return false; // Already degenerate, not fitting
  solve();  // Recalculate mean position if out of date
  double maxSq=0.;
  Detection* worst=nullptr;
  Vector2 err;
  for (auto i : elist) {
    if (!isFit(i)) continue;
    err[0] =  i->xw;
    err[1] =  i->yw;
    err -= xyMean;
    // Calculate deviation per DOF
    double devSq = (err.transpose() * i->invCovClip * err);
    devSq /= 2.;
    if ( devSq > sigThresh*sigThresh && devSq > maxSq) {
      // Mark this as point to clip
      worst = i;
      maxSq = devSq;
    }
  }

  if (worst) {
#ifdef DEBUG
    cerr << "clipped " << worst->catalogNumber
	 << " / " << worst->objectNumber 
	 << " at " << worst->xw << "," << worst->yw 
	 << " resid " << setw(6) << (worst->xw-xyMean[0])*DEGREE/ARCSEC
	 << " " << setw(6) << (worst->yw-xyMean[1])*DEGREE/ARCSEC
	 << " resid/sig " << np.sqrt(maxSq)
	 << endl;
#endif
      if (deleteDetection) {
	elist.remove(worst);
	delete worst;
      }	else {
	worst->isClipped = true;
      }
      nFit--;
      isPrepared = false;
      return true;
  } else {
    // Nothing to clip
    return false;
  }
}

double
Match::chisq(int& dofAccum, double& maxDeviateSq) const {
  double chi=0.;
  if (dof<0) return chi; // No info from degenerate fits
  solve();
  Vector2 dxy;
  for (auto i : elist) {
    if (!isFit(i)) continue;
    dxy[0] = i->xw;
    dxy[1] = i->yw;
    dxy -= xyMean;
    double cc = dxy.transpose() * i->invCov * dxy;
    maxDeviateSq = MAX(cc/2. , maxDeviateSq);
    chi += cc;
  }
  dofAccum += dof;
  return chi;
}


/////////////////////////////////////////////////////////////////////
// Routines for matches allowing proper motion and parallax
/////////////////////////////////////////////////////////////////////

PMMatch::PMMatch(Detection *e): Match(e), parallaxPrior(0.), pm(0.) {}
				  
void
PMMatch::prepare() const {
  // Calculate the total PM Fisher matrix and also
  // the sums over the PMDetection priors
  pmFisher.setZero();
  priorMean.setZero();
  priorChisq = 0.;
  dof = 0;
  
  // First put any parallax prior into Fisher
  if (parallaxPrior > 0.) {
    pmFisher(PAR,PAR) = 1./(parallaxPrior*parallaxPrior);
  }

  // Now sum Fisher over all Detections
  PMProjector m(0.);  // (x,y) = m * pm
  m(0,X0) = 1.;
  m(1,Y0) = 1.;
  for (auto i : elist) {
    if (!isFit(i)) continue;
    m.setZero();
    if (auto ii = dynamic_cast<const PMDetection*>(i)) {
      // This detection has full PM covariance of its own
      pmFisher += ii->invCovPM;
      dof += 5;
    } else {
      // Regular single-epoch detection
      i->fillProjector(m,true);
      pmFisher += m.transpose() * i->invCov * m;
      dof += 2;
    }
  }
  pmCov = pmFisher.inverse(); // ??? More stable?
  dof -= 5;  // 5 parameter fit here.
  
  // Go through again and accumulate the PMDetection
  // contributions to chisq
  PMSolution tmp;  // work vector
  for (auto i : elist) {
    if (!isFit(i)) continue;
    m.setZero();
    if (auto ii = dynamic_cast<const PMDetection*>(i)) {
      tmp = ii->invCovPM * ii->pmMean;
      priorChisq += ii->pmMean.transpose() * tmp;
      priorMean += tmp;
    }
  }
  tmp = pmCov * priorMean;
  priorChisq -= priorMean.transpose() * tmp;  // Constant term of chisq
  priorMean = tmp;         // linear term of chisq

  isPrepared = true;
}

void
PMMatch::solve() const {
  prepare();
  PMProjector m(0.);  // (x,y) = m * pm
  m(0,X0) = 1.;
  m(1,Y0) = 1.;
  PMSolution beta(0.);
  Vector2 xy;
  // chisq = pm^T * pmFisher * pm - 2 beta^T * pm + const
  // so pm = pmFisher^-1 * beta
  // Note that our priorMean vector is in fact the contribution
  // to this from the PMDetection priors, so we only
  // accumulate beta from non-PM Detections
  for (auto i : elist) {
    if (!isFit(i)) continue;
    if (dynamic_cast<const PMDetection*>(i)) {
      continue;  // already cached
    } else {
      // Regular single-epoch detection
      i->fillProjector(m,true);
      xy[0] = i->xw;
      xy[1] = i->yw;
      beta += m.transpose() * (i->invCov * xy);
    }
  }
  // Solve the system now
  pm = pmCov * beta + priorMean;
}

double
PMMatch::chisq(int& dofAccum, double& maxDeviateSq) const {
  double chi=0.;
  if (dof<0) return chi; // No info from degenerate fits
  solve();  // Get best PM for current parameters
  PMProjector m(0.);  // (x,y) = m * pm
  m(0,X0) = 1.;
  m(1,Y0) = 1.;
  Vector2 dxy;
  PMSolution dpm;
  double cc;
  for (auto i : elist) {
    if (!isFit(i)) continue;
    if (auto ii = dynamic_cast<const PMDetection*>(i)) {
      dpm = ii->pmMean - pm;
      cc = dpm.transpose() * ii->invCovPM * dpm;
      dof += 5;
      maxDeviateSq = MAX(cc/5. , maxDeviateSq);
    } else {    
      // Regular single-epoch detection
      dxy[0] = i->xw;
      dxy[1] = i->yw;
      i->fillProjector(m,true);
      dxy -= m * pm;
      cc = dxy.transpose() * i->invCov * dxy;
      dof += 2;
      maxDeviateSq = MAX(cc/2. , maxDeviateSq);
    }
    chi += cc;
  }
  dof -= 5;
  return chi;
}

bool
PMMatch::sigmaClip(double sigThresh,
		   bool deleteDetection) {
  // Only clip the worst outlier at most
  if (dof<0) return false; // Already unusable for fitting
  solve();   // Get PM solution
  double maxSq=0.;
  Detection* worst=nullptr;
  PMSolution dpm;
  double devSq;
  Vector2 err;
  for (auto i : elist) {
    if (!isFit(i)) continue;
    if (auto ii = dynamic_cast<const PMDetection*>(i)) {
      dpm = ii->pmMean - pm;
      devSq = dpm.transpose() * ii->invCovPM * dpm;
      devSq /= 5.;
    } else {
      err[0] =  i->xw;
      err[1] =  i->yw;
      err -= predict(i);
      double devSq = err.transpose() * i->invCovClip * err;
      devSq /= 2.;
    }
    if ( devSq > sigThresh*sigThresh && devSq > maxSq) {
      // Mark this as point to clip
      worst = i;
      maxSq = devSq;
    }
  }

  if (worst) {
#ifdef DEBUG
    if (auto ii = dynamic_cast<const PMDetection* worst) {
      cerr << "clipped PM" << worst->catalogNumber
	   << " / " << worst->objectNumber 
	   << " at " << worst->xw << "," << worst->yw 
	   << " resid/sig " << np.sqrt(maxSq) 
	   << endl;
    } else {
      cerr << "clipped " << worst->catalogNumber
	   << " / " << worst->objectNumber 
	   << " at " << worst->xw << "," << worst->yw 
	   << " resid " << setw(6) << (worst->xw-predict(worst)[0])*DEGREE/ARCSEC
	   << " " << setw(6) << (worst->yw-predict(worst)[1])*DEGREE/ARCSEC
	   << " resid/sig " << np.sqrt(maxSq) 
	   << endl;
    }
#endif
      if (deleteDetection) {
	elist.remove(worst);
	delete worst;
      }	else {
	worst->isClipped = true;
      }
      nFit--;
      isPrepared = false;
      return true;
  } else {
    // Nothing to clip
    return false;
  }
}

Vector2
PMMatch::predict(const Detection* d) const {
  if (!d)
    throw AstrometryError("PMMatch::predict called without Detection");
  solve();   // Update full PM solution
  PMProjector m;
  d->fillProjector(m);
  Vector2 out = m * pm;
  return out;
}

Matrix22
PMMatch::predictFisher(const Detection* d)  const {
  if (!d)
    throw AstrometryError("PMMatch::predict called without Detection");
  prepare();
  PMProjector m;
  d->fillProjector(m);
  Matrix22 out = (m * pmCov * m.transpose()).inverse();
  return out;
}

int
PMMatch::accumulateChisq(double& chisq,
			 DVector& beta,
			 SymmetricUpdater& updater,
			 bool reuseAlpha) {
  int nP = beta.size();

  // No contributions to fit for <2 detections:
  if (dof<0) return 0; // No contribution for degenerate fits

  // First loop updates mapping and accumulates derivs.
  // This invalidates the centroid solution
  isSolved = false;
  
  // Save away the derivatives
  vector<DMatrix*> dXYdP(elist.size(), nullptr);

  // And accumulate quantity related to deriv of PM
  DMatrix dPMdP(5,nP,0.);

  // Projector array to use for all points
  PMProjector m(0.);  // (x,y) = m * pm
  m(0,X0) = 1.;
  m(1,Y0) = 1.;
  // And a coordinate holder
  Vector2 xy;

  int ipt=0;  // Index for detections
  for (auto i = elist.begin(); i!=elist.end(); ++i, ++ipt) {
    if (!isFit(*i)) continue;
    if (dynamic_cast<const PMDetection*>(*i)) continue;

    // Regular single-epoch detection
    (*i)->fillProjector(m,true);
    int npi = (*i)->map->nParams();
    double xw, yw;
    if (npi>0) {
      dXYdP[ipt] = new DMatrix(2,npi);
      (*i)->map->toWorldDerivs((*i)->xpix, (*i)->ypix,
			       xw, yw,
			       *dXYdP[ipt],
			       (*i)->color);
    } else {
      // Remap but no derivs needed
      (*i)->map->toWorld((*i)->xpix, (*i)->ypix, 
			 xw, yw,
			 (*i)->color);
    }

    (*i)->xw = xw;
    (*i)->yw = yw;
    // For calculating best-fit PM
    xy[0] = (*i)->xw;
    xy[1] = (*i)->yw;
    beta += m.transpose() * (*i)->invCov * xy;
  }
  // Solve the system now
  pm = pmCov * beta + priorMean;
  isSolved = true;
  
  // Next loop through detections will
  // accumulate chisq, derivs of PM,
  // and chisq derivs that depend on
  // only one detection.

  Matrix22 invCov;
  Vector2 xyErr;
  PMSolution pmErr;
  PMProjector cInvM;
  
  map<int, iRange> mapsTouched;
  ipt = 0;
  for (auto i = elist.begin(); i!=elist.end(); ++i, ++ipt) {
    if (!isFit(*i)) continue;

    if (auto ii = dynamic_cast<const PMDetection*> (*i)) {
      // PMDetection has no parameter dependence, just chisq
      pmErr = ii->pmMean - pm;
      chisq += pmErr.transpose() * ii->invCovPM * pmErr;
    } else {
      // First calculate chisq contribution
      xyErr[0] = (*i)->xw;
      xyErr[1] = (*i)->yw;
      xyErr -= m * pm;
      chisq += xyErr.transpose() * (*i)->invCov * xyErr; 

      // Then contributions to chisq derivs
      int istart=0;  // running index into the Detection's params
      for (int iMap=0; iMap<(*i)->map->nMaps(); iMap++) {
	int ip=(*i)->map->startIndex(iMap);
	int np=(*i)->map->nSubParams(iMap);
	if (np==0) continue;
	int mapNumber = (*i)->map->mapNumber(iMap);
	// Keep track of parameter ranges we've messed with:
	mapsTouched[mapNumber] = iRange(ip,np);

#ifdef USE_TMV
	tmv::ConstVectorView<double> dxy1=dXYdP[ipt]->subMatrix(0,2,istart,istart+np);
#elif defined USE_EIGEN
	DMatrix dxy1 = dXYdP[ipt]->block(0,istart,2,np);
#endif
	// Contribution to first derivative of chisq:
	DMatrix cInvD = (*i)->invCov * dxy1;
	beta.subVector(ip, ip+np) -= xyErr.transpose() * cInvD;

	// Contribution to derivatives of PM:
	dPMdP.subMatrix(0, 5, ip, ip+np) += m.transpose() * cInvD;

	if (!reuseAlpha) {
	  // Increment the alpha matrix
	  // First the blocks astride the diagonal:
	  updater.update(mapNumber, ip, dxy1, (*i)->invCov);

	  // Then loop through all blocks below diagonal
	  int istart2 = istart+np;
	  for (int iMap2=iMap+1; iMap2<(*i)->map->nMaps(); iMap2++) {
	    int ip2=(*i)->map->startIndex(iMap2);
	    int np2=(*i)->map->nSubParams(iMap2);
	    int mapNumber2 = (*i)->map->mapNumber(iMap2);
	    if (np2==0) continue;
#ifdef USE_TMV
	    tmv::ConstVectorView<double> dxy2=dXYdP[ipt]->subMatrix(0,2,istart2,istart2+np2);
#elif defined USE_EIGEN
	    DMatrix dxy2 = dXYdP[ipt]->block(0,istart2,2,np2);
#endif
	    // Now update below diagonal
	    updater.update(mapNumber2, ip2, dxy2,
			   (*i)->invCov,
			   mapNumber,  ip,  dxy1);
	    istart2+=np2;
	  }
	} // End inner loop of dXYdP^T * invC * dXYdP for t
	istart+=np;
      } // outer parameter segment loop
      if (dXYdP[ipt]) {delete dXYdP[ipt]; dXYdP[ipt]=nullptr;}
    } // endif for PMDetections
  } // object loop

  if (!reuseAlpha) {
    // Subtract terms with derivatives of PM,
    /* but without touching the entire alpha matrix every time:
       alpha -=  dPMdP^T * pmCov * dPMdP);
    */

    // Do updates parameter block by parameter block
    PMCovariance negCov = -pmCov;
    for (auto m1 = mapsTouched.begin();
	 m1!=mapsTouched.end();
	 ++m1) {
      int map1 = m1->first;
      int i1 = m1->second.startIndex;
      int n1 = m1->second.nParams;
      if (n1<=0) continue;
      DMatrix dxy1 = dPMdP.subMatrix(0,5, i1, i1+n1);
      updater.update(map1, i1, dxy1, negCov);

      auto m2 = m1;
      ++m2;
      for ( ; m2 != mapsTouched.end(); ++m2) {
	int map2 = m2->first;
	int i2 = m2->second.startIndex;
	int n2 = m2->second.nParams;
	if (n2<=0) continue;
	DMatrix dxy2 = dPMdP.subMatrix(0,5,i2, i2+n2);
	updater.update(map2, i2, dxy2,
		       negCov,
		       map1, i1, dxy1);
      }
    }
  } // Finished putting terms from mean into alpha

  return dof;
}



/////////////////////////////////////////////////////////////////////
// Parameter adjustment operations
/////////////////////////////////////////////////////////////////////
void
CoordAlign::operator()(const DVector& p, double& chisq,
		       DVector& beta, DMatrix& alpha,
		       bool reuseAlpha) {
  int nP = pmc.nParams();
  Assert(p.size()==nP);
  Assert(beta.size()==nP);
  Assert(alpha.rows()==nP);
  setParams(p);
  double newChisq=0.;
  beta.setZero();
  if (!reuseAlpha) alpha.setZero();
  int matchCtr=0;

  const int NumberOfLocks = 2000;
  SymmetricUpdater updater(alpha, pmc.nFreeMaps(), NumberOfLocks);

#ifdef _OPENMP
  const int chunk=100;
  vector<Match*> vi(mlist.size());

#pragma omp parallel reduction(+:newChisq) 
  {
#pragma omp single
    {
      // Each thread will be handed chunk matches to process.
      // It will take nRounds chunks per thread, plus possibly a
      // further chunk for some of the threads, to complete.
      // I will reorder the matches in the vi vector that is
      // to be split into threads, because we want the threads
      // working on widely separated parts of the match list
      // to make it less likely they are working on the same catalogs
      // concurrently which would result in collisions writing to the
      // alpha array.
      const int nThreads=omp_get_num_threads();
      int nRounds = mlist.size()/chunk/nThreads;
      int iThread=0;
      int iChunk=0;
      int iRound=0;
      int j=0;  
      for (auto i : mlist) {
	vi[j] = i;
	j++;
	if (j%chunk==0 && iRound < nRounds) {
	  // We have just finished writing a chunk's worth of
	  // matches to the array.  Next chunk will be for a different
	  // thread.

	  // Note that if we have advanced to the last (partial) round,
	  // we just assign the matches consecutively.
	  iThread++;
	  if (iThread >= nThreads) {
	    // Now on to the next round of threads
	    iRound++;
	    iThread=0;
	  }
	  // Next thread will get matches that are nRounds*chunk removed
	  // from previous one.
	  if (iRound<nRounds) j = (iThread*nRounds+iRound)*chunk;
	}
      }
    }

  // Each thread accumulates its own beta
    DVector newBeta(beta.size(), 0.);
#pragma omp for schedule(dynamic,chunk) 
    for (int i=0; i<vi.size(); i++) {
      auto m = vi[i];
      if ( m->getReserved() ) continue;	//skip reserved objects
      m->accumulateChisq(newChisq, newBeta, updater, reuseAlpha);
    }
#pragma omp critical(beta)
    beta += newBeta;
  }
#else
  // Without OPENMP, just loop through all matches:
  for (auto i : mlist) {
    Match* m = i;
    if (matchCtr%10000==0) cerr << "# accumulating chisq at match # " 
				<< matchCtr 
				<< endl;
    matchCtr++;
    if ( m->getReserved() ) continue;	//skip reserved objects
    m->accumulateChisq(newChisq, beta, updater, reuseAlpha);
  }
#endif
  chisq = newChisq;

  if (!reuseAlpha) {
    // Code to spot unconstrained parameters:
    set<string> newlyFrozenMaps;
    for (int i = 0; i<alpha.rows(); i++) {
      bool blank = true;
      for (int j=0; j<alpha.cols(); j++) 
	if ( (i>=j && alpha(i,j)!=0.) || (i<j && alpha(j,i)!=0.)) {
	  // Note that we are checking in lower triangle only
	  blank = false;
	  break;
	}
      if (blank) {
	string badAtom = pmc.atomHavingParameter(i);
	// Is it a newly frozen parameter?
	if (!frozenMaps.count(badAtom) || !frozenMaps[badAtom].count(i))
	  newlyFrozenMaps.insert(badAtom);
	// Add to (or make) a list of the frozen parameters in this atom
	frozenMaps[badAtom].insert(i);
	alpha(i,i) = 1.;
	beta[i] = 0.;
      } else {
	// Something is weird if a frozen parameter is now constrained
	if (frozenParameters.count(i)>0) {
	  string badAtom = pmc.atomHavingParameter(i);
	  FormatAndThrow<AstrometryError>() << "Frozen parameter " << i
					    << " in map " << badAtom
					    << " became constrained??";
	}
      }
    } // End alpha row loop

    for (auto badAtom : newlyFrozenMaps) {
      // Print message about freezing parameters
      int startIndex, nParams;
      pmc.parameterIndicesOf(badAtom, startIndex, nParams);
      cerr << "Freezing " << frozenMaps[badAtom].size()
	   << " of " << nParams
	   << " parameters in map " << badAtom;
      if (frozenMaps[badAtom].size() < nParams) {
	// Give the parameter indices
	cerr << " (";
	for (auto i : frozenMaps[badAtom])
	  cerr << i - startIndex << " ";
	cerr << ")";
      }
      cerr << endl;
    }
  } // End degenerate parameter check
}

double
CoordAlign::fitOnce(bool reportToCerr, bool inPlace) {
  DVector p = getParams();
  // First will try doing Newton iterations, keeping a fixed Hessian.
  // If it increases chisq or takes too long to converge, we will 
  // go into the Marquardt solver.

  {
    // Get chisq, beta, alpha at starting parameters
    double oldChisq = 0.;
    int nP = pmc.nParams();
    DVector beta(nP, 0.);
    DMatrix alpha(nP, nP, 0.);

    Stopwatch timer;
    timer.start();
    (*this)(p, oldChisq, beta, alpha);
    timer.stop();
    if (reportToCerr) cerr << "..fitOnce alpha time " << timer << endl;
    timer.reset();
    timer.start();

    // Precondition alpha if wanted
    bool precondition = true;
    int N = alpha.cols();
    DVector ss(N, 1.);  // Values to scale parameters by
    if (precondition) {
      for (int i=0; i<N; i++) {
	if (alpha(i,i)<0.) {
	  cerr << "Negative alpha diagonal " << alpha(i,i)
	       << " at " << i << endl;
	  exit(1);
	} 
	if (alpha(i,i) > 0.)
	  ss[i] = 1./sqrt(alpha(i,i));
	// Scale row / col of lower triangle, hitting diagonal twice
	for (int j=0; j<=i; j++) 
	  alpha(i,j) *= ss(i);
	for (int j=i; j<N; j++) 
	  alpha(j,i) *= ss(i);
      }
    }

    // Do the Cholesky decomposition, set flag if it fails
    bool choleskyFails = false;
#ifdef USE_TMV
    // Make a symmetric view into alpha for Cholesky:
    auto symAlpha = tmv::SymMatrixViewOf(alpha,tmv::Lower);
    symAlpha.divideUsing(tmv::CH);
    if (inPlace) symAlpha.divideInPlace();
    try {
      symAlpha.setDiv();
    } catch (tmv::Error& m) {
      choleskyFails=true;
    }
#elif defined USE_EIGEN
    if (inPlace)
      throw AstrometryError("Do not currently support in-place Cholesky for Eigen");
    // The difficulty is that the type of an in-place LLT is different from
    // the type for a "normal" LLT.  I can't figure out how to allow both
    // possibilities without making two branches of the whole iteration code.
    auto llt = alpha.llt();
    if (llt.info()==Eigen::NumericalIssue)
      choleskyFails = true;
#endif
    
    // If the Cholesky decomposition failed for non-pos-def matrix, then
    // as a diagnostic we will do an SVD and report the nature of degeneracies.
    // Then exit with an error.
    // And since alpha is symmetric, we will use eigenvector routines to do the
    // SVD, which is what we want anyway to find the non-pos-def values

    if (choleskyFails) {
      cerr << "Caught exception" << endl;
      if (inPlace) {
	cerr << "Cannot describe degeneracies while dividing in place" << endl;
	exit(1);
      }
      int N = alpha.cols();
      set<int> degen;
      DMatrix U(N,N);
      DVector S(N);

#ifdef USE_TMV
      tmv::Eigen(symAlpha, U, S);
#elif defined USE_EIGEN
      {
	Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eig(alpha);
	U = eig.eigenvectors();
	S = eig.eigenvalues();
      }
#endif
      // Both packages promise to return eigenvalues in increasing
      // order, but let's not depend on that.  Report largest/smallest
      // abs values of eval's, and print them all
      int imax = S.size()-1; // index of largest, smallest eval
      double smax=abs(S[imax]);
      int imin = 0;
      double smin=abs(S[imin]);
      for (int i=0; i<U.cols(); i++) {
	cerr << i << " Eval: " << S[i] << endl;
	double s= abs(S[i]);
	if (s>smax) {
	  smax = s;
	  imax = i;
	}
	if (s<smin) {
	  smin = s;
	  imin = i;
	}
	if (S[i]<1e-6) degen.insert(i);
      }
      cerr << "Largest abs(eval): " << smax << endl;
      cerr << "Smallest abs(eval): " << smin << endl;
      degen.insert(imin);
      // Find biggest contributors to non-positive (or marginal) eigenvectors
      const int ntop=MIN(N,20);
      for (int isv : degen) {
	  cerr << "--->Eigenvector " << isv << " eigenvalue " << S(isv) << endl;
	  // Find smallest abs coefficient
	  int imin = 0;
	  for (int i=0; i<U.rows(); i++)
	    if (abs(U(i,isv)) < abs(U(imin,isv)))
	      imin = i;
	  vector<int> top(ntop,imin);
	  for (int i=0; i<U.rows(); i++) {
	    for (int j=0; j<ntop; j++) {
	      if (abs(U(i,isv)) >= abs(U(top[j],isv))) {
		// Push smaller entries to right
		for (int k=ntop-1; k>j; k--)
		  top[k] = top[k-1];
		top[j] = i;
		break;
	      }
	    }
	  }
	  for (int j : top) {
	    string badAtom = pmc.atomHavingParameter(j);
	    int startIndex, nParams;
	    pmc.parameterIndicesOf(badAtom, startIndex, nParams);
	    cerr << "Coefficient " << U(j, isv) 
		 << " at parameter " << j 
		 << " Map " << badAtom 
		 << " " << j - startIndex << " of " << nParams
		 << endl;
	  }
	}
	exit(1);
    }


    // Now attempt Newton iterations to solution, with fixed alpha
    const int MAX_NEWTON_STEPS = 8;
    int newtonIter = 0;
    for (int newtonIter = 0; newtonIter < MAX_NEWTON_STEPS; newtonIter++) {
      if (precondition) {
	beta = ElemProd(beta,ss);
      }

#ifdef USE_TMV
      beta /= symAlpha;
	
#elif defined USE_EIGEN
      beta = llt.solve(beta);
#endif
      if (precondition) {
	beta = ElemProd(beta,ss);
      }
	
      timer.stop();
      if (reportToCerr) cerr << "..solution time " << timer << endl;
      timer.reset();
      timer.start();
      DVector newP = p + beta;
      setParams(newP);
      // Get chisq at the new parameters
      remap();
      int dof;
      double maxDev;
      double newChisq = chisqDOF(dof, maxDev);
      timer.stop();
      //**if (reportToCerr) {
      if (true) {
	cerr << "....Newton iteration #" << newtonIter << " chisq " << newChisq 
	     << " / " << dof 
	     << " in time " << timer << " sec"
	     << endl;
      }
      timer.reset();
      timer.start();

      // Give up on Newton if chisq went up non-trivially
      if (newChisq > oldChisq * 1.0001) break;
      else if ((oldChisq - newChisq) < oldChisq * relativeTolerance) {
	// Newton has converged, so we're done.
	return newChisq;
      }
      // Want another Newton iteration, but keep alpha as before
      p = newP;
      (*this)(p, oldChisq, beta, alpha, true);
    }
    // If we reach this place, Newton is going backwards or nowhere, slowly.
    // So just give it up.
  }

  // ??? Signal that alpha should be fixed for all iterations of Marquardt?
  Marquardt<CoordAlign> marq(*this);
  marq.setRelTolerance(relativeTolerance);
  marq.setSaveMemory();
  double chisq = marq.fit(p, DefaultMaxIterations, reportToCerr);
  setParams(p);
  remap();
  return chisq;
}

void
CoordAlign::remap() {
  for (auto i : mlist)
    i->remap();
}

int
CoordAlign::sigmaClip(double sigThresh, bool doReserved, bool clipEntireMatch) {
  int nclip=0;
  // ??? parallelize this loop!
  cerr << "## Sigma clipping...";
  Stopwatch timer;
  timer.start();

  for (auto i : mlist) {
    // Skip this one if it's reserved and doReserved=false,
    // or vice-versa
    if (doReserved ^ i->getReserved()) continue;
    if ( i->sigmaClip(sigThresh)) {
      nclip++;
      if (clipEntireMatch) i->clipAll();
    }
  }
  timer.stop();
  cerr << " done in " << timer << " sec" << endl;
  return nclip;
}

double
CoordAlign::chisqDOF(int& dof, double& maxDeviate, 
		     bool doReserved) const {
  dof=0;
  maxDeviate=0.;
  double chisq=0.;
  // ??? parallelize this loop!
  for (auto i : mlist) {
    // Skip this one if it's reserved and doReserved=false,
    // or vice-versa
    if (doReserved ^ i->getReserved()) continue;
    chisq += i->chisq(dof, maxDeviate);
  }
  maxDeviate = sqrt(maxDeviate);
  if (!doReserved) dof -= pmc.nParams();
  return chisq;
}

void 
CoordAlign::count(long int& mcount, long int& dcount, 
		  bool doReserved, int minMatches) const {
  mcount = 0;
  dcount = 0;
  for (auto i : mlist) {
    if ((i->getReserved() ^ doReserved) 
	|| i->fitSize() < minMatches) {
      continue;
    }
    mcount++;
    dcount+=i->fitSize();
  }
}

void
CoordAlign::count(long int& mcount, long int& dcount,
                  bool doReserved, int minMatches, long catalog) const {
  mcount = 0;
  dcount = 0;
  for (auto i : mlist) {
    if ((i->getReserved() ^ doReserved)
        || i->fitSize() < minMatches) continue;

    int ddcount=0;
    for(auto d : *i) {
      if(!(d->isClipped) && d->catalogNumber==catalog)
	ddcount++;
    }

    if(ddcount>0) {
      mcount++;
      dcount+=ddcount;
    }
  }
}
