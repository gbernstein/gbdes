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

// Static member of the PMMatch object:
PMSolution PMMatch::priorFisher(0.);

// And its setter
void PMMatch::setPrior(double pmPrior, double parallaxPrior) {
  priorFisher.setZero();
  priorFisher[PM::VX] = priorFisher[PM::VY] = pow(pmPrior*PM_UNIT/(WCS_UNIT/TDB_UNIT), -2.);
  priorFisher[PM::PAR] = pow(parallaxPrior*PARALLAX_UNIT / WCS_UNIT, -2.);
}


void
Detection::buildProjector(double pmTDB,	       // Time in years from PM reference epoch
			  const Vector3& xObs, // Observatory position, barycentric ICRS, in AU
			  SphericalCoords* fieldProjection) {

  if (pmProj) {
    pmProj->setZero();
  } else {
    pmProj = new PMProjector(0.);
  }
  PMProjector& m = *pmProj;  // Just making a shorter name
  m(0,X0) = 1.;
  m(1,Y0) = 1.;

  // Get the ICRS direction of this detection, and derivative between
  // world coords and ICRS values
  fieldProjection->setLonLat(xw*WCS_UNIT,yw*WCS_UNIT);
  SphericalICRS icrs(*fieldProjection);  // Convert to ICRS

  Matrix22 dWdICRS;  // Coordinate derivative - radian units
  fieldProjection->convertFrom(icrs, dWdICRS);
  // Apply cos(dec) terms to ICRS RA derivatives:
  double ra,dec;
  icrs.getRADec(ra,dec);
  dWdICRS.col(0) /= cos(dec);
  
  // Use Astrometry classes to get observatory position components
  // to E/N of line of sight
  // First make a reference frame aligned to line of sight
  ReferenceFrame frame( CartesianICRS(0.,0.,0.),
			Orientation(icrs));
  // Then convert observatory position to this frame
  CartesianICRS obs(xObs);
  CartesianCustom cc(obs, frame);
  // Pull out components directed along ICRS E and N:
  Vector2 icrsParallax = cc.getVector().subVector(0,2);
  // Transform into WCS - parallax is already in WCS_UNIT
  Vector2 wcsParallax = -(dWdICRS * icrsParallax);
  m.col(PAR) = wcsParallax;

  // Convert proper motion from ICRS directions to WCS
  Matrix22 pm = dWdICRS * pmTDB;
 // (pm should already be in WCS_UNIT / TDB_UNIT)
  m.subMatrix(0,2,VX,VY+1) = pm;
}

Vector2
Detection::residWorld() const {
  Vector2 out;
  Assert(itsMatch);
  itsMatch->solve();
  out[0] = xw;
  out[1] = yw;
  out -= itsMatch->predict(this);
  return out * WCS_UNIT / RESIDUAL_UNIT;
}

Vector2
Detection::residPix() const {
  Assert(itsMatch);
  itsMatch->solve();
  Vector2 xyMean = itsMatch->predict(this);
  double xMeanPix=xpix, yMeanPix=ypix;
  Assert(map);
  try {
    map->toPix(xyMean[0], xyMean[1],
	       xMeanPix, yMeanPix, color);
  } catch (astrometry::AstrometryError& e) {
    std::cerr << "WARNING: toPix failure in map " << map->getName()
	 << " at world coordinates (" << xyMean[0] << "," << xyMean[1] << ")"
	 << " approx pixel coordinates (" << xpix << "," << ypix << ")"
	 << endl;
    // Just return zero residual
    return Vector2(0.);
  }
  Vector2 out;
  out[0] = xpix - xMeanPix;
  out[1] = ypix - yMeanPix;
  return out;
}

double
Detection::trueChisq() const {
  // Calculate chisq relative to best prediction of itsMatch,
  // using full fitting covariance but *before* any weighting factor is applied.
  Assert(itsMatch);
  itsMatch->solve();
  Vector2 xyW;
  xyW[0] = xw;
  xyW[1] = yw;
  xyW -= itsMatch->predict(this);
  double chisq = xyW.transpose() * invCov * xyW;
  return chisq;
}

double
PMDetection::trueChisq() const {
  // Get full 5d chisq for a PMDetection
  auto pmMatch = dynamic_cast<const PMMatch*> (itsMatch);
  Assert(pmMatch);
  // Have the match solve for best fit if not in that state
  pmMatch->solve();
  // Now get chisq
  PMSolution dpm = pmMean - pmMatch->getPM();
  double chisq = dpm.transpose() * pmInvCov * dpm;
  return chisq;
}

bool
Match::isFit(const Detection & e) {
  // A Detection will contribute to fit if it is not clipped,
  // and if it's either a full PMDetection or has finite position cov
  return e.fitWeight>0.;
}

Match::Match(unique_ptr<Detection> e): elist(), nFit(0),dof(0),
			    isReserved(false), isPrepared(false),
			    isMappedFit(false), isMappedAll(false),
			    isSolved(false), trivialWeights(false) {
  add(std::move(e));
}

void
Match::add(unique_ptr<Detection> e) {
  isPrepared = false;
  elist.push_back(std::move(e));
  e->itsMatch = this;
  e->isClipped = false;
  isMappedFit = false;
  isMappedAll = false;
  isSolved = false;
}

void
Match::remove(Detection const & e) {
  isPrepared = false;
  isSolved = false;
  elist.remove_if([&e](unique_ptr<Detection> const & d) { return &e == d.get(); });
  // no need to unset itsMatch, because the Detection has been deleted.
}

list<unique_ptr<Detection>>::iterator
Match::erase(list<unique_ptr<Detection>>::iterator i) {
  isPrepared = false;
  isSolved = false;
  return elist.erase(i);
}


void
Match::clear() {
  isPrepared = false;
  elist.clear();
  nFit = 0;
  isMappedFit = true;
  isMappedAll = true;
  isSolved = false;
}

void
Match::clipAll() {
  isPrepared = false;
  isSolved = true;  // No such thing as solution anymore.
  for (auto const & i : elist)
    if (i)
      i->clip();
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
  // centroidF is the matrix used in solution for centroid.
  // But we will also calculate the actual expected variance
  // of the centroid, which will differ if the weights are not
  // all equal (or zero).
  Matrix22 centroidCov(0.);
  dof = 0;
  nFit = 0;
  trivialWeights = true;
  for (auto const & i : elist) {
    if (dynamic_cast<const PMDetection*> (i.get())) {
      throw AstrometryError("PMDetection included in a non-PM Match");
      // Any PMDetection inputs should have been converted to
      // fixed-epoch Detection for a non-PM match.
    }
    if (!isFit(*i)) continue;
    dof += 2;
    nFit++;
    if (i->fitWeight != 1.) {
      trivialWeights = false;  // Need explicit weights in formulae
      centroidF += i->invCov * i->fitWeight;
    } else {
      centroidF += i->invCov;
    }
      
  }
  dof -= 2;  // Remove 2 parameters in this fit

  if (dof<0) {
    // Degenerate fit, will never yield a solution of interest
    for (auto const & i : elist)
      i->expectedTrueChisq = 0.;
    centroidF.setZero();
    trueCentroidCov.setZero();
    isPrepared = true;
    nFit = 0;
    return;
  }
	
  Matrix22 invF = centroidF.inverse();
  // Calculate expected centroid 
  if (trivialWeights) {
    trueCentroidCov = invF;
  } else {
    // Altered expectation for cov because of weighting
    // From notes of 28 Apr 2019,
    Matrix22 tmp(0.);
    for (auto const & i : elist) {
      if (!isFit(*i)) continue;
      tmp += i->invCov * (i->fitWeight*i->fitWeight);
    }
    trueCentroidCov = invF * tmp * invF;
  }

  // Now get expected chisq for each detection, and true centroid cov
  // From notes of 28 Apr 2019,
  // <chisq_i> = 2 - 2 * w_i * Tr( invCov_i * invF) + Tr(invCov_i * centroidCov)
  // where invF is the inverse of centroidF
  for (auto const & i : elist) {
    if (i->invCov(0,0)<=0.) {
      // No error on this detection so no chisq
      i->expectedTrueChisq = 0.;
    } else if (i->fitWeight==0.) {
      // Has uncertainties but was not included in the fit
      i->expectedTrueChisq = 2 + traceATB(i->invCov, trueCentroidCov);
    } else if (trivialWeights) {
      // included in fit, unit weight gauranteed
      i->expectedTrueChisq = 2 - traceATB(i->invCov, trueCentroidCov);
    } else {
      // Include non-trivial weight
      i->expectedTrueChisq = 2 - (2*i->fitWeight)* traceATB(i->invCov, invF)
	+ traceATB(i->invCov, trueCentroidCov);
    }
  }
  // Calculate expected chisq for each Detection:
  isPrepared = true;
}

void
Match::solve() const {
  // Recalculate the centroid, if needed.
  prepare();     // Make sure Fisher matrices are current
  remap(false);  // Make sure world coords of fitted Detections are current
  if (isSolved) return;
  isSolved = true;
  if (dof<0) {
    // Do not execute degenerate fit
    xyMean.setZero();
    return;
  }
  Vector2 sumwx(0.);
  Vector2 xy;
  for (auto const & i : elist) {
    if (!isFit(*i)) continue;
    xy[0] = i->xw;
    xy[1] = i->yw;
    sumwx += (i->invCov * xy) * i->fitWeight;
  }
  xyMean = centroidF.inverse() * sumwx;
}

void
Match::remap(bool doAll) const { 
  // No need to do anything if all is current
  if (isMappedAll) return;
  // Or if Fit are current and we only need them
  if (!doAll && isMappedFit) return;
  
  for (auto const & i : elist) {
    if (isFit(*i)) {
      if (!isMappedFit) {
	/**/ if (!i->map) std::cerr << "Missing map at " << i->catalogNumber << endl;
	i->map->toWorld(i->xpix, i->ypix,
			i->xw, i->yw,
			i->color);
	isSolved = false;  // Solution has changed
      }
    } else if (doAll) {
      // Only remap non-fits if requested
      i->map->toWorld(i->xpix, i->ypix,
		      i->xw, i->yw,
		      i->color);
    }
  }
  // Update state
  isMappedFit = true;
  if (doAll) isMappedAll = true;
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

  prepare();

  // No contributions to fit for degenerate/perfect fit
  if (dof<=0) return 0;

  // Update mapping and save derivatives for each detection:
  vector<DMatrix*> dXYdP(elist.size(), nullptr);
  int ipt=0;
  for (auto i = elist.begin(); i!=elist.end(); ++i, ++ipt) {
    if (!isFit(**i)) continue;
    int npi = (*i)->map->nParams();
    double xw, yw;
    if (npi>0) {
      dXYdP[ipt] = new DMatrix(2,npi);
      (*i)->map->toWorldDerivs((*i)->xpix, (*i)->ypix,
			       (*i)->xw, (*i)->yw,
			       *dXYdP[ipt],
			       (*i)->color);
    } else if (!isMappedFit) {
      // If there are no free parameters, we only
      // need to remap if the maps have changed.
      (*i)->map->toWorld((*i)->xpix, (*i)->ypix, 
			 (*i)->xw, (*i)->yw,
			 (*i)->color);
    }
  }
  isMappedFit = true;  // This is now true if it wasn't before.

  // Get new centroid
  solve();
  
  DMatrix dxyMean(2,nP, 0.);
  Matrix22 invCovW;  // inverse covariance * weight
  Vector2 err;
  map<int, iRange> mapsTouched;
  ipt = 0;
  for (auto i = elist.begin(); i!=elist.end(); ++i, ++ipt) {
    if (!isFit(**i)) continue;
    invCovW = (*i)->invCov * (*i)->fitWeight;
    err[0] = (*i)->xw;
    err[1] = (*i)->yw;
    err -= xyMean;  
    chisq += err.transpose() * invCovW * err; 

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
      beta.subVector(ip, ip+np) -= dxy1.transpose() * invCovW * err;

      // Contribution to derivatives of the mean position:
      // (actually \sum_i C_i^{-1} dx_i/dp)
      dxyMean.subMatrix(0,2, ip, ip+np) += invCovW * dxy1;

      if (!reuseAlpha) {
	// Increment the alpha matrix
	// First the blocks astride the diagonal:
	updater.update(mapNumber, ip, dxy1, invCovW);

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
			 invCovW,
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
  solve();  // Recalculate mean position if out of date
  if (dof<=0) return false; // Already degenerate, not fitting

  // Only clip the worst outlier at most
  double maxSq=0.;
  iterator worst = elist.end();
  Vector2 err;
  for (auto i = elist.begin(); i != elist.end(); ++i) {
    if (!isFit(**i)) continue;
    err[0] =  (*i)->xw;
    err[1] =  (*i)->yw;
    err -= xyMean;
    // Calculate deviation relative to expected chisq
    double devSq = (err.transpose() * (*i)->invCov * err);
    devSq /= (*i)->expectedTrueChisq;
    if ( devSq > sigThresh*sigThresh && devSq > maxSq) {
      // Mark this as point to clip
      worst = i;
      maxSq = devSq;
    }
  }

  if (worst != elist.end()) {
#ifdef DEBUG
    cerr << "clipped " << (*worst)->catalogNumber
	 << " / " << (*worst)->objectNumber 
	 << " at " << (*worst)->xw << "," << (*worst)->yw 
	 << " resid " << setw(6) << ((*worst)->xw-xyMean[0])*WCS_UNIT / RESIDUAL_UNIT;
	 << " " << setw(6) << ((*worst)->yw-xyMean[1])**WCS_UNIT / RESIDUAL_UNIT;
	 << " resid/sig " << np.sqrt(maxSq)
	 << endl;
#endif
      if (deleteDetection) {
	elist.erase(worst);
      }	else {
	(*worst)->clip();
      }
      isPrepared = false;
      isSolved = false;  // Note that world coords are still valid
      return true;
  } else {
    // Nothing to clip
    return false;
  }
}

double
Match::chisq(int& dofAccum, double& maxDeviateSq,
	     bool dump) const {
  solve();
  double chi=0.;
  if (dof<=0) return chi; // No info from degenerate/perfect fits
  Vector2 dxy;
  for (auto const & i : elist) {
    if (!isFit(*i)) continue;
    dxy[0] = i->xw;
    dxy[1] = i->yw;
    dxy -= xyMean;
    double cc = dxy.transpose() * i->invCov * dxy;
    maxDeviateSq = MAX(cc/i->expectedTrueChisq , maxDeviateSq);
    chi += cc * i->fitWeight;  // Note returning the fitting chisq
  }
  dofAccum += dof;
  return chi;
}


/////////////////////////////////////////////////////////////////////
// Routines for matches allowing proper motion and parallax
/////////////////////////////////////////////////////////////////////

const double PM_PRIOR = 100. * RESIDUAL_UNIT / WCS_UNIT; /**/

PMMatch::PMMatch(unique_ptr<Detection> e): Match(std::move(e)), pm(0.) {}
				  
void
PMMatch::prepare() const {
  if (isPrepared) return;
  isSolved = false;
  
  // Calculate the total PM Fisher matrix and also
  // the sums over the PMDetection priors
  pmFisher.setZero();
  pmInvFisher.setZero();
  pmTrueCov.setZero();
  priorMean.setZero();
  dof = 0;
  nFit = 0;
  trivialWeights = true;
  
  // First put any prior into Fisher
  pmFisher.diagonal() = priorFisher;

  // Now sum Fisher over all Detections
  PMProjector m;  // (x,y) = m * pm
  for (auto const & i : elist) {
    if (!isFit(*i)) continue;
    m.setZero();
    if (auto ii = dynamic_cast<PMDetection*>(i.get())) {
      // This detection has full PM covariance of its own
      dof += 5;
      nFit++;
      if (ii->fitWeight==1.) {
	//**/cerr << "PMDetection for Fisher:\n" << ii->pmInvCov << endl;
	pmFisher += ii->pmInvCov;
      } else {
	trivialWeights = false;
	pmFisher += ii->pmInvCov * ii->fitWeight;
      }
    } else {
      // Regular single-epoch detection
      m = i->getProjector();
      dof += 2;
      nFit++;
      if (i->fitWeight==1.) {
	pmFisher += m.transpose() * i->invCov *m;
	//**/cerr << "Detection for Fisher:\n" << pmFisher << endl;
      } else {
	trivialWeights=false;
	pmFisher += m.transpose() * i->invCov * (i->fitWeight *m);
      }
    }
  }
  dof -= 5;  // 5 parameter fit here.

  if (dof<0) {
    // Cannot fit a model to these data.
    // Will not be used for anything
    nFit = 0;
    pmInvFisher.setZero();
    pmFisher.setZero();
    pmTrueCov.setZero();
    for (auto const & i : elist)
      i->expectedTrueChisq = 0.;
    isPrepared = true;
    return;
  }

#ifdef USE_EIGEN
  {
    auto llt = pmFisher.llt();
    PMCovariance eye;
    eye.setIdentity();
    pmInvFisher = llt.solve(eye);
  }
#else
  pmInvFisher = pmFisher.inverse(); // ??? More stable?
#endif
  // Go through again and accumulate the PMDetection
  // contributions to chisq that have pmMean in them
  {
    PMSolution tmp;  // work vector
    for (auto const & i : elist) {
      if (!isFit(*i)) continue;
      m.setZero();
      if (auto ii = dynamic_cast<const PMDetection*>(i.get())) {
	tmp = ii->pmInvCov * ii->pmMean;
	tmp *= ii->fitWeight;
	priorMean += tmp;
      }
    }

    priorMean = pmInvFisher * priorMean;    // linear term of chisq
  }
  
  // Calculate the covariance matrix for fitted PM
  if (trivialWeights) {
    pmTrueCov = pmInvFisher;
  } else {
    PMCovariance tmp(0.);
    for (auto const & i : elist) {
      if (!isFit(*i)) continue;
      // Follow notes of 28 Apr 2019
      if (auto ii = dynamic_cast<const PMDetection*>(i.get())) {
	// This detection has full PM covariance of its own
	tmp += ii->pmInvCov * (ii->fitWeight*ii->fitWeight);
      } else {
	m = i->getProjector();
	tmp += m.transpose() * i->invCov * ((i->fitWeight*i->fitWeight) *m);
      }
    }
    pmTrueCov = pmInvFisher * tmp * pmInvFisher;
  }

  // Now get expected chisq for each detection
  // From notes of 28 Apr 2019
  // For a PMDetection,
  // <chisq_i> = 5 - 2*w_i * Tr( pmInvCov_i * pmInvFisher)
  //             + Tr(pmInvCov_i * pmTrueCov
  // for a PMDetection, where pmInvFisher is from the fitter process (with weights),
  // and pmTrueCov is the expected covariance of the fitted PM.
  // For a regular Detection with projection matrix m, the answer is
  // <chisq_i> = 2 - 2*w_i Tr( m^T * invCov_i * m * pmInvFisher)
  //             + Tr(m^T * invCov_i * m * pmTrueCov
  
  Matrix22 covMean = centroidF.inverse();
  for (auto const & i : elist) {
    if  (auto ii = dynamic_cast<PMDetection*>(i.get())) {
      ii->expectedTrueChisq = 5. - 2*ii->fitWeight*traceATB(ii->pmInvCov, pmInvFisher)
	+ traceATB(ii->pmInvCov, pmTrueCov);
    } else {
      //Plain old Detection needs projection matrix
      if (i->invCov(0,0)<=0.) {
	// No error on this detection so no chisq
	i->expectedTrueChisq = 0.;
      } else {
	PMCovariance tmp = i->getProjector().transpose() * i->invCov * i->getProjector();
	if (i->fitWeight==0.) {
	  i->expectedTrueChisq = 2. + traceATB(tmp,pmTrueCov);
	} else if (trivialWeights) {
	  // Unity weight
	  i->expectedTrueChisq = 2. - traceATB(tmp,pmTrueCov);
	} else {
	  // Allow non-unity weight
	  i->expectedTrueChisq = 2. - 2*i->fitWeight*traceATB(tmp,pmInvFisher)
	    + traceATB(tmp,pmTrueCov);
	}
      }
    } 
  }

  isPrepared = true;
}

void
PMMatch::solve() const {
  prepare();
  remap(false);  // Remap the fitted Detections if needed
  if (isSolved) return;

  if (dof<0) {
    // Insufficient data for a solution
    pm.setZero();
    isSolved = true;
    return;
  }
  
  PMProjector m(0.);  // (x,y) = m * pm
  PMSolution beta(0.);
  Vector2 xy;
  // chisq = pm^T * pmFisher * pm - 2 beta^T * pm + const
  // so pm = pmFisher^-1 * beta
  // Note that our priorMean vector is in fact the contribution
  // to this from the PMDetection priors, so we only
  // accumulate beta from non-PM Detections
  for (auto const & i : elist) {
    if (!isFit(*i)) continue;
    if (dynamic_cast<const PMDetection*>(i.get())) {
      continue;  // already cached
    } else {
      // Regular single-epoch detection
      m = i->getProjector();
      xy[0] = i->xw;
      xy[1] = i->yw;
      beta += i->fitWeight * m.transpose() * (i->invCov * xy);
    }
  }
  // Solve the system now
  pm = pmInvFisher * beta + priorMean;
  isSolved = true;
}

double
PMMatch::chisq(int& dofAccum, double& maxDeviateSq,
	       bool dump) const {
  solve();  // Get best PM for current parameters
  double chi=0.;
  if (dof<=0) return chi; // No info from degenerate/perfect fits
  PMProjector m;  // (x,y) = m * pm
  Vector2 dxy;
  PMSolution dpm;
  double cc;
  int j=0;
  for (auto const & i : elist) {
    if (dump) cerr << " Det " << j << " weight " << i->fitWeight << endl;
    j++;
    if (!isFit(*i)) continue;
    if (auto ii = dynamic_cast<const PMDetection*>(i.get())) {
      dpm = ii->pmMean - pm;
      cc = dpm.transpose() * ii->pmInvCov * dpm;
      if (dump) cerr << "   PM cc " << cc << endl;
    } else {    
      // Regular single-epoch detection
      dxy[0] = i->xw;
      dxy[1] = i->yw;
      m = i->getProjector();
      dxy -= m * pm;
      cc = dxy.transpose() * i->invCov * dxy;
      if (dump) cerr << "      cc " << cc 
		     << " dxy " << dxy[0]*WCS_UNIT/RESIDUAL_UNIT
		     << " " << dxy[1]*WCS_UNIT/RESIDUAL_UNIT
		     << endl;
      if (dump) cerr << "m:  " << m << endl;
    }
    chi += cc * i->fitWeight;
    maxDeviateSq = MAX(cc/i->expectedTrueChisq , maxDeviateSq);
  }
  // Add parallax prior
  chi += pm.cwiseProduct(pm).dot(priorFisher);

  dofAccum += dof;
  return chi;
}

bool
PMMatch::sigmaClip(double sigThresh,
		   bool deleteDetection) {
  solve();   // Get PM solution
  if (dof<=0) return false; // Already unusable for fitting

  // Only clip the worst outlier at most
  double maxSq=0.;
  iterator worst = elist.end();
  PMSolution dpm;
  double devSq;
  Vector2 err;
  for (auto i = elist.begin(); i != elist.end(); ++i) {
    if (!isFit(**i)) continue;
    if (auto ii = dynamic_cast<const PMDetection*>(i->get())) {
      dpm = ii->pmMean - pm;
      devSq = dpm.transpose() * ii->pmInvCov * dpm;
    } else {
      err[0] =  (*i)->xw;
      err[1] =  (*i)->yw;
      err -= predict(i->get());
      devSq = err.transpose() * (*i)->invCov * err;
    }
    devSq /= (*i)->expectedTrueChisq;
    if ( devSq > sigThresh*sigThresh && devSq > maxSq) {
      // Mark this as point to clip
      worst = i;
      maxSq = devSq;
    }
  }

  if (worst != elist.end()) {
#ifdef DEBUG
    if (auto ii = dynamic_cast<const PMDetection*>((*worst)->get())) {
      cerr << "clipped PM " << (*worst)->catalogNumber
	   << " / " << (*worst)->objectNumber 
	   << " at " << (*worst)->xw << "," << (*worst)->yw 
	   << " resid/sig " << sqrt(maxSq) 
	   << endl;
    } else {
      cerr << "clipped " << (*worst)->catalogNumber
	   << " / " << (*worst)->objectNumber 
	   << " at " << (*worst)->xw << "," << (*worst)->yw 
	   << " resid " << setw(6) << ((*worst)->xw-predict(worst->get())[0]) * WCS_UNIT/RESIDUAL_UNIT
	   << " " << setw(6) << ((*worst)->yw-predict(worst)[1]) * WCS_UNIT/RESIDUAL_UNIT
	   << " resid/sig " << sqrt(maxSq) 
	   << endl;
    }
#endif
    if (deleteDetection) {
      elist.erase(worst);
    } else {
      (*worst)->clip();
    }
    isPrepared = false;
    isSolved = false;
    return true;
  } else {
    // Nothing to clip
    return false;
  }
}

Vector2
PMMatch::predict(const Detection * d) const {
  if (!d)
    throw AstrometryError("PMMatch::predict called without Detection");
  solve();   // Update full PM solution
  if (dof<0) {
    // No solution possible, return zeros
    return Vector2(0.);
  }
  Vector2 out;
  if (auto dd = dynamic_cast<const PMDetection*> (d)) {
    // Report plain old coords for PMDetection
    out[0] = d->xw;
    out[1] = d->yw;
  } else {
    out = d->getProjector() * pm;
  }
  return out;
}

Matrix22
PMMatch::predictFisher(const Detection* d)  const {
  if (!d)
    throw AstrometryError("PMMatch::predict called without Detection");
  prepare();
  if (dof<0) {
    // Will not get any solution
    return Matrix22(0.);
  }
  PMProjector m = d->getProjector();
  Matrix22 out = (m * pmTrueCov * m.transpose()).inverse();
  return out;
}

int
PMMatch::accumulateChisq(double& chisq,
			 DVector& beta,
			 SymmetricUpdater& updater,
			 bool reuseAlpha) {
  int nP = beta.size();

  prepare();

  if (dof<=0) return 0; // No contribution for degenerate fits

  // First loop updates mapping and accumulates derivs.
  // Save derivatives here
  vector<DMatrix*> dXYdP(elist.size(), nullptr);

  // And accumulate a quantity related to deriv of PM
  // /sum_i (projector_i)^T * invCov_i * dXY_i/dp
  DMatrix dPMdP(5,nP,0.);

  // Projector array to use for all points
  PMProjector m;  // (x,y) = m * pm
  // And a coordinate holder
  Vector2 xy;

  // Acquire derivatives for all fitted points
  // (and update mapping for no-free-param cat if not done yet)
  int ipt=0;
  for (auto i = elist.begin(); i!=elist.end(); ++i, ++ipt) {
    if (!isFit(**i)) continue;
    if (dynamic_cast<const PMDetection*>(i->get())) continue; // No free params in PMDetections

    // Regular single-epoch detection
    int npi = (*i)->map->nParams();
    double xw, yw;
    if (npi>0) {
      dXYdP[ipt] = new DMatrix(2,npi);
      (*i)->map->toWorldDerivs((*i)->xpix, (*i)->ypix,
			       (*i)->xw, (*i)->yw,
			       *dXYdP[ipt],
			       (*i)->color);
    } else if (!isMappedFit) {
      // The world coordinates only need
      // updating if not already mapped.
      // Remap but no derivs needed
      (*i)->map->toWorld((*i)->xpix, (*i)->ypix, 
			 (*i)->xw, (*i)->yw,
			 (*i)->color);
    }
  }
  isMappedFit = true;  // We have just remapped everything
  // Solve the system now
  solve();
  
  // Next loop through detections will
  // accumulate chisq, derivs of PM,
  // and chisq derivs that depend on
  // only one detection.

  Vector2 xyErr;
  Matrix22 invCovW;
  PMSolution pmErr;
  PMProjector cInvM;
  
  map<int, iRange> mapsTouched;
  ipt = 0;
  for (auto i = elist.begin(); i!=elist.end(); ++i, ++ipt) {
    if (!isFit(**i)) continue;

    if (auto ii = dynamic_cast<const PMDetection*> (i->get())) {
      // PMDetection has no parameter dependence, just chisq
      pmErr = ii->pmMean - pm;
      double tmp = (pmErr.transpose() * ii->pmInvCov * pmErr);
      chisq += ii->fitWeight*tmp;
    } else {
      // First calculate chisq contribution
      invCovW = (*i)->invCov * (*i)->fitWeight;
      xyErr[0] = (*i)->xw;
      xyErr[1] = (*i)->yw;
      m = (*i)->getProjector();
      xyErr -= m * pm;
      chisq += xyErr.transpose() * invCovW * xyErr; 

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
	DMatrix cInvD = invCovW * dxy1;
	beta.subVector(ip, ip+np) -= cInvD.transpose() * xyErr;

	// Contribution to derivatives of PM:
	dPMdP.subMatrix(0, 5, ip, ip+np) += m.transpose() * cInvD;

	if (!reuseAlpha) {
	  // Increment the alpha matrix
	  // First the blocks astride the diagonal:
	  updater.update(mapNumber, ip, dxy1, invCovW);

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
			   invCovW,
			   mapNumber,  ip,  dxy1);
	    istart2+=np2;
	  }
	} // End inner loop of dXYdP^T * invC * dXYdP for this segment
	istart+=np;
      } // outer parameter segment loop
      if (dXYdP[ipt]) {delete dXYdP[ipt]; dXYdP[ipt]=nullptr;}
    } // endif for PMDetections
  } // object loop

  // Add parallax prior
  chisq += pm.cwiseProduct(pm).dot(priorFisher);

  if (!reuseAlpha) {
    // Subtract terms with derivatives of PM,
    /* but without touching the entire alpha matrix every time:
       alpha -=  dPMdP^T * pmInvFisher * dPMdP);
    */

    // Do updates parameter block by parameter block
    PMCovariance negCov = -pmInvFisher;
    for (auto m1 = mapsTouched.begin();
	 m1!=mapsTouched.end();
	 ++m1) {
      int map1 = m1->first;
      int i1 = m1->second.startIndex;
      int n1 = m1->second.nParams;
      if (n1<=0) continue;
      DMatrix dxy1 = dPMdP.subMatrix(0,5, i1, i1+n1);
      if (reuseAlpha)
	continue;

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
    } // Finished putting terms from mean into alpha
  }

  return dof;
}



/////////////////////////////////////////////////////////////////////
// Parameter adjustment operations
/////////////////////////////////////////////////////////////////////
void
CoordAlign::setParams(const DVector& p) {
  // First, change parameters in the map collection
  pmc.setParams(p);
  // Then alert all Matches that their world coordinates are crap
  for (auto m : mlist)
    m->mapsHaveChanged();
}

void
CoordAlign::operator()(const DVector& p, double& chisq,
		       DVector& beta, DMatrix& alpha,
		       bool reuseAlpha) {
  int nP = pmc.nParams();
  Assert(p.size()==nP);
  Assert(beta.size()==nP);
  Assert(alpha.rows()==nP);
  setParams(p);  // Note this will force remapping
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
    typedef Eigen::LLT<Eigen::Ref<typename DMatrix::Base>,Eigen::Lower> InplaceLLT;
    typedef Eigen::LLT<typename DMatrix::Base,Eigen::Lower> LLT;
    InplaceLLT* inplaceLLT = nullptr;
    LLT* llt = nullptr;
    if (inPlace) {
      inplaceLLT = new InplaceLLT(alpha);
      if (inplaceLLT->info()==Eigen::NumericalIssue)
	choleskyFails = true;
    } else {
      llt = new LLT(alpha);
      if (llt->info()==Eigen::NumericalIssue)
	choleskyFails = true;
    }
#endif
    
    // If the Cholesky decomposition failed for non-pos-def matrix, then
    // as a diagnostic we will do an SVD and report the nature of degeneracies.
    // Then exit with an error.
    // And since alpha is symmetric, we will use eigenvector routines to do the
    // SVD, which is what we want anyway to find the non-pos-def values

    if (choleskyFails) {
      cerr << "Caught exception during Cholesky" << endl;
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
      if (inPlace)
	beta = inplaceLLT->solve(beta);
      else
	beta = llt->solve(beta);
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
      int dof;
      double maxDev;
      double newChisq = chisqDOF(dof, maxDev);
      timer.stop();
      if (reportToCerr) {
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
#ifdef USE_EIGEN
	if (llt) delete llt;
	if (inplaceLLT) delete inplaceLLT;
#endif
	return newChisq;
      }
      // Want another Newton iteration, but keep alpha as before
      p = newP;
      (*this)(p, oldChisq, beta, alpha, true);
    }
    // If we reach this place, Newton is going backwards or nowhere, slowly.
    // So just give it up.
#ifdef USE_EIGEN
    if (llt) delete llt;
    if (inplaceLLT) delete inplaceLLT;
#endif
  }

  // ??? Signal that alpha should be fixed for all iterations of Marquardt?
  Marquardt<CoordAlign> marq(*this);
  marq.setRelTolerance(relativeTolerance);
  marq.setSaveMemory();  
  double chisq = marq.fit(p, DefaultMaxIterations, reportToCerr);
  setParams(p);
  {
    int dof;
    double maxDev;
    chisq = chisqDOF(dof, maxDev);
  }
  return chisq;
}

void
CoordAlign::remap(bool doAll) {
  for (auto i : mlist)
    i->remap(doAll);
}

int
CoordAlign::sigmaClip(double sigThresh, bool doReserved, bool clipEntireMatch,
		      bool logging) {
  int nclip=0;

  Stopwatch timer;
  timer.start();

#ifdef _OPENMP
  const int chunk=100;
  // Make a vector version on mlist since openMP requires
  // integers in a for loop.
  int n = mlist.size();
  vector<Match*> mvec(n);
  std::copy(mlist.begin(), mlist.end(), mvec.begin());
#pragma omp parallel for schedule(dynamic,chunk) reduction(+:nclip)
  for (auto ii=0; ii<n; ++ii) {
    auto i = mvec[ii];
#else
    for (auto i : mlist) {
#endif
    // Skip this one if it's reserved and doReserved=false,
    // or vice-versa
    if (doReserved ^ i->getReserved()) continue;
  if ( i->sigmaClip(sigThresh)) {
      nclip++;
      if (clipEntireMatch) i->clipAll();
    }
  }
  timer.stop();
  if (logging)
    cerr << "-->Sigma clipping done in " << timer << " sec" << endl;
  return nclip;
}

double
CoordAlign::chisqDOF(int& dof, double& maxDeviate, 
		     bool doReserved) const {
  dof=0;
  maxDeviate=0.;
  double chisq=0.;

#ifdef _OPENMP
  const int chunk=100;
  // Make a vector version of mlist since openMP requires
  // integers in a for loop.
  int n = mlist.size();
  vector<Match*> mvec(n);
  std::copy(mlist.begin(), mlist.end(), mvec.begin());
#pragma omp parallel for schedule(dynamic,chunk) reduction(+:dof) reduction(+:chisq)
  for (auto ii=0; ii<n; ++ii) {
    auto i = mvec[ii];
#else
  for (auto i : mlist) {
#endif
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
    if ( (i->getReserved() ^ doReserved) ||
	 i->getDOF() < 0 || i->fitSize() < minMatches) {
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
    if ((i->getReserved() ^ doReserved) ||
	i->getDOF() < 0 || i->fitSize() < minMatches) continue;

    int ddcount=0;
    for(auto const & d : *i) {
      if(!(d->isClipped) && d->catalogNumber==catalog)
	ddcount++;
    }

    if(ddcount>0) {
      mcount++;
      dcount+=ddcount;
    }
  }
}
