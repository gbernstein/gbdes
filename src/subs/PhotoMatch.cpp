// Photometric matching and fitting classes.

#include "PhotoMatch.h"
#include <list>
using std::list;
#include <set>
using std::set;
#include "StringStuff.h"
#include <fstream>
#include "AstronomicalConstants.h"
#include <algorithm>
#include "Stopwatch.h"
#include "Units.h"

// Need Eigenvalue routines
#ifdef USE_EIGEN
#include "Eigen/Eigenvalues"
#endif

#ifdef _OPENMP
#include <omp.h>
#endif

// #define DEBUG
#include "Marquardt.h"

using namespace photometry;

double
Detection::residMag() const {
  Assert(itsMatch);
  itsMatch->solve();
  return (magOut - itsMatch->getMean()) / MMAG;
}

double
Detection::trueChisq() const {
  Assert(itsMatch);  
  itsMatch->solve();
  double resid = magOut - itsMatch->getMean();
  return resid * resid * invVar;
}

bool
Match::isFit(const Detection & e) {
  // A Detection will contribute to fit if it has nonzero weight and is not clipped.
  return e.fitWeight>0.;
}

Match::Match(unique_ptr<Detection> e): elist(), nFit(0), dof(0),
			    isReserved(false), isPrepared(false),
			    isMappedFit(false), isMappedAll(false),
			    isSolved(false), trivialWeights(true) {
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
Match::remove(const Detection & e) {
  isPrepared = false;
  isSolved = false;
  elist.remove_if([&e](unique_ptr<Detection> const & d) { return &e == d.get(); });
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
  for (auto const & i : elist)
    if (i)
      i->clip();
  nFit = 0;
}

void
Match::prepare() const {
  // Recalculate, if needed, quantities that do not depend on the
  // detection position values - namely the (inverse) variance of
  // the centroid.
  if (isPrepared) return;
  isSolved = false;   // Change here means centroid solution is invalid
  meanF = 0.;
  trueMeanVar = 0.;
  dof = 0;
  nFit = 0;
  trivialWeights = true;
  for (auto const & i : elist) {
    if (!isFit(*i)) continue;
    dof += 1;
    nFit++;
    if (i->fitWeight != 1.) {
      trivialWeights = false;  // Need explicit weights in formulae
      meanF += i->invVar * i->fitWeight;
    } else {
      meanF += i->invVar;
    }
      
  }
  dof -= 1;  // Remove 1 DOF for freedom in mean

  if (dof<0) {
    // Degenerate fit, will never yield a solution of interest
    for (auto const & i : elist)
      i->expectedTrueChisq = 0.;
    meanF = 0.;
    trueMeanVar = 0.;
    isPrepared = true;
    nFit = 0;
    return;
  }
	
  double invF = 1./meanF;
  // Calculate expected centroid 
  if (trivialWeights) {
    trueMeanVar = invF;
  } else {
    // Altered expectation for cov because of weighting
    // From notes of 28 Apr 2019,
    double tmp = 0.;
    for (auto const & i : elist) {
      if (!isFit(*i)) continue;
      tmp += i->invVar * (i->fitWeight*i->fitWeight);
    }
    trueMeanVar = invF * tmp * invF;
  }

  // Now get expected chisq for each detection, accounting for weighting
  // From notes of 28 Apr 2019,
  // where invF is the inverse of centroidF
  for (auto const & i : elist) {
    if (i->invVar<=0.) {
      // No error on this detection so no chisq
      i->expectedTrueChisq = 0.;
    } else if (i->fitWeight==0.) {
      // Has uncertainties but was not included in the fit
      i->expectedTrueChisq = 1. + i->invVar * trueMeanVar;
    } else if (trivialWeights) {
      // included in fit, unit weight gauranteed
      i->expectedTrueChisq = 1. - i->invVar * trueMeanVar;
    } else {
      // Include non-trivial weight
      i->expectedTrueChisq = 1 - 2.*i->fitWeight*i->invVar*invF
	+ i->invVar*trueMeanVar;
    }
  }
  // Calculate expected chisq for each Detection:
  isPrepared = true;
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
	i->magOut = i->map->forward(i->magIn, i->args);
	isSolved = false;  // Solution has changed
      }
    } else if (doAll) {
      // Only remap non-fits if requested
      i->magOut = i->map->forward(i->magIn, i->args);
    }
  }
  // Update state
  isMappedFit = true;
  if (doAll) isMappedAll = true;
}

void
Match::solve() const {
  // Recalculate the mean, if needed.
  prepare();     // Make sure Fisher matrices are current
  remap(false);  // Make sure world coords of fitted Detections are current
  if (isSolved) return;
  isSolved = true;
  if (dof<0) {
    // Do not execute degenerate fit
    meanMag = 0.;
    return;
  }
  double sumwm = 0.;
  for (auto const & i : elist) {
    if (!isFit(*i)) continue;
    sumwm += i->invVar * i->magOut * i->fitWeight;
  }
  meanMag = meanF * sumwm;
}

///////////////////////////////////////////////////////////
// Coordinate-matching routines
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
  
  // No contributions to fit for <2 detections:
  if (dof<=0) return 0;

  // Update mapping and save derivatives for each detection:
  vector<DVector*> di(elist.size());
  int ipt=0;
  for (auto i = elist.begin(); i!=elist.end(); ++i, ++ipt) {
    if (!isFit(**i)) continue;
    int npi = (*i)->map->nParams();
    if (npi>0) {
      di[ipt] = new DVector(npi);
      (*i)->magOut = (*i)->map->forwardDerivs((*i)->magIn, (*i)->args,
					      *di[ipt]);
    } else if (!isMappedFit) {
      // If there are no free parameters, we only
      // need to remap if the maps have changed.
      (*i)->magOut = (*i)->map->forward((*i)->magIn, (*i)->args);
    }
  }
  isMappedFit = true;  // This is now true if it wasn't before.

  // Get new mean mag
  solve();
  
  DVector dMean(nP, 0.); // Derive of mean wrt parameters
  // Places to put mag resid and invVar * weight for each Detection
  double resid;
  double invVarW;  

  map<int, iRange> mapsTouched;
  ipt = 0;
  for (auto i = elist.begin(); i!=elist.end(); ++i, ++ipt) {
    if (!isFit(**i)) continue;
    invVarW = (*i)->invVar * (*i)->fitWeight;
    resid = (*i)->magOut - meanMag;
    chisq += resid*resid*invVarW;

    // Accumulate derivatives:
    int istart=0;
    for (int iMap=0; iMap<(*i)->map->nMaps(); iMap++) {
      int ip=(*i)->map->startIndex(iMap);
      int np=(*i)->map->nSubParams(iMap);
      if (np==0) continue;
      int mapNumber = (*i)->map->mapNumber(iMap);
      // Keep track of parameter ranges we've messed with:
      mapsTouched[mapNumber] = iRange(ip,np);
      DVector dm1=di[ipt]->subVector(istart,istart+np);
      beta.subVector(ip, ip+np) -= (invVarW * resid)*dm1;

      // Derivatives of the mean position:
      dMean.subVector(ip, ip+np) += invVarW*dm1;

      if (!reuseAlpha) {
	// Increment the alpha matrix
	// First the blocks astride the diagonal:
	updater.rankOneUpdate(mapNumber, ip, dm1, invVarW);

	// Put the wti factor into dm:
	int istart2 = istart+np;
	for (int iMap2=iMap+1; iMap2<(*i)->map->nMaps(); iMap2++) {
	  int ip2=(*i)->map->startIndex(iMap2);
	  int np2=(*i)->map->nSubParams(iMap2);
	  int mapNumber2 = (*i)->map->mapNumber(iMap2);
	  if (np2==0) continue;
	  DVector dm2=di[ipt]->subVector(istart2,istart2+np2);
	  // Update below the diagonal:
	  updater.rankOneUpdate(mapNumber2, ip2, dm2,
				mapNumber,  ip,  dm1,
				invVarW);
	  istart2+=np2;
	}
	istart+=np;
      }
    } // outer parameter segment loop
    if (di[ipt]) {delete di[ipt]; di[ipt]=nullptr;}
  } // object loop


  if (!reuseAlpha) {
    // Subtract effects of derivatives on mean
    /*  We want to do this, but without touching the entire alpha matrix every time:
	alpha -=  (dMean ^ dMean)/wt;
    */

    // Do updates parameter block by parameter block
    double negVarMean = -1./meanF;
    for (auto m1=mapsTouched.begin();
	 m1 != mapsTouched.end();
	 ++m1) {
      int map1 = m1->first;
      int i1 = m1->second.startIndex;
      int n1 = m1->second.nParams;
      if (n1<=0) continue;
      DVector dm1 = dMean.subVector(i1, i1+n1);
      // Update astride diagonal:
      updater.rankOneUpdate(map1, i1, dm1, negVarMean);

      auto m2 = m1;
      ++m2;
      for (; m2 != mapsTouched.end(); ++m2) {
	int map2 = m2->first;
	int i2 = m2->second.startIndex;
	int n2 = m2->second.nParams;
	if (n2<=0) continue;
	DVector dm2 = dMean.subVector(i2, i2+n2);
	updater.rankOneUpdate(map2, i2, dm2,
			      map1, i1, dm1,
			      negVarMean);
      }
    }
  }
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

  for (auto i = elist.begin(); i != elist.end(); ++i) {
    if (!isFit(**i)) continue;
    double resid = (*i)->magOut - meanMag;
    double devSq = resid*resid*(*i)->invVar;
    // Calculate deviation relative to expected chisq
    devSq /= (*i)->expectedTrueChisq;
    
    if ( devSq > sigThresh*sigThresh && devSq > maxSq) {
      // Mark this as point to clip
      worst = i;
      maxSq = devSq;
    }
  }

  if (worst != elist.end()) {
#ifdef DEBUG
    cerr << "clipped " << worst->catalogNumber
	 << " / " << worst->objectNumber 
	 << " at " << worst->magOut
	 << " resid " << setw(6) << (worst->magOut-meanMag)
	 << " sigma " << 1./sqrt(worst->invVar)
	 << endl;
#endif
      if (deleteDetection) {
	elist.erase(worst);
      }	else {
	(*worst)->isClipped = true;
      }
      isPrepared = false;
      isSolved = false;  // Note that magOut's are still valid
      return true;
  } else {
    // Nothing to clip
    return false;
  }
}

double
Match::chisq(int& dofAccum, double& maxDeviateSq) const {
  solve();
  double chi=0.;
  if (dof<=0) return chi;
  for (auto const & i : elist) {
    if (!isFit(*i)) continue;
    double resid = i->magOut - meanMag;
    double cc = resid * resid * i->invVar;
    maxDeviateSq = MAX(cc/i->expectedTrueChisq, maxDeviateSq);
    chi += cc * i->fitWeight;  // Returning the fitting chisq here
  }
  dofAccum += dof;
  return chi;
}

/////////////////////////////////////////////////////////////////////
// Parameter adjustment operations
/////////////////////////////////////////////////////////////////////

void
PhotoAlign::setParams(const DVector& p) {
  // First, change parameters in the map collection
  pmc.setParams(p.subVector(0,pmc.nParams()));
  // Then alert all Matches that their world coordinates are crap
  for (auto const & m : mlist)
    m->mapsHaveChanged();
  int startIndex = pmc.nParams();
  for (auto i : priors) {
    i->mapsHaveChanged();
    if (i->isDegenerate()) continue;
    i->setParams(p.subVector(startIndex, startIndex + i->nParams()));
    startIndex += i->nParams();
  }
}

DVector
PhotoAlign::getParams() const {
  DVector p(nParams(), -888.);
  p.subVector(0,pmc.nParams()) = pmc.getParams();
  int startIndex = pmc.nParams();
  for (auto i : priors) {
    if (i->isDegenerate()) continue;
    p.subVector(startIndex, startIndex + i->nParams()) = i->getParams();
    startIndex += i->nParams();
  }
  return p;
}

void
PhotoAlign::operator()(const DVector& p, double& chisq,
		       DVector& beta, DMatrix& alpha,
		       bool reuseAlpha) {
  countPriorParams();
  int nP = nParams();
  Assert(p.size()==nP);
  Assert(beta.size()==nP);
  Assert(alpha.rows()==nP);
  setParams(p);
  double newChisq=0.;
  beta.setZero();
  if (!reuseAlpha) alpha.setZero();
  int matchCtr=0;

  const int NumberOfLocks = 2000;
  SymmetricUpdater updater(alpha, maxMapNumber, NumberOfLocks);

#ifdef _OPENMP
  const int chunk=200;
  vector<Match*> vi(mlist.size());

#pragma omp parallel reduction(+:newChisq) 
  {
#pragma omp single
    {
      // Distribute the matches in chunks that are scattered
      // around the array to discourage parameter collisions
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
	  iThread++;
	  if (iThread >= nThreads) {
	    iRound++;
	    iThread=0;
	  }
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
  for (auto const & i : mlist) {
    Match* m = i.get();
    if (matchCtr%10000==0) cerr << "# accumulating chisq at match # " 
				<< matchCtr  //**<< " newChisq " << newChisq
				<< endl;
    matchCtr++;
    if ( m->getReserved() ) continue;	//skip reserved objects
    m->accumulateChisq(newChisq, beta, updater, reuseAlpha);
  }
#endif
  chisq = newChisq;

  // Now need to accumulate the contributions from the priors.
  // A single thread will do.
  for (auto i : priors)
    i->accumulateChisq(chisq, beta, updater, reuseAlpha);

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
	string badAtom="";
	bool badIsMap; // Is the bad parameter in a map or in a prior?
	if (i < pmc.nParams()) {
	  badAtom = pmc.atomHavingParameter(i);
	  badIsMap = true;
	} else {
	  // Look among the priors for this parameter
	  for (auto iprior : priors) {
	    if ( i >= iprior->startIndex() &&
		 i < iprior->startIndex() + iprior->nParams()) {
	      badAtom = iprior->getName();
	      badIsMap = false;
	      break;
	    }
	  }
	}
	if (badAtom.empty()) {
	  FormatAndThrow<PhotometryError>() << "Could not locate parent map for "
					    << " degenerate parameter " << i;
	}
	// Add to (or make) a list of the frozen parameters in this atom
	// Is it a newly frozen parameter?
	if (!frozenMaps.count(badAtom) || !frozenMaps[badAtom].count(i))
	  newlyFrozenMaps.insert(badAtom);
	// Add to (or make) a list of the frozen parameters in this atom
	frozenMaps[badAtom].insert(i);
	// Fudge matrix to freeze parameter:
	alpha(i,i) = 1.;
	beta[i] = 0.;
      } else {
	// Something is weird if a frozen parameter is now constrained
	if (frozenParameters.count(i)>0) {
	  if (i < pmc.nParams()) {
	    string badAtom = pmc.atomHavingParameter(i);
	    FormatAndThrow<PhotometryError>() << "Frozen parameter " << i
					      << " in map " << badAtom
					      << " became constrained??";
	  } else {
	    FormatAndThrow<PhotometryError>() << "Frozen parameter " << i
					      << " in prior became constrained??";
	  }
	}
      }
    } // End alpha row loop

    for (auto badAtom : newlyFrozenMaps) {
      // Print message about freezing parameters
      int startIndex, nParams;
      if (pmc.mapExists(badAtom)) {
	// Message for a map parameter:
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
      } else {
	// Message for a prior parameter
	cerr << "Freezing " << frozenMaps[badAtom].size()
	     << " parameters in map " << badAtom
	     << endl;  // ??? Could get more specific here.
      }
    }
  } // End degenerate parameter check
}

double
PhotoAlign::fitOnce(bool reportToCerr, bool inPlace) {
  DVector p = getParams();
  // First will try doing Newton iterations, keeping a fixed Hessian.
  // If it increases chisq or takes too long to converge, we will 
  // go into the Marquardt solver.

  {
    // Get chisq, beta, alpha at starting parameters
    double oldChisq = 0.;
    int nP = p.size();
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
	    if (j < pmc.nParams()) {
	      string badAtom = pmc.atomHavingParameter(j);
	      int startIndex, nParams;
	      pmc.parameterIndicesOf(badAtom, startIndex, nParams);
	      cerr << "Coefficient " << U(j, isv) 
		   << " at parameter " << j 
		   << " Map " << badAtom 
		   << " " << j - startIndex << " of " << nParams
		   << endl;
	    } else {
	      cerr << "Coefficient " << U(j, isv) 
		   << " at parameter " << j 
		   << " in priors"
		   << endl;
	  }
	}
	exit(1);
      }
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
  
  Marquardt<PhotoAlign> marq(*this);
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

/**
void
PhotoAlign::remap(bool doAll) const {
  for (auto i : mlist)
    i->remap(doAll);

  for (auto i : priors) 
    i->remap();
}
**/

int
PhotoAlign::sigmaClip(double sigThresh, bool doReserved, bool clipEntireMatch,
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
    for (auto const & i : mlist) {
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

int
PhotoAlign::sigmaClipPrior(double sigThresh, bool clipEntirePrior) {
  int nClip = 0;
  for (auto i : priors) 
    if (i->sigmaClip(sigThresh)) {
      // Clipped something.
      nClip++;
      if (clipEntirePrior) i->clipAll();
    }
  return nClip;
}

double
PhotoAlign::chisqDOF(int& dof, double& maxDeviate, 
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
  for (auto const & i : mlist) {
#endif
    if (doReserved ^ i->getReserved()) continue;
    chisq += i->chisq(dof, maxDeviate);
  }
  maxDeviate = sqrt(maxDeviate);
  if (!doReserved) {
    // If doing the fitted objects, include prior and adjust DOF for fit
    dof -= pmc.nParams();
    for (auto i : priors) 
      chisq += i->chisq(dof);
  }
  return chisq;
}

void 
PhotoAlign::count(long int& mcount, long int& dcount, 
		  bool doReserved, int minMatches) const {
  mcount = 0;
  dcount = 0;
  for (auto const & i: mlist) {
    if ((i->getReserved() ^ doReserved) 
	|| i->getDOF()<0 || i->fitSize() < minMatches) continue;
    mcount++;
    dcount+=i->fitSize();
  }
}
 
void
PhotoAlign::count(long int& mcount, long int& dcount,
                  bool doReserved, int minMatches, long catalog) const {
  mcount = 0;
  dcount = 0;
  for (auto const & i : mlist) {
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

void
PhotoAlign::countPriorParams() {
  // Reassign all counts/pointers for parameters of priors.
  // Degenerate priors will be skipped.
  int startIndex = pmc.nParams();
  int mapNumber = pmc.nFreeMaps();
  for (auto i : priors) {
    if (i->isDegenerate()) continue;
    i->globalStartIndex = startIndex;
    startIndex += i->nParams();
    i->globalMapNumber = mapNumber++;
  }
  nPriorParams = startIndex - pmc.nParams();
  maxMapNumber = mapNumber;
}

