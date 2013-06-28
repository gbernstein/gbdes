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

#ifdef _OPENMP
#include <omp.h>
#endif

// #define DEBUG
#include "Marquardt.h"

using namespace photometry;


/////////////////////////////////////////////////////////////
///// Class that regulates rank-one updates to small portions of
///// the giant alpha symmetric matrix, allowing finer-grained locking
/////////////////////////////////////////////////////////////
AlphaUpdater::AlphaUpdater(tmv::SymMatrix<double>& alpha_,
			   int nMaps_,
			   int nLocks): alpha(alpha_), nMaps(nMaps_) {
#ifdef _OPENMP
  // We will set up locks for altering subregions of the matrix.
  // The regions will be divided into sets of parameters from
  // a range of PixelMaps.  nMaps is the total number of maps
  // with free parameters and nLocks is the number of distinct lockable regions
  // we want to have.  Matrix is symmetric so we are dividing a triangular region.
  // First decide how many map elements per lock
  int nBlocks = nMaps*(nMaps+1)/2;
  if (nLocks > nBlocks || nLocks<=0) {
    blockLength = 1;
    nLocks = nBlocks;
  } else {
    blockLength = nBlocks / nLocks;
    nLocks = (nBlocks+blockLength+-1) / blockLength;
  }
  locks.resize(nLocks);
  // build array of locks
  for (int i=0; i<locks.size(); i++)
    omp_init_lock(&locks[i]);
#endif
}

AlphaUpdater::~AlphaUpdater() {
#ifdef _OPENMP    // Release locks
  for (int i=0; i<locks.size(); i++)
    omp_destroy_lock(&locks[i]);
#endif
}


// Update submatrix with corner at (startIndex1,startIndex2) with += scalar * v1 ^ v2
// map[12] are sequence numbers of map components in the PixelMapCollection, used
// to divide the matrix into blocks for locking in multithreaded case
void 
AlphaUpdater::rankOneUpdate(int map1, int startIndex1, tmv::ConstVectorView<double> v1, 
			    int map2, int startIndex2, tmv::ConstVectorView<double> v2,
			    double scalar) {
  if (v1.size() <= 0 || v2.size() <=0) return;	// Nothing to add.
  bool diagonal = (map1 == map2);
  bool swapVectors = (map1 < map2);

#ifdef _OPENMP
  // Get the block index corresponding to these 2 keys
  int block = (swapVectors ? (map2*(map2+1)/2 + map1) : (map1*(map1+1)/2 + map2))
    / blockLength;
  Assert(block < locks.size());
  // Set lock for its region - will wait here if busy!
  omp_set_lock(&locks[block]);
  // Pre-multiply the vectors?
#endif
  if (diagonal) {
    alpha.subSymMatrix(startIndex1, startIndex1 + v1.size()) += scalar * v1 ^ v1;
  } else if (swapVectors) {
    alpha.subMatrix(startIndex2, startIndex2+v2.size(),
		    startIndex1, startIndex1+v1.size()) += scalar * v2 ^ v1;
  } else {
    alpha.subMatrix(startIndex1, startIndex1+v1.size(),
		    startIndex2, startIndex2+v2.size()) += scalar * v1 ^ v2;
  }
#ifdef _OPENMP
  omp_unset_lock(&locks[block]);
#endif
}

bool
Match::isFit(const Detection* e) {
  // A Detection will contribute to fit if it has nonzero weight and is not clipped.
  return !(e->isClipped) && e->wt>0.;
}

Match::Match(Detection *e): elist(), nFit(0), isReserved(false) {
  add(e);
}

void
Match::add(Detection *e) {
  elist.push_back(e);
  e->itsMatch = this;
  e->isClipped = false;
  if ( isFit(e)) nFit++;
}

void
Match::remove(Detection* e) {
  elist.remove(e);  
  e->itsMatch = 0;
  if ( isFit(e) ) nFit--;
}

list<Detection*>::iterator 
Match::erase(list<Detection*>::iterator i,
	     bool deleteDetection) {
  if (isFit(*i)) nFit--;
  if (deleteDetection) delete *i;
  return elist.erase(i);
}

void
Match::clear(bool deleteDetections) {
  for (list<Detection*>::iterator i=elist.begin();
       i!=elist.end();
       ++i) {
    if (deleteDetections)
      delete *i;
    else
      (*i)->itsMatch = 0;
  }
  elist.clear();
  nFit = 0;
}

void
Match::clipAll() {
  for (list<Detection*>::iterator i=elist.begin();
       i!=elist.end();
       ++i) 
    (*i)->isClipped = true;
  nFit = 0;
}

// Update and return count of fittable objects:
void
Match::countFit() {
  nFit = 0;
  for (list<Detection*>::iterator i=elist.begin();
       i!=elist.end();
       ++i)
    if (isFit(*i)) nFit++;
}

void
Match::getMean(double& m) const {
  double wt;
  getMean(m,wt);
}

void
Match::getMean(double& m, double& wt) const {
  if (nFit<1) {
    wt=m=0.;
    return;
  }
  double sum_mw=0.;
  double sum_w =0.;
  for (list<Detection*>::const_iterator i=elist.begin();
       i!=elist.end();
       ++i) {
    if ((*i)->isClipped) continue;
    double m = (*i)->magOut;
    double w = (*i)->wt;
    sum_mw += m*w;
    sum_w += w;
  }
  m = sum_mw/sum_w;
  wt = sum_w;
}


///////////////////////////////////////////////////////////
// Coordinate-matching routines
///////////////////////////////////////////////////////////

void
Match::remap() {
  for (list<Detection*>::iterator i=elist.begin(); 
       i!=elist.end(); 
       ++i) 
    (*i)->magOut = (*i)->map->forward((*i)->magIn, (*i)->args);
}

// A structure that describes a range of parameters
struct iRange {
  iRange(int i=0, int n=0): startIndex(i), nParams(n) {}
  int startIndex;
  int nParams;
};

int
Match::accumulateChisq(double& chisq,
		       DVector& beta,
		       AlphaUpdater& updater) {
		       //**		       tmv::SymMatrix<double>& alpha) {
  double mean;
  double wt;
  int nP = beta.size();

  // No contributions to fit for <2 detections:
  if (nFit<=1) return 0;

  // Update mapping and save derivatives for each detection:
  vector<DVector*> di(elist.size());
  int ipt=0;
  for (list<Detection*>::iterator i=elist.begin(); 
       i!=elist.end(); 
       ++i, ipt++) {
    if (!isFit(*i)) continue;
    int npi = (*i)->map->nParams();
    if (npi>0) {
      di[ipt] = new DVector(npi);
      (*i)->magOut = (*i)->map->forwardDerivs((*i)->magIn, (*i)->args,
					      *di[ipt]);
    } else {
      (*i)->magOut = (*i)->map->forward((*i)->magIn, (*i)->args);
    }
  }

  getMean(mean,wt);
  DVector dmean(nP, 0.);

  map<int, iRange> mapsTouched;
  ipt = 0;
  for (list<Detection*>::iterator i=elist.begin(); 
       i!=elist.end(); 
       ++i, ipt++) {
    if (!isFit(*i)) continue;
    double wti=(*i)->wt;
    double mi=(*i)->magOut;

    chisq += (mi-mean)*(mi-mean)*wti;

    // Accumulate derivatives:
    int istart=0;
    for (int iMap=0; iMap<(*i)->map->nMaps(); iMap++) {
      int ip=(*i)->map->startIndex(iMap);
      int np=(*i)->map->nSubParams(iMap);
      if (np==0) continue;
      int mapNumber = (*i)->map->mapNumber(iMap);
      // Keep track of parameter ranges we've messed with:
      mapsTouched[mapNumber] = iRange(ip,np);
      tmv::VectorView<double> dm=di[ipt]->subVector(istart,istart+np);
      beta.subVector(ip, ip+np) -= wti*(mi-mean)*dm;

      // Derivatives of the mean position:
      dmean.subVector(ip, ip+np) += wti*dm;

      // The alpha matrix: a little messy because of symmetry
      //**      alpha.subSymMatrix(ip, ip+np) += wti*(dm ^ dm));
      updater.rankOneUpdate(mapNumber, ip, dm, 
			    mapNumber, ip, dm, wti);

      int istart2 = istart+np;
      for (int iMap2=iMap+1; iMap2<(*i)->map->nMaps(); iMap2++) {
	int ip2=(*i)->map->startIndex(iMap2);
	int np2=(*i)->map->nSubParams(iMap2);
	int mapNumber2 = (*i)->map->mapNumber(iMap2);
	if (np2==0) continue;
	tmv::VectorView<double> dm2=di[ipt]->subVector(istart2,istart2+np2);
	// Note that subMatrix here will not cross diagonal:
	//**	  alpha.subMatrix(ip, ip+np, ip2, ip2+np2) += wti*(dm ^ dm2);
	updater.rankOneUpdate(mapNumber2, ip2, dm2, 
			      mapNumber,  ip,  dm, wti);
	istart2+=np2;
      }
      istart+=np;
    } // outer parameter segment loop
    if (di[ipt]) {delete di[ipt]; di[ipt]=0;}
  } // object loop

  // Subtract effects of derivatives on mean
  /*  We want to do this, but without touching the entire alpha matrix every time:
  alpha -=  (dmean ^ dmean)/wt;
  */

  // Do updates parameter block by parameter block
  for (map<int,iRange>::iterator m1=mapsTouched.begin();
       m1 != mapsTouched.end();
       ++m1) {
    int map1 = m1->first;
    int i1 = m1->second.startIndex;
    int n1 = m1->second.nParams;
    if (n1<=0) continue;
    for (map<int, iRange>::iterator m2=m1;
	 m2 != mapsTouched.end();
	 ++m2) {
      int map2 = m2->first;
      int i2 = m2->second.startIndex;
      int n2 = m2->second.nParams;
      if (n2<=0) continue;
      tmv::VectorView<double> dm1 = dmean.subVector(i1, i1+n1);
      tmv::VectorView<double> dm2 = dmean.subVector(i2, i2+n2);
      updater.rankOneUpdate(map2, i2, dm2,
			    map1, i1, dm1, -1./wt);
    }
  }

  return nFit-1;
}

bool
Match::sigmaClip(double sigThresh,
		 bool deleteDetection) {
  // Only clip the worst outlier at most
  double mean;
  if (nFit<=1) return false;
  double maxSq=0.;
  list<Detection*>::iterator worst=elist.end();
  getMean(mean);
  for (list<Detection*>::iterator i=elist.begin(); 
       i!=elist.end(); ++i) {
    if (!isFit(*i)) continue;
    double mi = (*i)->magOut;
    double clipsq = (*i)->clipsq;
    double devSq = (mi-mean)*(mi-mean)*clipsq;
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
	 << " at " << (*worst)->magOut
	 << " resid " << setw(6) << ((*worst)->magOut-mean)
	 << " sigma " << 1./sqrt((*worst)->clipsqx)
	 << endl;
#endif
      if (deleteDetection) {
	delete *worst;
	elist.erase(worst);
      }	else {
	(*worst)->isClipped = true;
      }
      nFit--;
      return true;
  } else {
    // Nothing to clip
    return false;
  }
}

double
Match::chisq(int& dof, double& maxDeviateSq) const {
  double mean;
  double chi=0.;
  if (nFit<=1) return chi;
  getMean(mean);
  for (list<Detection*>::const_iterator i=elist.begin(); 
       i!=elist.end(); 
       ++i) {
    if (!isFit(*i)) continue;
    double mi = (*i)->magOut;
    double wti = (*i)->wt;
    double cc = (mi-mean)*(mi-mean)*wti;
    maxDeviateSq = MAX(cc , maxDeviateSq);
    chi += cc;
  }
  dof += nFit-1;
  return chi;
}

void
PhotoAlign::operator()(const DVector& p, double& chisq,
		       DVector& beta, tmv::SymMatrix<double>& alpha) {
  countPriorParams();
  int nP = nParams();
  Assert(p.size()==nP);
  Assert(beta.size()==nP);
  Assert(alpha.nrows()==nP);
  setParams(p);
  double newChisq=0.;
  beta.setZero();
  alpha.setZero();
  int matchCtr=0;

  const int NumberOfLocks = 2000;
  AlphaUpdater updater(alpha, pmc.nFreeMaps(), NumberOfLocks);

#ifdef _OPENMP
  const int chunk=10000;
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
      for (list<Match*>::const_iterator i=mlist.begin();
	   i != mlist.end();
	   ++i) {
	vi[j] = *i;
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
  for (int j=0; j<vi.size(); j++) {
    Match* m = vi[j];

    //**    if (j%chunk==0) {
    if (false) {
#pragma omp critical (output)
      cerr << "# Thread " << omp_get_thread_num()
	   << "# accumulating chisq at match # " << j 
	   << endl;
    }
    if ( m->getReserved() ) continue;	//skip reserved objects
    //**    m->accumulateChisq(newChisq, newBeta, alpha);
    m->accumulateChisq(newChisq, newBeta, updater);
  }
#pragma omp critical(beta)
  beta += newBeta;
  }
#else
  // This is the match loop for single threading:
  for (list<Match*>::const_iterator i=mlist.begin();
       i != mlist.end();
       ++i) {
    Match* m = *i;
    if (matchCtr%10000==0) cerr << "# accumulating chisq at match # " 
				<< matchCtr  //**<< " newChisq " << newChisq
				<< endl;
    matchCtr++;
    if ( m->getReserved() ) continue;	//skip reserved objects
    m->accumulateChisq(newChisq, beta, updater);
  }
#endif
  chisq = newChisq;

  // Now need to accumulate the contributions from the priors.  A single thread will do.
  for (list<PhotoPrior*>::iterator i = priors.begin();
       i != priors.end();
       ++i) 
    (*i)->accumulateChisq(chisq, beta, updater);

  {
    //***** Code that will be useful in spotting unconstrained parameters:
    for (int i = 0; i<alpha.nrows(); i++) {
      bool blank = true;
      for (int j=0; j<alpha.ncols(); j++)
	if (alpha(i,j)!=0.) {
	  blank = false;
	  break;
	}
      if (blank) 
	cerr << "***No constraints on row " << i << endl;
    }
    for (int j=0; j<alpha.ncols(); j++) {
      bool blank = true;
      for (int i = 0; i<alpha.nrows(); i++) 
	if (alpha(i,j)!=0.) {
	  blank = false;
	  break;
	}
      if (blank) 
	cerr << "***No constraints on col " << j << endl;
    }
  }
}

double
PhotoAlign::fitOnce(bool reportToCerr) {
  DVector p = getParams();
  Marquardt<PhotoAlign> marq(*this);
  marq.setRelTolerance(relativeTolerance);
  marq.setSaveMemory();
  double chisq = marq.fit(p, DefaultMaxIterations, reportToCerr);
  setParams(p);
  remap();
  return chisq;
}

void
PhotoAlign::remap() {
  for (list<Match*>::const_iterator i=mlist.begin();
       i != mlist.end();
       ++i)
    (*i)->remap();

  for (list<PhotoPrior*>::iterator i = priors.begin();
       i != priors.end();
       ++i) 
    (*i)->remap();
}

int
PhotoAlign::sigmaClip(double sigThresh, bool doReserved, bool clipEntireMatch) {
  int nclip=0;
  // ??? parallelize this loop!
  cerr << "## Sigma clipping...";
  Stopwatch timer;
  timer.start();

  for (list<Match*>::const_iterator i=mlist.begin();
       i != mlist.end();
       ++i) {
    // Skip this one if it's reserved and doReserved=false,
    // or vice-versa
    if (doReserved ^ (*i)->getReserved()) continue;
    if ( (*i)->sigmaClip(sigThresh)) {
      nclip++;
      if (clipEntireMatch) (*i)->clipAll();
    }
  }
  timer.stop();
  cerr << " done in " << timer << " sec" << endl;
  return nclip;
}

int
PhotoAlign::sigmaClipPrior(double sigThresh, bool clipEntirePrior) {
  int nClip = 0;
  for (list<PhotoPrior*>::iterator i = priors.begin();
       i != priors.end();
       ++i) 
    if ((*i)->sigmaClip(sigThresh)) {
      // Clipped something.
      nClip++;
      if (clipEntirePrior) (*i)->clipAll();
    }
  return nClip;
}

double
PhotoAlign::chisqDOF(int& dof, double& maxDeviate, 
		     bool doReserved) const {
  dof=0;
  maxDeviate=0.;
  double chisq=0.;
  // ??? parallelize this loop!
  for (list<Match*>::const_iterator i=mlist.begin();
       i != mlist.end();
       ++i) {
    // Skip this one if it's reserved and doReserved=false,
    // or vice-versa
    if (doReserved ^ (*i)->getReserved()) continue;
    chisq += (*i)->chisq(dof, maxDeviate);
  }
  maxDeviate = sqrt(maxDeviate);
  if (!doReserved) {
    // If doing the fitted objects, include prior and adjust DOF for fit
    dof -= pmc.nParams();
    for (list<PhotoPrior*>::const_iterator i = priors.begin();
	 i != priors.end();
	 ++i) {
      chisq += (*i)->chisq(dof);
    }
  }
  return chisq;
}

void 
PhotoAlign::count(long int& mcount, long int& dcount, 
		  bool doReserved, int minMatches) const {
  mcount = 0;
  dcount = 0;
  for (list<Match*>::const_iterator i=mlist.begin();
       i != mlist.end();
       ++i) {
    if (((*i)->getReserved() ^ doReserved) 
	|| (*i)->fitSize() < minMatches) continue;
    mcount++;
    dcount+=(*i)->fitSize();
  }
}
void
PhotoAlign::count(long int& mcount, long int& dcount,
                  bool doReserved, int minMatches, long catalog) const {
  mcount = 0;
  dcount = 0;
  for (list<Match*>::const_iterator i=mlist.begin();
       i != mlist.end();
       ++i) {
    if (((*i)->getReserved() ^ doReserved)
        || (*i)->fitSize() < minMatches) continue;

    int ddcount=0;
    for(list<Detection*>::iterator d=(*i)->begin(); d != (*i)->end(); ++d) {
      if(!((*d)->isClipped) && (*d)->catalogNumber==catalog)
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
  for (list<PhotoPrior*>::const_iterator i = priors.begin();
       i != priors.end();
       ++i) {
    if ((*i)->isDegenerate()) continue;
    (*i)->globalStartIndex = startIndex;
    startIndex += (*i)->nParams();
    (*i)->globalMapNumber = mapNumber++;
  }
  nPriorParams = startIndex - pmc.nParams();
}

void
PhotoAlign::setParams(const DVector& p) {
  pmc.setParams(p.subVector(0,pmc.nParams()));
  int startIndex = pmc.nParams();
  for (list<PhotoPrior*>::iterator i = priors.begin();
       i != priors.end();
       ++i) {
    if ((*i)->isDegenerate()) continue;
    (*i)->setParams(p.subVector(startIndex, startIndex+(*i)->nParams()));
    startIndex += (*i)->nParams();
  }
}

DVector
PhotoAlign::getParams() const {
  DVector p(nParams(), -888.);
  p.subVector(0,pmc.nParams()) = pmc.getParams();
  int startIndex = pmc.nParams();
  for (list<PhotoPrior*>::const_iterator i = priors.begin();
       i != priors.end();
       ++i) {
    if ((*i)->isDegenerate()) continue;
    p.subVector(startIndex, startIndex+(*i)->nParams()) = (*i)->getParams();
    startIndex += (*i)->nParams();
  }
}
