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

#ifdef _OPENMP
#include <omp.h>
#endif

// #define DEBUG
#include "Marquardt.h"

using namespace astrometry;


/////////////////////////////////////////////////////////////
///// Class that regulates rank-one updates to small portions of
///// the giant alpha symmetric matrix, allowing finer-grained locking
/////////////////////////////////////////////////////////////
class astrometry::AlphaUpdater {
public:
  AlphaUpdater(tmv::SymMatrix<double>& alpha_,
	       const PixelMapCollection& pmc_,
	       int nLocks=0): alpha(alpha_), pmc(pmc_) {
#ifdef _OPENMP
    // Decide how many map elements per lock
    nMaps = pmc.nFreeMaps();
    int nBlocks = nMaps*(nMaps+1)/2;
    if (nLocks > nBlocks || nLocks<=0) {
      blockLength = 1;
      nLocks = nBlocks;
    } else {
      blockLength = nBlocks / nLocks;
      nLocks = (nBlocks+blockLength-1) / blockLength;
    }
    locks.resize(nLocks);
    // build array of locks
    for (int i=0; i<locks.size(); i++)
      omp_init_lock(&locks[i]);
#endif
  }
  ~AlphaUpdater() {
#ifdef _OPENMP    // Release locks
    for (int i=0; i<locks.size(); i++)
      omp_destroy_lock(&locks[i]);
#endif
  }
  // Update submatrix with corner at (startIndex1,startIndex2) with += scalar * v1 ^ v2
  // map[12] are sequence numbers of map components in the PixelMapCollection, used
  // to divide the matrix into blocks for locking in multithreaded case
  void rankOneUpdate(int map1, int startIndex1, tmv::VectorView<double>& v1, 
		     int map2, int startIndex2, tmv::VectorView<double>& v2,
		     double scalar=1.) {
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
private:
  tmv::SymMatrix<double>& alpha;
  const PixelMapCollection& pmc;
#ifdef _OPENMP
  int nLocks;
  int blockLength;
  int nMaps;
  // Lock array
  vector<omp_lock_t> locks;
#endif
};

bool
Match::isFit(const Detection* e) {
  // A Detection will contribute to fit if it has nonzero weight and is not clipped.
  return !(e->isClipped) && !( e->wtx==0. && e->wty==0.);
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
Match::countFit() {
  nFit = 0;
  for (list<Detection*>::iterator i=elist.begin();
       i!=elist.end();
       ++i)
    if (isFit(*i)) nFit++;
}

void
Match::clipAll() {
  for (list<Detection*>::iterator i=elist.begin();
       i!=elist.end();
       ++i) 
    (*i)->isClipped = true;
  nFit = 0;
}

void
Match::centroid(double& x, double& y) const {
  double wtx, wty;
  centroid(x,y,wtx,wty);
}

void
Match::centroid(double& x, double& y,
		double& wtx, double& wty) const {
  if (nFit<1) {
    wtx=wty=x=y=0.;
    return;
  }
  double swx, swy;
  swx = swy =0.;
  wtx = wty = 0.;
  for (list<Detection*>::const_iterator i=elist.begin();
       i!=elist.end();
       ++i) {
    if (!isFit(*i)) continue;
    double xx = (*i)->xw;
    double yy = (*i)->yw;
    double wx = (*i)->wtx;
    double wy = (*i)->wty;
    swx += xx*wx;
    wtx += wx;
    swy += yy*wy;
    wty += wy;
  }
  x = swx/wtx;
  y = swy/wty;
}


///////////////////////////////////////////////////////////
// Coordinate-matching routines
///////////////////////////////////////////////////////////

void
Match::remap() {
  for (list<Detection*>::iterator i=elist.begin(); 
       i!=elist.end(); 
       ++i) 
    (*i)->map->toWorld((*i)->xpix, (*i)->ypix,
		       (*i)->xw, (*i)->yw);
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
		       AlphaUpdater& updater,
		       bool reuseAlpha) {
  double xmean, ymean;
  double xW, yW;
  int nP = beta.size();

  // No contributions to fit for <2 detections:
  if (nFit<=1) return 0;

  // Update mapping and save derivatives for each detection:
  vector<DMatrix*> dxyi(elist.size(), 0);
  int ipt=0;
  for (list<Detection*>::iterator i=elist.begin(); 
       i!=elist.end(); 
       ++i, ipt++) {
    if (!isFit(*i)) continue;
    int npi = (*i)->map->nParams();
    double xw, yw;
    if (npi>0) {
      dxyi[ipt] = new DMatrix(2,npi);
      (*i)->map->toWorldDerivs((*i)->xpix, (*i)->ypix,
			       xw, yw,
			       *dxyi[ipt]);
    } else {
      (*i)->map->toWorld((*i)->xpix, (*i)->ypix, xw, yw);
    }
    (*i)->xw = xw;
    (*i)->yw = yw;
  }

  centroid(xmean,ymean, xW, yW);
  DVector dxmean(nP, 0.);
  DVector dymean(nP, 0.);


  map<int, iRange> mapsTouched;
  ipt = 0;
  for (list<Detection*>::iterator i=elist.begin(); 
       i!=elist.end(); 
       ++i, ipt++) {
    if (!isFit(*i)) continue;
    double wxi=(*i)->wtx;
    double wyi=(*i)->wty;
    double xi=(*i)->xw;
    double yi=(*i)->yw;

    chisq += 
      (xi-xmean)*(xi-xmean)*wxi
      + (yi-ymean)*(yi-ymean)*wyi;

    // Accumulate derivatives:
    int istart=0;
    for (int iMap=0; iMap<(*i)->map->nMaps(); iMap++) {
      int ip=(*i)->map->startIndex(iMap);
      int np=(*i)->map->nSubParams(iMap);
      if (np==0) continue;
      int mapNumber = (*i)->map->mapNumber(iMap);
      // Keep track of parameter ranges we've messed with:
      mapsTouched[mapNumber] = iRange(ip,np);
      tmv::VectorView<double> dx=dxyi[ipt]->row(0,istart,istart+np);
      tmv::VectorView<double> dy=dxyi[ipt]->row(1,istart,istart+np);
      beta.subVector(ip, ip+np) -= wxi*(xi-xmean)*dx;
      beta.subVector(ip, ip+np) -= wyi*(yi-ymean)*dy;

      // Derivatives of the mean position:
      dxmean.subVector(ip, ip+np) += wxi*dx;
      dymean.subVector(ip, ip+np) += wyi*dy;

      if (!reuseAlpha) {
	// Increment the alpha matrix
	//**      alpha.subSymMatrix(ip, ip+np) += wxi*(dx ^ dx));
	//**      alpha.subSymMatrix(ip, ip+np) += wxi*(dy ^ dy));
	updater.rankOneUpdate(mapNumber, ip, dx, 
			      mapNumber, ip, dx, wxi);
	updater.rankOneUpdate(mapNumber, ip, dy, 
			      mapNumber, ip, dy, wyi);

	int istart2 = istart+np;
	for (int iMap2=iMap+1; iMap2<(*i)->map->nMaps(); iMap2++) {
	  int ip2=(*i)->map->startIndex(iMap2);
	  int np2=(*i)->map->nSubParams(iMap2);
	  int mapNumber2 = (*i)->map->mapNumber(iMap2);
	  if (np2==0) continue;
	  tmv::VectorView<double> dx2=dxyi[ipt]->row(0,istart2,istart2+np2);
	  tmv::VectorView<double> dy2=dxyi[ipt]->row(1,istart2,istart2+np2);
	  // Note that subMatrix here will not cross diagonal:
	  //**	  alpha.subMatrix(ip, ip+np, ip2, ip2+np2) += wxi*(dx ^ dx2);
	  //**	  alpha.subMatrix(ip, ip+np, ip2, ip2+np2) += wyi*(dy ^ dy2);
	  updater.rankOneUpdate(mapNumber2, ip2, dx2, 
				mapNumber,  ip,  dx, wxi);
	  updater.rankOneUpdate(mapNumber2, ip2, dy2, 
				mapNumber,  ip,  dy, wyi);
	  istart2+=np2;
	}
      }
      istart+=np;
    } // outer parameter segment loop
    if (dxyi[ipt]) {delete dxyi[ipt]; dxyi[ipt]=0;}
  } // object loop

  if (!reuseAlpha) {
    // Subtract effects of derivatives on mean
    /*  We want to do this, but without touching the entire alpha matrix every time:
	alpha -=  (dxmean ^ dxmean)/xW;
	alpha -=  (dymean ^ dymean)/yW;
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
	tmv::VectorView<double> dx1 = dxmean.subVector(i1, i1+n1);
	tmv::VectorView<double> dx2 = dxmean.subVector(i2, i2+n2);
	tmv::VectorView<double> dy1 = dymean.subVector(i1, i1+n1);
	tmv::VectorView<double> dy2 = dymean.subVector(i2, i2+n2);
	updater.rankOneUpdate(map2, i2, dx2,
			      map1, i1, dx1, -1./xW);
	updater.rankOneUpdate(map2, i2, dy2,
			      map1, i1, dy1, -1./yW);
      }
    }
  } // Finished putting terms from mean into alpha

  return 2*(nFit-1);
}

bool
Match::sigmaClip(double sigThresh,
		 bool deleteDetection) {
  // Only clip the worst outlier at most
  double xmean, ymean;
  if (nFit<=1) return false;
  double maxSq=0.;
  list<Detection*>::iterator worst=elist.end();
  centroid(xmean,ymean);
  for (list<Detection*>::iterator i=elist.begin(); 
       i!=elist.end(); ++i) {
    if (!isFit(*i)) continue;
    double xi = (*i)->xw;
    double yi = (*i)->yw;
    double clipx = (*i)->clipsqx;
    double clipy = (*i)->clipsqy;
    double devSq = (xi-xmean)*(xi-xmean)*clipx
      + (yi-ymean)*(yi-ymean)*clipy;
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
	 << " resid " << setw(6) << ((*worst)->xw-xmean)*DEGREE/ARCSEC
	 << " " << setw(6) << ((*worst)->yw-ymean)*DEGREE/ARCSEC
	 << " sigma " << DEGREE/ARCSEC/sqrt((*worst)->clipsqx)
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
  double xmean, ymean;
  double chi=0.;
  if (nFit<=1) return chi;
  centroid(xmean,ymean);
  for (list<Detection*>::const_iterator i=elist.begin(); 
       i!=elist.end(); 
       ++i) {
    if (!isFit(*i)) continue;
    double xi = (*i)->xw;
    double yi = (*i)->yw;
    double wxi = (*i)->wtx;
    double wyi = (*i)->wty;
    double cc = (xi-xmean)*(xi-xmean)*wxi
      + (yi-ymean)*(yi-ymean)*wyi;
    maxDeviateSq = MAX(cc , maxDeviateSq);
    chi += cc;
  }
  dof += 2*(nFit-1);
  return chi;
}

void
CoordAlign::operator()(const DVector& p, double& chisq,
		       DVector& beta, tmv::SymMatrix<double>& alpha,
		       bool reuseAlpha) {
  int nP = pmc.nParams();
  Assert(p.size()==nP);
  Assert(beta.size()==nP);
  Assert(alpha.nrows()==nP);
  setParams(p);
  double newChisq=0.;
  beta.setZero();
  alpha.setZero();
  int matchCtr=0;

  const int NumberOfLocks = 2000;
  AlphaUpdater updater(alpha, pmc, NumberOfLocks);

#ifdef _OPENMP
  const int chunk=2000;
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
      for (list<Match*>::const_iterator i=mlist.begin();
	   i != mlist.end();
	   ++i) {
	vi[j] = *i;
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
    m->accumulateChisq(newChisq, newBeta, updater, reuseAlpha);
  }
#pragma omp critical(beta)
  beta += newBeta;
  }
#else
  for (list<Match*>::const_iterator i=mlist.begin();
       i != mlist.end();
       ++i) {
    Match* m = *i;
    if (matchCtr%10000==0) cerr << "# accumulating chisq at match # " 
				<< matchCtr  //**<< " newChisq " << newChisq
				<< endl;
    matchCtr++;
    if ( m->getReserved() ) continue;	//skip reserved objects
    m->accumulateChisq(newChisq, beta, updater, reuseAlpha);
  }
#endif
  chisq = newChisq;

  if (!reuseAlpha) {
    //***** Code that will be useful in spotting unconstrained parameters:
    for (int i = 0; i<alpha.nrows(); i++) {
      bool blank = true;
      for (int j=0; j<alpha.ncols(); j++)
	if (alpha(i,j)!=0.) {
	  blank = false;
	  break;
	}
      if (blank) {
	cerr << "***No constraints on row " << i << endl;
	string badAtom = pmc.atomHavingParameter(i);
	cerr << "Serialized version of the degenerate map:" << endl;
	pmc.writeMap(cerr, badAtom);
	cerr << "Startindex is " << pmc.issueMap(badAtom)->startIndex(0) <<endl;
      }
    }
  }
}

double
CoordAlign::fitOnce(bool reportToCerr) {
  DVector p = getParams();
  // First will try doing Newton iterations, keeping a fixed Hessian.
  // If it increases chisq or takes too long to converge, we will 
  // go into the Marquardt solver.

  {
    // Get chisq, beta, alpha at starting parameters
    double oldChisq = 0.;
    int nP = pmc.nParams();
    DVector beta(nP, 0.);
    tmv::SymMatrix<double> alpha(nP, 0.);

    Stopwatch timer;
    timer.start();
    (*this)(p, oldChisq, beta, alpha);
    timer.stop();
    /**/cerr << "fitOnce starting chisq " << oldChisq 
	     << " alpha time " << timer << endl;
    timer.reset();
    timer.start();

    // Set alpha for in-place Cholesky
    alpha.divideUsing(tmv::CH);   // Use Cholesky decomposition for fastest inversion
    alpha.divideInPlace();  // In-place avoids doubling memory needs

    const int MAX_NEWTON_STEPS = 8;
    int newtonIter = 0;
    for (int newtonIter = 0; newtonIter < MAX_NEWTON_STEPS; newtonIter++) {
      beta /= alpha;
      DVector newP = p + beta;
      setParams(newP);
      // Get chisq at the new parameters
      remap();
      int dof;
      double maxDev;
      double newChisq = chisqDOF(dof, maxDev);
      timer.stop();
      /**/cerr << "Newton iteration #" << newtonIter << " yields chisq " << newChisq 
	       << " / " << dof 
	       << " in time " << timer << " sec"
	       << endl;
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
  for (list<Match*>::const_iterator i=mlist.begin();
       i != mlist.end();
       ++i)
    (*i)->remap();
}

int
CoordAlign::sigmaClip(double sigThresh, bool doReserved, bool clipEntireMatch) {
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

double
CoordAlign::chisqDOF(int& dof, double& maxDeviate, 
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
  if (!doReserved) dof -= pmc.nParams();
  return chisq;
}

void 
CoordAlign::count(long int& mcount, long int& dcount, 
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
CoordAlign::count(long int& mcount, long int& dcount,
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


