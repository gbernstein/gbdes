// $Id: Match.cpp,v 1.18 2012/08/13 17:54:08 dgru Exp $

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
#include <Marquardt.h>

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
      nLocks = (nBlocks+blockLength+-1) / blockLength;
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

Match::Match(Detection *e): elist(1,e), nFit(1), isReserved(false) {
  e->itsMatch = this;
}
void
Match::add(Detection *e) {
  elist.push_back(e);
  e->itsMatch = this;
  e->isClipped = false;
  nFit++;
}

void
Match::remove(Detection* e) {
  elist.remove(e);  
  e->itsMatch = 0;
  // Recount fittable objects
  countFit();
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
int
Match::countFit() {
  nFit = 0;
  for (list<Detection*>::iterator i=elist.begin();
       i!=elist.end();
       ++i) 
    if ( !(*i)->isClipped && !((*i)->wtx==0. && (*i)->wty==0.)) nFit++;
  return nFit;
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
    if ((*i)->isClipped) continue;
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
		       AlphaUpdater& updater) {
		       //**		       tmv::SymMatrix<double>& alpha) {
  double xmean, ymean;
  double xW, yW;
  int nP = beta.size();
  //**  Assert(alpha.nrows()==nP);

  // No contributions to fit for <2 detections:
  if (nFit<=1) return 0;

  // Update mapping and save derivatives for each detection:
  vector<DMatrix*> dxyi(elist.size(), 0);
  int ipt=0;
  for (list<Detection*>::iterator i=elist.begin(); 
       i!=elist.end(); 
       ++i, ipt++) {
    if ((*i)->isClipped) continue;
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
    if ((*i)->isClipped) continue;
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

      // The alpha matrix: a little messy because of symmetry
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
      istart+=np;
    } // outer parameter segment loop
    if (dxyi[ipt]) {delete dxyi[ipt]; dxyi[ipt]=0;}
  } // object loop

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

  return 2*(nFit-1);
}

bool
Match::sigmaClip(double sigThresh,
		 bool deleteDetection) {
  // ??? Only clip the worst outlier at most
  double xmean, ymean;
  if (nFit<=1) return false;
  double maxSq=0.;
  list<Detection*>::iterator worst=elist.end();
  centroid(xmean,ymean);
  for (list<Detection*>::iterator i=elist.begin(); 
       i!=elist.end(); ++i) {
    if ((*i)->isClipped) continue;
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
    if ((*i)->isClipped) continue;
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
		       DVector& beta, tmv::SymMatrix<double>& alpha) {
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
  for (list<Match*>::const_iterator i=mlist.begin();
       i != mlist.end();
       ++i) {
    Match* m = *i;
    if (matchCtr%10000==0) cerr << "# accumulating chisq at match # " 
				<< matchCtr  //**<< " newChisq " << newChisq
				<< endl;
    matchCtr++;
    if ( m->getReserved() ) continue;	//skip reserved objects
    //**    m->accumulateChisq(newChisq, beta, alpha);
    m->accumulateChisq(newChisq, beta, updater);
  }
#endif
  chisq = newChisq;

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
CoordAlign::fitOnce(bool reportToCerr) {
  DVector p = getParams();
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


