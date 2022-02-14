/** A class that permits efficient parallel construction of large symmetric matrix
    as sum of outer products, when each outer product exists over short ranges
    of indices.  The updates are multithreaded when OpenMP is available.
    The index range is divided into disjoint "segments" and there is a locking mechanism
    to avoid having more than one thread trying to update a given submatrix at 
    the same time.
    It is required that the segment numbers be orders so that the indices 0...i1 
    are in segment 0, i1+1..i2 are in segment 1, etc.
    The matrix being managed is a regular matrix (either TMV or Eigen) but we
    are only going to fill in the LOWER triangle elements (i,j) with j<=i.
**/
#ifndef SYMMETRICUPDATER_H
#define SYMMETRICUPDATER_H

#ifdef _OPENMP
#include <omp.h>
#endif

#include "LinearAlgebra.h"
#include <vector>

namespace linalg {
  class SymmetricUpdater {
  public:
  SymmetricUpdater(linalg::Matrix<double>& alpha_,
		   int nSegs_,
		   int nLocks=0): alpha(alpha_), nSegs(nSegs_) {
#ifdef _OPENMP
      // Decide how many segments per lock
      int nBlocks = nSegs*(nSegs+1)/2;
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

    ~SymmetricUpdater() {
#ifdef _OPENMP    // Release locks
      for (int i=0; i<locks.size(); i++)
	omp_destroy_lock(&locks[i]);
#endif
    }

    // Update submatrix with corner at (start1,start2)
    // as alpha += scalar * v1 * v2^T
    // seg[12] are numbers of the index segments in which v1 and v2 must lie.
    // First version must have seg1!=seg2, i.e. does not cross diagonal
    // Order will be swapped if needed to put this in lower triangle of alpha.
    template<class V>
    void rankOneUpdate(int seg1, int start1, const V& v1, 
		       int seg2, int start2, const V& v2,
		       double scalar=1.) {
      if (v1.size() <= 0 || v2.size() <=0) return; // Nothing to add.
      if (seg1==seg2)
	throw std::runtime_error("SymmetricUpdater got 2 equal segments");
      if (seg1 < seg2) {
	rankOneUpdate(seg2,start2,v2,
		      seg1,start1,v1, scalar);  // Swap so seg1>=seg2.
	return;
      }
#ifdef _OPENMP
      // Get the block index corresponding to these 2 keys
      int block = (seg1*(seg1+1)/2 + seg2)/ blockLength;
      Assert(block < locks.size());
      // Set lock for its region - will wait here if busy!
      omp_set_lock(&locks[block]);
#endif
      // Use the native outer-product language for each package in
      // case there are efficiencies.
#ifdef USE_TMV
      if (scalar==1.) 
	alpha.subMatrix(start1, start1+v1.size(),
			start2, start2+v2.size()) += v1 ^ v2;
      else
	alpha.subMatrix(start1, start1+v1.size(),
			start2, start2+v2.size()) += scalar * v1 ^ v2;
#elif defined USE_EIGEN
      if (scalar==1.) 
	alpha.subMatrix(start1, start1+v1.size(),
			start2, start2+v2.size()) += v1 * v2.transpose();
      else
	alpha.subMatrix(start1, start1+v1.size(),
			start2, start2+v2.size()) += scalar * v1 * v2.transpose();
#endif
#ifdef _OPENMP
      omp_unset_lock(&locks[block]);
#endif
    }
    // Second version does += scalar * v * v^T, which crosses diagonal
    template<class V>
    void rankOneUpdate(int seg, int start, const V& v, 
		       double scalar=1.) {
#ifdef _OPENMP
      // Get the block index corresponding to these 2 keys
      int block = (seg*(seg+3)/2)/ blockLength;
      Assert(block < locks.size());
      // Set lock for its region - will wait here if busy!
      omp_set_lock(&locks[block]);
#endif
      // Use each package's special routine
#ifdef USE_TMV
      tmv::SymMatrix<double> tmp = v ^ v;
      if (scalar!=1.) tmp *= scalar;
      tmv::SymMatrixViewOf(alpha.subMatrix(start, start+v.size(),
					   start, start+v.size()),
			   tmv::Lower) += tmp;
#elif defined USE_EIGEN
      alpha.subMatrix(start, start+v.size(),start, start+v.size()).
	template selfadjointView<Eigen::Lower>().
	rankUpdate(v,scalar);
#endif
#ifdef _OPENMP
      omp_unset_lock(&locks[block]);
#endif
    }

    // Update submatrix with corner at (start1,start2)
    // as alpha +=  M1^T * SM * M2 
    // where SM is an m x m symmetric matrix,
    // M1 is m x n1, M2 is m x n2, and
    // (start1,start1+n1) is disjoint from (start2,start2+n2)
    // so that update does not cross diagonal.
    // seg[12] are numbers of the index segments in which v1 and v2 must lie.
    // First version must have seg1!=seg2, i.e. does not cross diagonal
    // Order will be swapped if needed to put this in lower triangle of alpha,
    // i.e. we will use M2^T * SM * M1 as needed.
    template<class M1, class SM, class M2>
    void update(int seg1, int start1, const M1& m1, 
		const SM& sm,
		int seg2, int start2, const M2& m2) {
      int n1 = m1.cols();
      int n2 = m2.cols();
      if (n1 <= 0 || n2 <=0) return; // Nothing to add.
      if (seg1==seg2)
	throw std::runtime_error("SymmetricUpdater got 2 equal segments");
      if (seg1 < seg2) {
	update(seg2,start2,m2,
	       sm,
	       seg1,start1,m1);  // Swap so seg1>=seg2.
	return;
      }
#ifdef _OPENMP
      // Get the block index corresponding to these 2 keys
      int block = (seg1*(seg1+1)/2 + seg2)/ blockLength;
      Assert(block < locks.size());
      // Set lock for its region - will wait here if busy!
      omp_set_lock(&locks[block]);
#endif

      alpha.subMatrix(start1, start1+n1, start2, start2+n2) +=
	m1.transpose() * sm * m2;

#ifdef _OPENMP
      omp_unset_lock(&locks[block]);
#endif
    }
  
    // This updates with M^T * SM * M,
    // where SM is m x m symmetric, and
    // M is m x n matrix.  So we cross diagonal
    template<class M, class SM>
    void update(int seg, int start, const M& m, 
		const SM& sm) {
      int n = m.cols();
      if (n<=0) return;
      
#ifdef _OPENMP
      // Get the block index corresponding to these 2 keys
      int block = (seg*(seg+3)/2)/ blockLength;
      Assert(block < locks.size());
      // Set lock for its region - will wait here if busy!
      omp_set_lock(&locks[block]);
#endif
      // Use package's symmetric-matrix routine if possible
#ifdef USE_TMV
      tmv::SymMatrix<double> tmp = m.transpose() * sm * m;
      tmv::SymMatrixViewOf(alpha.subMatrix(start, start+n,
					   start, start+n),
			   tmv::Lower) += tmp;
#elif defined USE_EIGEN
      DMatrix tmp = m.transpose() * sm * m; // Can we save 2x???
      alpha.subMatrix(start, start+n,start, start+n)
	+= tmp;
#endif
#ifdef _OPENMP
      omp_unset_lock(&locks[block]);
#endif
    }

  private:
    linalg::Matrix<double>& alpha;
    int nSegs;
#ifdef _OPENMP
    int nLocks;
    int blockLength;
    // Lock array
    std::vector<omp_lock_t> locks;
#endif
  };

} // end namespace linalg

#endif  // SYMMETRICUPDATER_H
