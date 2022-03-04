// Implementation of the Marquardt-Levenberg nonlinear fitting algorithm
// in C++, and using Jarvis TMV
// The algorithm is taken from Numerical Recipes code, but has been 
// rewritten to handle more general kinds of chi-squared.

// It's a template so that classes instead of just functions can be used
// in the evaluation loops.  Not too worried about code bloat since it's
// not likely this will have multiple instantiations in one program.

// ??? Make the convergence decision a separate class and implement
// a default for it.

#ifndef MARQUARDT_H
#define MARQUARDT_H

//#define DEBUG

#include <stdexcept>
#include "LinearAlgebra.h"
#include "Std.h"
#include "Stopwatch.h"

#ifdef USE_TMV
#include "TMV_Sym.h"
#elif defined USE_EIGEN
#include "Eigen/Cholesky"
#endif


// An exception class:
class NoConverge: public std::runtime_error {
public:
 NoConverge(const string m): std::runtime_error("Marquardt did not converge; "+m)
  {}
};

const int DefaultMaxIterations=100;
const double DefaultAbsTolerance=0.05;
const double DefaultRelTolerance=1e-4;
const double MaxLambda=1e6;

template <class P>
class MarquardtPointSum;

// Marquardt minimizer will call the object of class T with args
// (const Vector& a, double& chisq, Vector& beta, Matrix& alpha))
// to build the alpha/beta/chisq from parameters a.  
// P is precision of calculation, defaults to double.
// Only the lower triangle of alpha is used, symmetry assumed.
// MarquandtPointSum class is below, when you have a point function.

// Call method setSaveMemory(true) to instruct routine to do matrix
// solution in-place and not store a copy of the best-previous alpha
// matrix.  This means we'll have to recompute it if we have a non-improving
// iteration or want to know alpha after fitting.

template <class T=MarquardtPointSum<double>, class P=double>
class Marquardt {
public: 
  typedef linalg::Vector<P> Vector;
  typedef linalg::Matrix<P> Matrix;
  
  explicit Marquardt(T& f):
    derivs(f), isFit(false), absTol(DefaultAbsTolerance),
    relTol(DefaultRelTolerance), bestA(0), bestAlpha(0),
    saveMemory(false) {}
  // No copying:
  Marquardt(const Marquardt& rhs) =delete;
  void operator=(const Marquardt& rhs) =delete;
  ~Marquardt() {
    if (bestAlpha) delete bestAlpha;
    if (bestA) delete bestA;
  }

  // Does the fit starting at a, returns chisq.  Set the flag to get progress/timing to cerr
  P fit(Vector& a, int maxIter=DefaultMaxIterations, bool progressToCerr=false);
  void setSaveMemory(bool saveMemory_=true) {saveMemory=saveMemory_;}

  //Return (reference to) inverse covariance matrix at last fit
  const Matrix& getAlpha();
  void setAbsTolerance(P t) {absTol=t;}
  void setRelTolerance(P t) {relTol=t;}
  P getAbsTolerance() const {return absTol;}
  P getRelTolerance() const {return relTol;}

private:
  // The object returning alpha/beta:
  T& derivs;

  // Convergence criteria:
  P absTol;
  P relTol;
  int    maxIter;

  // Saved best fit conditions:
  Vector *bestA;
  Matrix *bestAlpha;		//we're going to save our best alpha
  P  bestChisq;

  bool isFit;		//set flag if bestA is the best fit.
  bool lastDropWasSmall;//flag used to find 2 small drops in a row.
  bool saveMemory;	// Set to minimize memory footprint at cost of recomputations
  bool isConverged(P chisq);
};

template <class T, class P>
const linalg::Matrix<P>&
Marquardt<T,P>::getAlpha() {
  if (saveMemory || !bestAlpha) {
    // Recalculate bestAlpha with current best params if there is not a stored one:
    if (!bestA) throw std::runtime_error("Marquardt::getAlpha() called with no stored solution");
    int nparam = bestA->size();
    if (bestAlpha && bestAlpha->nrows()==nparam) {
      // Already have an array for bestAlpha, do nothing
    } else {
      // Make new matrix for bestAlpha
      if (bestAlpha) delete bestAlpha;
      bestAlpha = new Matrix(nparam,nparam);
    }
    Vector bestBeta(nparam);
    //Call subroutine that accumulates alpha, beta
    derivs(*bestA,bestChisq,bestBeta,*bestAlpha);
  } 
  // Copy the lower triangle into upper before making matrix public:
  for (int i=0; i<bestA->size(); i++) 
    for (int j=0; j<i; j++)
      (*bestAlpha)(j,i) = (*bestAlpha)(i,j);
  return *bestAlpha;
}

template <class T, class P>
P
Marquardt<T,P>::fit(Vector& a, int maxIter, bool progressToCerr) {
  P lambda;
  int nparam = a.size();

  Stopwatch timer;

  // These are used temporarily during fit():
  P chisq;
  Vector *bestBeta=0;
  // Discard any old info
  if (bestAlpha) delete bestAlpha;
  // get matrix to save results if we are not conserving memory
  if (!saveMemory) {
    bestAlpha = new Matrix(nparam,nparam);
    bestBeta = new Vector(nparam);
  }

  auto beta = new Vector(nparam);
  auto alpha= new Matrix(nparam,nparam);

  if (bestA) delete bestA;
  bestA = new Vector(a);
  isFit = false;
  lastDropWasSmall=false;

  //Make normal matrices first time through:
  if (progressToCerr) {
    timer.start();
    cerr << "## Marquart iteration 0 derivatives...";
  }
  if (saveMemory) 
    derivs(a,bestChisq,*beta,*alpha);
  else
    derivs(a,bestChisq,*bestBeta,*bestAlpha);
  if (progressToCerr) {
    timer.stop();
    cerr << "done in " << to_string(timer) << " sec" 
	 << " chisq=" << to_string(bestChisq) 
	 << endl;
  }

  lambda = 0.001;

  // If maxIter=0, the client just wanted to calculate the alpha.
  if (maxIter==0) {
    if (saveMemory) {
      bestAlpha = alpha;
    } else {
      delete alpha;
    }
    delete beta; 
    if (bestBeta) delete bestBeta;
    return bestChisq;
  }

  // Iterate to convergence:
  for (int i=0; i<maxIter; i++) {
    // Calculate the differential step
    if (!saveMemory) {
      *alpha = *bestAlpha;
      *beta =  *bestBeta;
    }

    for (int j=0; j<a.size(); j++)
      (*alpha)(j,j) *= (1+lambda);

    if (progressToCerr) {
      timer.reset();
      timer.start();
      cerr << "## Marquart iteration " << to_string(i) << " inversion...";
    }

    // Solve linear system and get new parameters.
#ifdef USE_TMV
    {
      // Make a symmetric view into alpha for Cholesky:
      auto s = tmv::SymMatrixViewOf(*alpha,tmv::Lower);
      s.divideUsing(tmv::CH);
      if (saveMemory) s.divideInPlace();
      *beta /= s;
    } 
    a = *bestA + *beta;
#elif defined USE_EIGEN
    if (saveMemory) {
      // In-place solution of alpha
      Eigen::LLT<Eigen::Ref<typename Matrix::Base>> llt(*alpha);
      a = *bestA + llt.solve(*beta);
    } else {
      a = *bestA + alpha->llt().solve(*beta);
    }
#endif
    
    if (progressToCerr) {
      timer.stop();
      cerr << "done in " << to_string(timer) << " sec" << endl;
    }


    // Get chisq and derivs at the new spot
    if (progressToCerr) {
      timer.reset();
      timer.start();
      cerr << "## Marquart iteration " << to_string(i+1) << " derivatives...";
    }

    derivs(a,chisq,*beta,*alpha);

    if (progressToCerr) {
      timer.stop();
      cerr << "done in " << to_string(timer) << " sec" 
	   << " lambda=" << to_string(lambda) 
	   << " chisq=" << to_string(chisq) 
	   << endl;
    }

#ifdef DEBUG
    cerr << "Marquardt Trial " << to_string(i) 
	 << " lambda=" << to_string(lambda) 
	 << " chisq=" << to_string(chisq) 
	 << " best chisq=" << to_string(bestChisq) 
	 << endl;
#endif

    // See if we are done
    if ( isConverged(chisq) ) {
      isFit = true;
      bestChisq = chisq;
      *bestA = a;
      if (!saveMemory) SWAP(bestAlpha, alpha);
      delete alpha;
      delete beta;
      if (bestBeta) delete bestBeta;
      return bestChisq;
    }

    // Now decide circumstances of next iteration
    if (chisq <= bestChisq) {
      // An improvement:
      bestChisq = chisq;
      *bestA = a;
      if (!saveMemory) {
	SWAP(alpha, bestAlpha);
	SWAP(beta, bestBeta);
      }
      lambda *= 0.1;	//Move toward quadratic solution
      // Too many times wandering down quadratic solutions, do
      // a restart with more steep descent.
      if (lambda<1e-9) lambda=1.;
    } else {
      lambda *=10.;	//Move toward steepest descent
      if (lambda > MaxLambda) {
	// We're going to declare a local min
	// if tiny steepest-descent steps are not improving:
	cerr << "WARNING: Marquardt lambda=" << to_string(lambda) << endl;
	isFit = true;
	bestChisq = chisq;
	*bestA = a;
	if (!saveMemory) SWAP(bestAlpha, alpha);
	delete alpha;
	delete beta;
	if (bestBeta) delete bestBeta;
	return bestChisq;
      } else if (saveMemory) {
	// Need to recalculate alpha & beta at the best spot because we wrecked it
	derivs(*bestA,chisq,*beta,*alpha);
      }
    }
  }
  // Get here after max iterations
  a = *bestA;
  delete alpha;
  delete beta;
  if (bestBeta) delete bestBeta;
  throw NoConverge("Too many iterations");
}

template <class T, class P>
bool
Marquardt<T,P>::isConverged(P chisq) {
  P diff = bestChisq-chisq;
  if (diff<0.) return false;	//chisq must not increase
  bool isSmall = (diff < absTol) 
    || (bestChisq>0. && diff/bestChisq<relTol);
  if (isSmall && lastDropWasSmall){
    return true;  //two small drops=converge
  }
  lastDropWasSmall=isSmall;
  return false;
}

// This class does what mrqmin() does in NR, namely sums of derivs, etc.,
// over a set of (x,y,sigma) data points. Precision is template arg.
// Create this object with the data and pointer to the derivative function.
template <class P=double>
class MarquardtPointSum {
  typedef linalg::Vector<P> Vector;
  typedef linalg::Matrix<P> Matrix;
private:
  const Vector& x;
  const Vector& y;
  const Vector& sig;
  Vector isig;
  void (*func)(const Vector& a,
	       const P xpt,
	       P& ypt,
	       Vector& dpt);
public:
  MarquardtPointSum(const Vector& _x, 
		    const Vector& _y, 
		    const Vector& _sig,
		    void (*f)(const Vector& a,
			      const P xpt,
			      P& ypt,
			      Vector& dpt)):
    x(_x), y(_y), sig(_sig), func(f), isig(x.size()) {
    if (y.size() < x.size() || sig.size() < x.size())
      throw std::runtime_error("Data or Sigma vector size mismatch in "
			       "MarquandtPointSum");
    for (int i=0; i<x.size(); i++) isig[i] = pow(sig[i],-2.);
  }
  void operator() (const Vector& a, 
		   P& chisq, 
		   Vector& beta, 
		   Matrix& alpha) const {
    chisq = 0;
    beta.setZero();
    alpha.setZero();
    // Define a symmetric matrix view to fill just lower triangle
#ifdef USE_TMV
    auto s = tmv::SymMatrixViewOf(alpha, tmv::Lower);
#elif defined USE_EIGEN
    auto s = alpha.template selfadjointView<Eigen::Lower>();
#endif
    Vector derivs(a.size());
    P yfit;
    for (int k=0; k<x.size(); k++) {
      func(a, x[k], yfit, derivs);
      yfit = y[k] - yfit;
      chisq += yfit*yfit*isig[k];
      beta += (yfit*isig[k])*derivs;
#ifdef USE_TMV
      tmv::Rank1Update<true>(isig[k], derivs, s);
#elif defined USE_EIGEN
      s.rankUpdate(derivs, isig[k]);
#endif
    }
  }
};

#endif   //MARQUARDT_H
