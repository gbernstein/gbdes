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
#include "UseTMV.h"
#include "Std.h"
#include "TMV_Sym.h"
#include "Stopwatch.h"

using tmv::Vector;
using tmv::Matrix;
using tmv::SymMatrix;

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
// MarquandtPointSum class is below, when you have a point function.

// Call method setSaveMemory(true) to instruct routine to do matrix
// solution in-place and not store a copy of the best-previous alpha
// matrix.  This means we'll have to recompute it if we have a non-improving
// iteration or want to know alpha after fitting.

template <class T=MarquardtPointSum<double>, class P=double>
class Marquardt {
public: 
  explicit Marquardt(T& f):
    derivs(f), isFit(false), absTol(DefaultAbsTolerance),
    relTol(DefaultRelTolerance), bestA(0), bestAlpha(0),
    saveMemory(false) {}
  ~Marquardt() {
    if (bestAlpha) delete bestAlpha;
    if (bestA) delete bestA;
  }

  // Does the fit starting at a, returns chisq.  Set the flag to get progress/timing to cerr
  P fit(Vector<P>& a, int maxIter=DefaultMaxIterations, bool progressToCerr=false);
  void setSaveMemory(bool saveMemory_=true) {saveMemory=saveMemory_;}

  //Return (pointer to) inverse covariance matrix at last fit
  const SymMatrix<P>* getAlpha();
  void setAbsTolerance(P t) {absTol=t;}
  void setRelTolerance(P t) {relTol=t;}
  P getAbsTolerance() const {return absTol;}
  P getRelTolerance() const {return relTol;}

private:
  // private copy/assignment to keep it from happening
  Marquardt(const Marquardt& m) {}
  void operator=(const Marquardt& rhs) {}
  // The object returning alpha/beta:
  T& derivs;

  // Convergence criteria:
  P absTol;
  P relTol;
  int    maxIter;

  // Saved best fit conditions:
  Vector<P> *bestA;
  SymMatrix<P> *bestAlpha;		//we're going to save our best alpha
  P  bestChisq;

  bool isFit;		//set flag if bestA is the best fit.
  bool lastDropWasSmall;//flag used to find 2 small drops in a row.
  bool saveMemory;	// Set to minimize memory footprint at cost of recomputations
  bool isConverged(P chisq);
};

template <class T, class P>
const SymMatrix<P>* 
Marquardt<T,P>::getAlpha() {
  if (saveMemory || !bestAlpha) {
    // Recalculate bestAlpha with current best params if there is not a stored one:
    if (!bestA) throw std::runtime_error("Marquardt::getAlpha() called with no stored solution");
    int nparam = bestA->size();
    if (bestAlpha && bestAlpha->size()==nparam) {
      bestAlpha->setZero();
    } else {
      if (bestAlpha) delete bestAlpha;	//get proper sized matrix
      bestAlpha = new SymMatrix<P>(nparam);
    }
    Vector<P> bestBeta(nparam);
    derivs(*bestA,bestChisq,bestBeta,*bestAlpha);	//build all the matrices
  } 
  return bestAlpha;
}

template <class T, class P>
P
Marquardt<T,P>::fit(Vector<P>& a, int maxIter, bool progressToCerr) {
  P lambda;
  int nparam = a.size();

  Stopwatch timer;

  // These are used temporarily during fit():
  P chisq;
  Vector<P> *bestBeta=0;
  // Discard any old info
  if (bestAlpha) delete bestAlpha;
  // get matrix to save results if we are not conserving memory
  if (!saveMemory) {
    bestAlpha = new SymMatrix<P>(nparam);
    bestBeta = new DVector(nparam);
  }

  Vector<P> *beta=new DVector(nparam);
  SymMatrix<P> *alpha=new SymMatrix<P>(nparam);
  alpha->divideUsing(tmv::CH);   // Use Cholesky decomposition for fastest inversion
  if (saveMemory) alpha->divideInPlace();  // In-place avoids doubling memory needs

  if (bestA) delete bestA;
  bestA = new DVector(a);
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
    derivs(a,bestChisq,*bestBeta,*bestAlpha);	//build all the matrices
  if (progressToCerr) {
    timer.stop();
    cerr << "done in " << timer << " sec" 
	 << " chisq=" << bestChisq 
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
      cerr << "## Marquart iteration " << i << " inversion...";
    }

    *beta /= *alpha;

    if (progressToCerr) {
      timer.stop();
      cerr << "done in " << timer << " sec" << endl;
    }

    // Seems necessary for TMV to work:
    if (saveMemory)  alpha->unsetDiv();

    a = *bestA + *beta;

    // Get chisq and derivs at the new spot
    if (progressToCerr) {
      timer.reset();
      timer.start();
      cerr << "## Marquart iteration " << i+1 << " derivatives...";
    }

    derivs(a,chisq,*beta,*alpha);

    if (progressToCerr) {
      timer.stop();
      cerr << "done in " << timer << " sec" 
	   << " lambda=" << lambda 
	   << " chisq=" << chisq 
	   << endl;
    }

#ifdef DEBUG
    cerr << "Marquardt Trial " << i 
	 << " lambda=" << lambda 
	 << " chisq=" << chisq 
	 << " best chisq=" << bestChisq 
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
	cerr << "WARNING: Marquardt lambda=" << lambda << endl;
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
  if (isSmall && lastDropWasSmall) return true;  //two small drops=converge
  lastDropWasSmall=isSmall;
  return false;
}

// This class does what mrqmin() does in NR, namely sums of derivs, etc.,
// over a set of (x,y,sigma) data points. Precision is template arg.
// Create this object with the data and pointer to the derivative function.
template <class P=double>
class MarquardtPointSum {
private:
  const Vector<P>& x;
  const Vector<P>& y;
  const Vector<P>& sig;
  Vector<P> isig;
  void (*func)(const Vector<P>& a,
	       const P xpt,
	       P& ypt,
	       Vector<P>& dpt);
public:
  MarquardtPointSum(const Vector<P>& _x, 
		    const Vector<P>& _y, 
		    const Vector<P>& _sig,
		    void (*f)(const Vector<P>& a,
			      const P xpt,
			      P& ypt,
			      Vector<P>& dpt)):
    x(_x), y(_y), sig(_sig), func(f), isig(x.size()) {
    if (y.size() < x.size() || sig.size() < x.size())
      throw std::runtime_error("Data or Sigma vector size mismatch in "
			       "MarquandtPointSum");
    for (int i=0; i<x.size(); i++) isig[i] = pow(sig[i],-2.);
  }
  void operator() (const Vector<P>& a, 
		   P& chisq, 
		   Vector<P>& beta, 
		   SymMatrix<P>& alpha) const {
    chisq = 0;
    beta.setZero();
    alpha.setZero();
    Vector<P> derivs(a.size());
    P yfit;
    for (int k=0; k<x.size(); k++) {
      func(a, x[k], yfit, derivs);
      yfit = y[k] - yfit;
      chisq += yfit*yfit*isig[k];
      for (int i=0; i<a.size(); i++) {
	beta[i] += yfit*derivs[i]*isig[k];
	for (int j=0; j<=i; j++) 
	  alpha(i,j)+=derivs[i]*derivs[j]*isig[k];
      }
    }
  }
};

#endif   //MARQUARDT_H
