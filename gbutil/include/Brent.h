// Template to find the minimum of a function.
// Use Numerical Recipes routines

#ifndef BRENT_H
#define BRENT_H

#include "Std.h"
#include <cmath>

namespace brent {

  class BrentError: public std::runtime_error {
  public:
    BrentError(const string m): std::runtime_error("Brent error: "+m) {}
  };

  const double defaultTolerance=1.e-7;
  const int defaultMaxSteps=40;

  template <class F, class T=double>
  class Brent {
  private:
    T	ax, bx, cx;
    T   fa, fb, fc;
    T	xTolerance;
    int	maxSteps;
    F&  func;
    bool isBracketed;
    static const double defaultTolerance;
    static const int defaultMaxSteps=40;
    void shift(T& a, T& b, T& c, T d) {a=b; b=c; c=d;}
    T sign(T a, T b) {return b>=0 ? abs(a) : -abs(a);}

  public:
    // Constructor where minimum is not known to be bracketed
    Brent(F& func_, T ax_=1., T bx_=0.):
      func(func_), ax(ax_), bx(bx_), xTolerance(defaultTolerance),
      maxSteps(defaultMaxSteps), isBracketed(false) {}
    // Constructor starting with a triple means that we know
    // there's a minimum between ax and cd
    Brent(F& func_, T ax_, T bx_, T cx_,
	  T fa_, T fb_, T fc_):
      func(func_), ax(ax_), bx(bx_), cx(cx_), 
      fa(fa_), fb(fb_), fc(fc_),
      xTolerance(defaultTolerance),
      maxSteps(defaultMaxSteps), isBracketed(true) {}
    void setMaxSteps(int m) {maxSteps=m;}
    T    getXTolerance() const {return xTolerance;}
    void setXTolerance(T tol) {xTolerance=tol;}
    // Search for bracket of minimum, starting with this range:
    void bracket(T ax, T bx);
    // Return location and value of minimum:
    T minimum(T& xmin);
  };

  template<class F, class T>
  const double 
  Brent<F,T>::defaultTolerance=1e-7;

  template<class F, class T>
  void
  Brent<F,T>::bracket(T ax_, T bx_) {
    // This is mnbrak.c from Numerical Recipes
    const double GOLD=1.618034;
    const double BIGGEST_STEP_RATIO=5.;	// How far to follow parabola
    const double SMALL = 1e-20;	// Worry about roundoff in double here
    ax = ax_;
    bx = bx_;
    isBracketed = false;
    T u,fu;
    
    fa=func(ax);
    fb=func(bx);
    int stepCount = 2;	// count function evalutions
    if (fb > fa) {
      SWAP(ax,bx);
      SWAP(fa,fb);
    }
    cx=bx + GOLD*(bx-ax);
    fc=func(cx);
    while (fb > fc) {
      if (stepCount > maxSteps) {
	isBracketed = false;
	return;
      }
      // Do parabolic solution
      double r=(bx-ax)*(fb-fc);
      double q=(bx-cx)*(fb-fa);
      double step=-1.;
      if ( abs(q-r)>SMALL) {
	step = 0.5*((ax-bx)*r-(cx-bx)*q) / (q-r) / (cx-bx);
	if (step > BIGGEST_STEP_RATIO) step = BIGGEST_STEP_RATIO;
      }
      double u=bx + (cx-bx)*step;
      if (step<=0.) {
	// parabola goes backwards - take golden step beyond cx
	u=cx+GOLD*(cx-bx);
	fu=func(u);
	stepCount++;
      } else if (step < 1.) {
	// Try point between b and c:
	fu=func(u);
	stepCount++;
	if (fu < fc) {
	  ax=bx;
	  bx=u;
	  fa=fb;
	  fb=fu;
	  break;
	} else if (fu > fb) {
	  cx=u;
	  fc=fu;
	  break;
	}
	// fu was between fb and fc: keep going with golden:
	u=cx+GOLD*(cx-bx);
	fu=func(u);
	stepCount++;
      } else {
	// u is beyond cx:
	fu=func(u);
	stepCount++;
	if (fu < fc) {
	  // take another golden step
	  shift(bx,cx,u,cx+GOLD*(cx-bx));
	  shift(fb,fc,fu,func(u));
	  stepCount++;
	}
      }
      // Throw away ax:
      shift(ax,bx,cx,u);
      shift(fa,fb,fc,fu);
    }
    // Success here
    isBracketed = true;
  };

  template<class F, class T>
  T
  Brent<F,T>::minimum(T& xmin) {
    // Brent's method, from Numberical Recipes
    int iter;
    T a,b,d,etemp,fu,fv,fw,fx,p,q,r,tol2,u,v,w,x,xm;
    T e=0.0;
    const double ZEPS=1e-10;	// sqrt of small step
    const double CGOLD=0.3819660;
    if (!isBracketed) bracket(ax, bx);
    if (!isBracketed || fb>fa || fb>fc)
      FormatAndThrow<BrentError>() << "minimum() inputs " << ax
				   << ", " << bx << ", " << cx
				   << " do not bracket min";
    a=(ax < cx ? ax : cx);
    b=(ax > cx ? ax : cx);
    x=w=v=bx;
    fw=fv=fx=fb;
    /*cerr << "a b c " << ax << " " << bx << " " << cx
      << " fa fb fc " << fa << " " << fb << " " << fc << endl; /**/

    for (iter=0;iter<maxSteps;iter++) {
      /*cerr << "  iter " << iter 
	       << " u v w x " << u << " " << v << " " << w << " " << x
	       << " fx " << fx
	       << endl; /**/
      xm=0.5*(a+b);
      tol2=2.0*xTolerance;
      if (abs(x-xm) <= (tol2-0.5*(b-a))) {
	xmin=x;
	return fx;
      }
      if (abs(e) > xTolerance) {
	r=(x-w)*(fx-fv);
	q=(x-v)*(fx-fw);
	p=(x-v)*q-(x-w)*r;
	q=2.0*(q-r);
	if (q > 0.0) p = -p;
	q=abs(q);
	etemp=e;
	e=d;
	if (abs(p) >= abs(0.5*q*etemp) || p <= q*(a-x) || p >= q*(b-x))
	  d=CGOLD*(e=(x >= xm ? a-x : b-x));
	else {
	  d=p/q;
	  u=x+d;
	  if (u-a < tol2 || b-u < tol2)
	    d=sign(xTolerance,xm-x);
	}
      } else {
	d=CGOLD*(e=(x >= xm ? a-x : b-x));
      }
      u=(abs(d) >= xTolerance ? x+d : x+sign(xTolerance,d));
      fu=func(u);
      if (fu <= fx) {
	if (u >= x) a=x; else b=x;
	shift(v,w,x,u);
	shift(fv,fw,fx,fu);
      } else {
	if (u < x) a=u; else b=u;
	if (fu <= fw || w == x) {
	  v=w;
	  w=u;
	  fv=fw;
	  fw=fu;
	} else if (fu <= fv || v == x || v == w) {
	  v=u;
	  fv=fu;
	}
      }
    }
    FormatAndThrow<BrentError>() << "Too many iterations " << iter;
    xmin=x;
    return fx;
  }

} // namespace brent
#endif
