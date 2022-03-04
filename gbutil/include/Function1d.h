// Useful unary functions; binary here also!

// Multiplying Function1d<>'s gives the ProductFunction1d<>, which is also a 
// Function1d<>.  
// Note that the ProductFunction1d requires that its component objects still
// exist.  Don't destroy components before Product.
// The exception is that intermediate Products can go away, so it's ok to
// write Product pp = f1 * f2 * f3, as the temporary object isn't needed.
#ifndef FUNCTION1D_H
#define FUNCTION1D_H

#include <vector>
#include <list>
using std::list;
using std::vector;

#include "Std.h"

#define BIG 1e30
#define BIGNEG -1e30

template<class V, class A>
class Function1d;

template<class V, class A>
class ProductFunction1d;

template<class V=double, class A=double>
class Function1d {
 public:
 virtual V operator()(const A a) const =0;
 virtual A argMin() const =0; //Range of function defined
 virtual A argMax() const =0;
 virtual A supportMin() const {return argMin();} //Range of function non-zero
 virtual A supportMax() const {return argMax();}
 ProductFunction1d<V,A> operator*(Function1d &rhs) const {
   return ProductFunction1d<V,A>(this,&rhs);
 }
};

template<class A1=double, class A2=double, class V=double>
class Function2d {
 public:
 //need this??
 // virtual V operator()(const A1 a1, const A2 a2) const =0;
};

// constant return
class Constant: public Function1d<> {
 public:
  Constant(double _c=1.): c(_c) {}
  double operator()(const double d) const {return c;}
  double argMin() const {return BIGNEG;}
  double argMax() const {return BIG;}
 private:
  double c;
};

// boxcar function - double precision
class Boxcar: public Function1d<> {
 public:
  Boxcar(double _floor, double _ceil, double _v=1.):
    floor(_floor), ceil(_ceil), val(_v) {}
  double operator()(const double a) const {return (a<ceil)*(a>floor)*val;}
  double argMin() const {return BIGNEG;}
  double argMax() const {return BIG;}
  double supportMin() const {return floor;}
  double supportMax() const {return ceil;}
 private:
  double floor;
  double ceil;
  double val;
};

// Function which is a product of other functions
template<class V=double, class A=double>
class ProductFunction1d: public Function1d<V,A> {
public:
  // Create with 2, 1, or 0 initial functions
  ProductFunction1d(const Function1d<V,A>* f1=0, const Function1d<V,A>* f2=0) {
    if (f1) Append(f1); if (f2) Append(f2);
  }
  void Append(const Function1d<V,A>* f) {
    if (flist.size()==0) {
      amin=f->argMin(); amax=f->argMax();
      smin=f->supportMin(); smax=f->supportMax();}
    else {
      amin= (f->argMin()>amin ? f->argMin() : amin);
      amax= (f->argMax()<amax ? f->argMax() : amax);
      smin= (f->supportMin()>smin ? f->supportMin() : smin);
      smax= (f->supportMax()<smax ? f->supportMax() : smax);
    }
    if ( const ProductFunction1d* fp
	 = dynamic_cast< const ProductFunction1d*> (f) ) {
      //Appending what's already a product; copy its list
      flist.insert(flist.end(),fp->flist.begin(),fp->flist.end());
    } else {
      //Appending a single function, add its pointer
      flist.push_back(f);
    }
  }
  V operator()(const A arg) const {
    typename list<const Function1d<V,A>*>::const_iterator i = flist.begin();
    if (i==flist.end()) return 0;
    V retval=(**i)(arg); ++i;
    for ( ; i!=flist.end(); ++i)
      retval *= (**i)(arg);
    return retval;
  }
  void operator*=(const Function1d<V,A> &f) {Append(&f);}
  A argMin() const {return amin;}
  A argMax() const {return amax;}
  A supportMin() const {return smin;}
  A supportMax() const {return smax;}
private:
  list<const Function1d<V,A>*> flist;
  A amin;
  A amax;
  A smin;
  A smax;
};

// General polynomial function
template<class V=double, class A=double>
class Polynomial: public Function1d<V,A> {
private:
  int order;
  vector<V> coeffs;
public:
  Polynomial(int _ord): order(_ord) {
    coeffs.resize(order+1); for (int i=0; i<=order; i++) coeffs[i]=0;
  }
  Polynomial(int _ord, vector<V>& cin): order(_ord) {
    if (cin.size()<order+1) 
      FormatAndThrow<std::runtime_error>() << "Input coefficient array size "
					   << cin.size() << " too small for Polynomial order"
					   << order;
    coeffs.resize(order+1); 
    for (int i=0; i<=order; i++) coeffs[i]=cin[i];
  }
  V operator()(const A arg) const {
    V sum=0; 
    for (int i=order; i>=0; i--) sum=sum*arg+coeffs[i];
    return sum;
  }
  double argMin() const {return BIGNEG;}
  double argMax() const {return BIG;}
  int NParams() const {return order+1;}
  double partial(const A arg, int iparam) const {
    if (iparam<0 || iparam>order)
      FormatAndThrow<std::runtime_error>() << "Invalid parameter number " << iparam
					   << " in Polynomial";
    return pow(arg, static_cast<double> (iparam));
  }
  void Set(const V c, int iparam) {
    if (iparam<0 || iparam>order)
      FormatAndThrow<std::runtime_error>() << "Invalid parameter number " << iparam
					   << " in Polynomial.Set()";
    coeffs[iparam]=c;
  }
  V Get(int iparam) const {
    if (iparam<0 || iparam>order)
      FormatAndThrow<std::runtime_error>() << "Invalid parameter number " << iparam
					   << " in Polynomial.Get()";
    return coeffs[iparam];
  }
};

#undef BIG
#undef BIGNEG

#endif
