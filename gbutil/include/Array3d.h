// Cubic 3-dimensional array template.  Just simple scalar
// operations and element access
//---------------------------------------------------------------------------
#ifndef Array3dH
#define Array3dH

#include "Std.h"
#include <vector>
using std::vector;

template <class T> class Cube {

public:
  Cube() : m(0),va(0) {}
  Cube(size_t mm, T val=0) 
    : m(mm),va(m*m*m, val) {}
  Cube(const Cube& rhs)
    : m(rhs.m),va(rhs.va) {}
  Cube<T>& operator=(const Cube<T>& rhs)
    { if (&rhs == this) return *this;
      m=rhs.m; va = rhs.m;
      return *this; }
  Cube<T>& operator=(const T rhs) {
    for (int i=0; i<va.size(); i++) va[i]=rhs; 
    return *this;}
  ~Cube() {}
  void resize(size_t mm) 
    { m = mm; va.resize(m*m*m);}
  size_t getM() const {return m;}

  T& operator()(size_t i,size_t j, size_t k) 
    { Assert(i<m && j<m && k<m);
      return va[index(i,j,k)]; }
  T operator()(size_t i,size_t j, size_t k) const
    { Assert(i<m && j<m && k<m);
      return va[index(i,j,k)]; }
  void Zero() {*this=0.;}

  Cube& operator*=(T x) {
    for (int i=0; i<va.size(); i++) va[i]*=x;
    return *this; }
  Cube& operator/=(T x) {
    for (int i=0; i<va.size(); i++) va[i]/=x;
    return *this; }

 private:
  size_t m;
  vector<T> va;
  size_t index(size_t i,size_t j, size_t k) const { return (i*m+j)*m+k; }

}; // Cube

#endif
