// Classes to represent lookup tables.
// A is the argument class, which must have ordering
// operations, and +-*/ for interpolation.
// D is the value class, which must have + and * operations
// to permit interpolation.
//
// Behavior with duplicate arguments is not well defined.  Current
// implementation with std::set will probably ignore beyond 1st one.
#ifndef TABLE_H
#define TABLE_H

#include "Function1d.h"
#include "Std.h"

#include <vector>
#include <set>
using std::vector;
#include <algorithm>
#include <string>
#include <sstream>
#include <stdexcept>

// Exception classes:
class TableError: public std::runtime_error {
 public:
  TableError(const string &m=""): std::runtime_error("Table Error: " +m) {}
};
class TableOutOfRange: public TableError {
public:
  TableOutOfRange(): TableError("Argument out of range") {}
};
class TableReadError: public TableError {
public:
  TableReadError(const string &c): 
    TableError("Data read error for line ->"+c) {}
};

// Table element:
template<class V=double, class A=double>
class TableEntry {
public:
  TableEntry(A a, V v): arg(a), val(v) {}
  A arg;
  V val;
  bool operator<(const TableEntry rhs) const {return arg<rhs.arg;}

  // Cache some things to speed up interpolations
  mutable A prevArg;  // arg of next entry in sorted table.
  mutable A inverseInterval;  // 1/(arg - prevArg)
  mutable V prevVal;   // value at next arg
  mutable V y2Prev;    // y2[prev element] * interval^2 / 6., used for spline.
  mutable V y2This;    // y2[this element] * interval^2 / 6., used for spline.
  bool contains(const A rhs) const {return rhs>=prevArg && rhs<=arg;}
};

// The Table itself:
template<class V=double, class A=double>
class Table: public Function1d<V,A> {
public:
  enum interpolant {linear, spline, floor, ceil};
  //Construct empty table
  Table(interpolant i=linear): v(), iType(i), isReady(false) {} 
  //Table from two arrays:
  Table(const A* argvec, const V* valvec, int N, interpolant in=linear) ;
  Table(const vector<A> &a, const vector<V> &v, interpolant in=linear) ;
  Table(istream &is, interpolant in=linear): v(), iType(in), isReady(false)
					     {read(is);}
  void clear() {v.clear(); isReady=false;}
  void read(istream &is);
  void addEntry(const A a, const V v) ; //new element for table.
  V operator() (const A a) const ;	//lookup & interp. function value.
  V lookup(const A a) const ;	//interp, but exception if beyond bounds
  int size() const {return v.size();}	//size of table
  A argMin() const { setup(); if (v.size()>0) return _argMin;
   else throw TableError("argMin for null Table");}	//Smallest argument
  A argMax() const { setup(); if (v.size()>0) return _argMax;
   else throw TableError("argMax for null Table");} 	//Largest argument

  template <class T>
  void TransformVal(T &xfrm) {
    for (iter p=v.begin(); p!=v.end(); ++p)
      p->val = xfrm(p->arg, p->val);
    isReady=false; setup();
  }

/**/void dump() const {for (citer p=v.begin(); p!=v.end(); ++p) 
		 cout << p->arg << " " << p->val << endl; }
private:
  typedef TableEntry<V,A> Entry;
  typedef typename set<Entry>::const_iterator citer;
  typedef typename set<Entry>::iterator iter;
  const interpolant iType;

  mutable std::set<Entry> v;
  mutable A _argMin;
  mutable A _argMax;
  mutable citer	lastEntry;	//Element used for last lookup into table.
  mutable bool  isReady;	//Flag if table has been prepped.
  mutable bool  equalSpaced;	//Flag set if arguments are nearly equally spaced.
  mutable vector<citer> directIndex;  // Direct access to equal-spaced elements
  mutable A     dx;		// ...in which case this is argument interval
  //get index to 1st element >= argument.  Can throw the exception here.
  citer upperIndex(const A a) const;
  void setup() const;	//Do any necessary preparation;
  //Interpolate value btwn p & --p:
  V interpolate(const A a, const citer p) const; 
};
#endif
