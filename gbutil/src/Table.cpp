#include "Table.h"
#include <cmath>

// Look up an index.  Use STL binary search; maybe faster to use
template<class V, class A>
typename Table<V,A>::citer Table<V,A>::upperIndex(const A a) const {
  setup();
  if (v.size()==0 || a<_argMin)  throw TableOutOfRange();
  // Go directly to index if arguments are regularly spaced.
  if (equalSpaced) {
    int index = static_cast<int> ( std::ceil( (a-_argMin) / dx) );
    if (index==0) return v.begin();
    if (index >= v.size()) throw TableOutOfRange();
    // check if we need to move ahead or back one step due to rounding errors
    citer p = directIndex[index];
    if (a > p->arg) { 
      ++p;
      if (p == v.end()) throw TableOutOfRange();
    } else if (a < p->prevArg) {
      --p;
    }
    return p;
  }


  // ??? might choose to skip this!
  // First see if the previous index is still ok
  if (lastEntry!=v.end() && lastEntry->contains(a)) return lastEntry;

  // This STL algorithm uses binary search to get 1st element >= ours.
  Entry e(a,0); 
  citer p = v.lower_bound(e);
  // bounds check
  if (p==v.end()) throw TableOutOfRange();
  lastEntry = p;
  return p;
}

//new element for table.
template<class V, class A>
void Table<V,A>::addEntry(const A _arg, const V _val) {
  v.insert(Entry(_arg,_val));
  isReady = false;
}

template<class V, class A>
Table<V,A>::Table(const A* argvec, const V* valvec, int N, 
		  interpolant in): v(), iType(in), isReady(false) {
  const A* aptr;
  const V* vptr;
  int i;
  for (i=0, aptr=argvec, vptr=valvec; i<N; i++, aptr++, vptr++) {
    v.insert(Entry(*aptr,*vptr));
  }
}

template<class V, class A>
Table<V,A>::Table(const vector<A> &aa, const vector<V> &vv, 
		  interpolant in): v(), iType(in), isReady(false) {
  if (vv.size()<aa.size()) 
    throw TableError("input vector lengths don't match");
  typename vector<A>::const_iterator aptr=aa.begin();
  typename vector<V>::const_iterator vptr=vv.begin();
  for (int i=0; i<aa.size(); i++, ++aptr, ++vptr) {
    v.insert(Entry(*aptr,*vptr));
  }
  isReady = false;
}

//lookup & interp. function value. - this one returns 0 out of bounds.
template<class V, class A>
V Table<V,A>::operator() (const A a) const {
  try {
    citer p1(upperIndex(a));
    return interpolate(a,p1);
  } catch (TableOutOfRange) {
    return static_cast<V> (0);
  }
}
//lookup & interp. function value.  Throws if out of bounds.
template<class V, class A>
V Table<V,A>::lookup(const A a) const {
  citer p1(upperIndex(a));
  return interpolate(a,p1);
}

template<class V, class A>
V Table<V,A>::interpolate(const A a, const citer p1) const { 
  setup();	//do any necessary prep
  // First case when there is for single-point table
  if (v.size()==1) return p1->val;
  else if (iType==linear) {
    if (p1==v.begin())  return p1->val;
    double frac=(a - p1->prevArg) * p1->inverseInterval;
    return frac*p1->val + (1-frac) * p1->prevVal;
  } else if (iType==spline) {
    if (p1==v.begin())  return p1->val;
    A aa=(p1->arg - a) * p1->inverseInterval;
    A bb=1.-aa;
    return aa*p1->prevVal + bb*p1->val +
      (aa*aa*aa-aa)*p1->y2Prev + (bb*bb*bb-bb)*p1->y2This;
  } else if (iType==floor) {
    if (p1==v.begin())  return p1->val;
    else return p1->prevVal;
  } else if (iType==ceil) {
    return p1->val;
  } else {
    throw TableError("interpolation method not yet implemented");
  }
}

template<class V, class A>
void Table<V,A>::read(istream &is) {
  string line;
  const string comments="#;!";	//starts comment
  V vv;
  A aa;
  while (is) {
    getline(is,line);
    // skip leading white space:
    int i;
    for (i=0;  isspace(line[i]) && i<line.length(); i++) ;
    // skip line if blank or just comment
    if (i==line.length()) continue;
    if (comments.find(line[i])!=string::npos) continue;
    // try reading arg & val from line:
    std::istringstream iss(line);
    iss >> aa >> vv;
    if (iss.fail()) throw TableReadError(line) ;
    addEntry(aa,vv);
  }
}

// Do any necessary setup of the table before using
template<class V, class A>
void Table<V,A>::setup() const {
  if (isReady) return;
  // Invalidate lastEntry since container may have changed
  lastEntry = v.end();
  equalSpaced = false;
  if (v.size()==0) return;  // Nothing to do without entries.

  const int n = v.size();

  _argMin = v.begin()->arg;
  {
    citer p = v.end();
    --p;
    _argMax = p->arg;
  }

  if (v.size()<=1) {
    // Nothing to do if the table is degenerate
    isReady = true;
    return;
  }

  // See if arguments are equally spaced
  // ...within this fractional error:
  const double tolerance = 0.01;
  dx = (_argMax - _argMin) / (v.size()-1);
  equalSpaced = true;
  {
    citer p = v.begin();
    ++p;
    for (int i=1; i<v.size(); i++, ++p) {
      if ( abs( (p->arg-_argMin)/dx - i) > tolerance) {
	equalSpaced = false;
	break;
      }
    }
  }
  if (equalSpaced) {
    // Set up indices to entries
    directIndex.clear();
    for (citer p = v.begin(); p != v.end(); ++p)
      directIndex.push_back(p);
  }

  vector<V> y2(n);
  if (iType==spline) {
    // Set up the 2nd-derivative table for splines
    // Derive this from Numerical Recipes spline.c
    if (n<2) throw TableError("Spline Table with only 1 entry");

    V  p,qn,sig,un;

    vector<V> u(n-1);
    vector<A> args(n);
    vector<V> vals(n);
    y2[0]=u[0]= static_cast<V>(0);
    citer p1 = v.begin();
    for (int i=0; i<n; i++, ++p1) {
      args[i] = p1->arg;
      vals[i] = p1->val;
    }

    for (int i=1;i<=n-2;i++) {
      sig=(args[i]-args[i-1])/(args[i+1] - args[i-1]);
      p=sig*y2[i-1]+2.0;
      y2[i]=(sig-1.0)/p;
      u[i]=(vals[i+1]-vals[i])/(args[i+1]-args[i]) 
	- (vals[i]-vals[i-1])/(args[i]-args[i-1]);
      u[i]=(6.0*u[i]/(args[i+1]-args[i-1])-sig*u[i-1])/p;
    }
	
    qn=un=0.0;

    y2[n-1]=(un-qn*u[n-2])/(qn*y2[n-2]+1.0);
    for (int k=n-2;k>=0;k--)
      y2[k]=y2[k]*y2[k+1]+u[k];

  }

  iter p0 = v.begin();
  iter p1 = v.begin();
  ++p1;
  for (int i=1; i<n; i++, ++p0, ++p1) {
    p1->prevArg = p0->arg;
    A h = p1->arg - p1->prevArg;
    p1->inverseInterval = 1./h;
    p1->prevVal = p0->val;
    if (iType==spline) {
      p1->y2Prev = y2[i-1]*h*h/6.;
      p1->y2This = y2[i]*h*h/6.;
    }
  }
  isReady = true;
  return;
}

template class Table<>;
