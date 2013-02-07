// $Id: Poly2d.cpp,v 1.4 2012/02/02 02:20:04 garyb Exp $
#include "Poly2d.h"
#include "StringStuff.h"

using namespace poly2d;


/////////////////////////////////////////////////////////////////////
// Serialization
/////////////////////////////////////////////////////////////////////
void
Poly2d::write(std::ostream& os) const {
  switch (otype) {
  case Sum:
    os << "Sum " << orderX << endl;
    break;
  case Each:
    os << "Each " << orderX << " " << orderY << endl;
    break;
  default:
    throw Poly2dError("Bad Poly2d OrderType in Poly2d::write()");
  }
  const int coeffsPerLine = 5;
  int iOnLine = 0;
  DVector coeffs = getC();
  for (int i=0; i<nCoeffs(); i++) {
    if (iOnLine >= coeffsPerLine) {
      os << endl;
      iOnLine = 0;
    }
    os << coeffs[i] << " " ;
    iOnLine++;
  }
  os << endl;
}

Poly2d*
Poly2d::create(std::istream& is) {
  string buffer;
  if (!stringstuff::getlineNoComment(is, buffer))
    throw Poly2dError("Missing data reading Poly2d::create()");
  string otypeString;
  bool orderTypeEach = false; // false means otype=Sum
  int ordX;
  int ordY;
  {
    bool badFormat=false;
    istringstream iss(buffer);
    if (!(iss >> otypeString)) {
      badFormat = true;
    } else if (stringstuff::nocaseEqual(otypeString,"Sum")) {
      if (! (iss >> ordX))
	badFormat = true;
    } else if (stringstuff::nocaseEqual(otypeString,"Each")) {
      if (! (iss >> ordX >> ordY))
	badFormat = true;
      else 
	orderTypeEach=true;
    } else {
      badFormat = true;
    }
    if (badFormat)
      throw Poly2dError("Bad order format for Poly2d::create: " + buffer);
  }
  Poly2d* poly = orderTypeEach ? new Poly2d(ordX, ordY) : new Poly2d(ordX);

  // Read coefficients:
  int iCoeff=0;
  DVector coeffs(poly->nCoeffs());
  while (iCoeff < coeffs.size()) {
    if (!(stringstuff::getlineNoComment(is, buffer)))
      throw Poly2dError("Out of input reading coefficients for Poly2d: " + buffer);
    double c;
    istringstream iss(buffer);
    do {
      iss >> c;
      if (!iss.fail()) {
	// Successful read
	coeffs[iCoeff] = c;
	iCoeff++;
      } else if (iss.eof()) {
	break; // just ran out of data on this line, get more
      } else {
	// Fail without eof means a true read error:
	throw Poly2dError("Format error reading coefficientss for Poly2d: " + buffer);
      }
    } while (iCoeff < coeffs.size());
  }
  poly->setC(coeffs);
  return poly;
}

/////////////////////////////////////////////////////////////////////
// Computations:
/////////////////////////////////////////////////////////////////////

DVector
Poly2d::powers(double z, int order) const {
  DVector p(order+1, 1.);
  for (int i=1; i<=order; i++)
    p[i] = z*p[i-1];
  return p;
}

DVector
Poly2d::derivs(double z, int order) const {
  DVector d(order+1, 0.);
  double zn=1.;
  for (int i=1; i<=order; i++) {
    d[i] = i*zn;
    zn *= z;
  }
  return d;
}

Poly2d::Poly2d(DMatrix m, OrderType otype_):
  orderX(m.nrows()-1), orderY(m.ncols()-1), otype(otype_), cm(m) {
  // If otype==Sum, check for square and zero the too-high-order elements
  if (otype==Sum) {
    if (orderX != orderY) throw Poly2dError("Constructor matrix not square.");
    for (int i=1; i<=orderX; i++)
      for (int j=orderY-i+1; j<=orderY; j++)
	cm(i,j)=0.;
  }
}

int
Poly2d::vectorIndex(int i, int j) const {
  if (i<0 || j<0) return -1;
  switch (otype) {
  case Sum: 
    // For OrderType=Sum, increase total
    // polynomial order, start each order with x^N:
    // 1, x, y, x^2, xy, y^2, x^3, xy^2, ...
    return ((i+j)*(i+j+1))/2+j;
  case Each: 
    // Standard order: for OrderType=Each, index of 
    // y will be rapidly varying, move along matrix rows:
    // 1., y, y^2, y^3 ... , x, xy, xy^2, ...
    return i*(orderY+1) + j;
  default: 
    throw Poly2dError("Bad OrderType in vectorIndex");
  }
}

void
Poly2d::powersOfIndex(int n, int& i, int& j) const {
  if (n<0) { i=j=-1; return;}
#ifdef _OPENMP 
  int N=0;
  switch (otype) {
  case Sum: 
    while (n > N) {
      ++N;
      n -= N;
    }
    i = N-n;
    j = n;
#else
  // Easiest to just keep a static array for the Sum case:
  // not thread safe!!
  static vector<int> ix(1,0);
  static vector<int> iy(1,0);

  switch (otype) {
  case Sum: 
    while (n >= ix.size()) {
      // Grow the lookup table as needed
      if (ix.back() > 0) {
	ix.push_back( ix.back()-1 );
	iy.push_back( iy.back()+1 );
      } else {
	ix.push_back( iy.back()+1 );
	iy.push_back(0);
      }
    }
    // Lookup values for this n.
    i = ix[n];
    j = iy[n];
#endif

    return;
  case Each:
    i = n / (orderY+1);
    j = n % (orderY+1);
    return;
  default: 
    throw Poly2dError("Bad OrderType in powersOfIndex");
  }
}

DVector
Poly2d::vectorFromMatrix(const DMatrix& m) const {
  Assert(m.nrows()==orderX+1 && m.ncols()==orderY+1);
  switch (otype) {
  case Sum:
    Assert(orderX==orderY);
    {
      int N = (orderX+1)*(orderX+2)/2;
      DVector v(N,0.);
      for (int i=0; i<=orderX; i++)
	for (int j=0; (i+j)<=orderX; j++)
	  v[vectorIndex(i,j)]=m(i,j);
      return v;
    }
  case Each:
    {
      int N = (orderX+1)*(orderY+1);
      DVector v(N,0.);
      for (int i=0; i<=orderX; i++)
	for (int j=0; j<=orderY; j++)
	  v[vectorIndex(i,j)]=m(i,j);
      return v;
    }
  default:
    throw Poly2dError("Bad OrderType in vectorFromMatrix");
  }
}

int
Poly2d::nCoeffs() const {
  switch (otype) {
  case Sum:
    Assert(orderX==orderY);
    return (orderX+1)*(orderX+2)/2;
  case Each:
    return (orderX+1)*(orderY+1);
  default:
    throw Poly2dError("Bad OrderType in nCoeffs()");
  }
}

void
Poly2d::fillFromVector(const DVector& v) {
  if (otype==Sum) Assert(v.size()==(orderX+1)*(orderX+2)/2);
  if (otype==Each) Assert(v.size()==(orderX+1)*(orderY+1));
  int i,j;
  for (int n=0; n<v.size(); n++) {
    powersOfIndex(n,i,j);
    cm(i,j)=v[n];
  }
}

