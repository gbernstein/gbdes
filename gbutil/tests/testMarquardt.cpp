// Test Marquardt by fitting to a polynomial.
// 

#include <iostream>
#include "LinearAlgebra.h"
#include "Marquardt.h"
#include "Random.h"


typedef linalg::Vector<double> Vector;

void poly(const Vector& a,
	  double x,
	  double& y,
	  Vector& derivs) {
  y = 0.;
  double xp = 1.;
  for (int i=0; i<a.size(); i++) {
    y += xp * a[i];
    derivs[i] = xp;
    xp *= x;
  }
  return;
}

int main(int argc,
	 char *argv[])
{
  int npts = 100;
  double dx = 0.1;
  double noise = 0.1;
  Vector a(3);
  a[0] = 5.; a[1] = 3.; a[2]= -8.;
  Vector x(100);
  Vector y(100);
  Vector sigma(100,noise);
  ran::GaussianDeviate<> gd(1224);

  for (int i=0; i<npts; i++) {
    double xx = dx * i;
    x[i] = xx;
    double xp = 1.;
    double yy = 0.;
    for (int j=0; j<a.size(); j++) {
      yy += a[j] * xp;
      xp *= xx;
    }
    y[i] = yy + noise*gd;
    sigma[i] = noise;
  }

  MarquardtPointSum<> mps(x,y,sigma,&poly);
  Vector atry(a * 0.5);
  Marquardt<> marq(mps);
  marq.setSaveMemory();  /*??? Try both with and without this*/ 
  try {
    double chisq = marq.fit(atry, 50, true);
    cout << "Chisq: " << chisq << endl;
    cout << "Params: " << atry << endl;
    linalg::Matrix<double> alpha = marq.getAlpha();
    cout << "Alpha: " << alpha << endl;
    auto covar = alpha.inverse();
    cout << "sigmas: ";
    for (int i=0; i<alpha.nrows(); i++) cout << " " << sqrt(covar(i,i));
    cout << endl;
  } catch (std::runtime_error& m) {
    cerr << m.what() << endl;
    exit(1);
  }

  exit(0);
}

      

  
