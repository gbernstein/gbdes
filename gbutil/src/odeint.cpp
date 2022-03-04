#include "odeint.h"
#include "Std.h"
#include <stdexcept>

using namespace ode;
using namespace linalg;

void
ode::rkck(const DVector& y,
	  const DVector& dydx,
	  double x,
	  double h,
	  DVector& yout,
	  DVector& yerr,
	  const ODEFunction& derivs)  {

  int n=y.size();
  Assert(n==dydx.size());
  Assert(n==yout.size());
  Assert(n==yerr.size());
  Assert(n==derivs.order());

  static double a2=0.2,a3=0.3,a4=0.6,a5=1.0,a6=0.875,b21=0.2,
		b31=3.0/40.0,b32=9.0/40.0,b41=0.3,b42 = -0.9,b43=1.2,
		b51 = -11.0/54.0, b52=2.5,b53 = -70.0/27.0,b54=35.0/27.0,
		b61=1631.0/55296.0,b62=175.0/512.0,b63=575.0/13824.0,
		b64=44275.0/110592.0,b65=253.0/4096.0,c1=37.0/378.0,
		c3=250.0/621.0,c4=125.0/594.0,c6=512.0/1771.0,
		dc5 = -277.0/14336.0;
  static double dc1=c1-2825.0/27648.0,dc3=c3-18575.0/48384.0,
                dc4=c4-13525.0/55296.0,dc6=c6-0.25;
  //  static float *ak2,*ak3,*ak4,*ak5,*ak6,*ytemp;
  DVector ytemp(n);
  for (int i=0; i<n; i++)
    ytemp[i]=y[i]+b21*h*dydx[i];
  DVector ak2 = derivs(x+a2*h,ytemp);
  for (int i=0; i<n; i++)
    ytemp[i]=y[i]+h*(b31*dydx[i]+b32*ak2[i]);
  DVector ak3 = derivs(x+a3*h,ytemp);
  for (int i=0; i<n; i++)
    ytemp[i]=y[i]+h*(b41*dydx[i]+b42*ak2[i]+b43*ak3[i]);
  DVector ak4 = derivs(x+a4*h,ytemp);
  for (int i=0; i<n; i++)
    ytemp[i]=y[i]+h*(b51*dydx[i]+b52*ak2[i]+b53*ak3[i]+b54*ak4[i]);
  DVector ak5 = derivs(x+a5*h,ytemp);
  for (int i=0; i<n; i++)
    ytemp[i]=y[i]+h*(b61*dydx[i]+b62*ak2[i]+b63*ak3[i]+b64*ak4[i]+b65*ak5[i]);
  DVector ak6 = derivs(x+a6*h,ytemp);
  for (int i=0; i<n; i++)
    yout[i]=y[i]+h*(c1*dydx[i]+c3*ak3[i]+c4*ak4[i]+c6*ak6[i]);
  for (int i=0; i<n; i++)
    yerr[i]=h*(dc1*dydx[i]+dc3*ak3[i]+dc4*ak4[i]+dc5*ak5[i]+dc6*ak6[i]);
}

#define SAFETY 0.9
#define PGROW -0.2
#define PSHRNK -0.25
#define ERRCON 1.89e-4

void 
ode::rkqs(DVector& y, 
	  DVector& dydx,
	  double& x,
	  double htry,
	  double eps,
	  const DVector yscal,
	  double& hdid,
	  double& hnext,
	  const ODEFunction& derivs)
{
  int n=derivs.order();
  Assert(n==y.size());
  Assert(n==dydx.size());
  Assert(n==yscal.size());
  DVector yerr(n);
  DVector ytemp(n);

  double h=htry;
  double errmax, xnew, htemp;
  for (;;) {
    /*cerr << "Step from " << x << " by " << h 
      << " values " << y[0] << " " << y[1] << endl;*/
    rkck(y,dydx,x,h,ytemp,yerr,derivs);
    errmax=0.0;
    for (int i=0; i<n; i++) errmax=MAX(errmax,fabs(yerr[i]/yscal[i]));
    errmax /= eps;
    if (errmax > 1.0) {
      htemp=SAFETY*h*pow(errmax,PSHRNK);
      h = (h>=0.0 ? MAX(htemp, 0.1*h) : MIN(htemp, 0.1*h));
      xnew=x+h;
      if (xnew == x) 
	throw std::runtime_error("stepsize underflow in rkqs");
      continue;
    } else {
      if (errmax > ERRCON) hnext=SAFETY*h*pow(errmax,PGROW);
      else hnext=5.0*h;
      x += (hdid=h);
      for (int i=0; i<n; i++) y[i]=ytemp[i];
      break;
    }
  }
}
#undef SAFETY
#undef PGROW
#undef PSHRNK
#undef ERRCON

/////////////////////////////////////////////////////////////////

#define MAXSTP 10000
#define TINY 1.0e-30

void 
ode::odeint(DVector& ystart,
	    DVector& xp,
	    DMatrix& yp,
	    double dxsav,
	    int& kount,
	    double x1,
	    double x2,
	    double eps,
	    double h1,
	    double hmin,
	    int& nok,
	    int& nbad,
	    const ODEFunction& derivs)
{
  int nvar = derivs.order();
  int kmax = xp.size();
  Assert (kmax==0 || nvar == yp.nrows());
  Assert (kmax <= yp.ncols());

  double xsav,hnext,hdid;

  Assert(nvar==ystart.size());

  DVector yscal(nvar);
  DVector y(nvar);

  double x=x1;
  double h= x2-x1>=0. ? abs(h1) : -abs(h1);
  kount=0;
  nok = nbad = 0;
  for (int i=0;i<nvar;i++) y[i]=ystart[i];
  if (kmax > 0) xsav=x-dxsav*2.0;
  for (int nstp=0;nstp<MAXSTP;nstp++) {
    DVector dydx = derivs(x,y);
    for (int i=0;i<nvar;i++)
      yscal[i]=abs(y[i])+abs(dydx[i]*h)+TINY;
    if (kmax > 0 && kount < kmax-1 && abs(x-xsav) > abs(dxsav)) {
      xp[kount]=x;
      for (int i=0;i<nvar;i++) yp(i,kount)=y[i];
      kount++;
      xsav=x;
    }
    if ((x+h-x2)*(x+h-x1) > 0.0) h=x2-x;
    rkqs(y,dydx,x,h,eps,yscal,hdid,hnext,derivs);
    if (hdid == h) ++nok; else ++nbad;
    if ((x-x2)*(x2-x1) >= 0.0) {
      for (int i=0;i<nvar;i++) ystart[i]=y[i];
      if (kmax) {
	xp[kount]=x;
	for (int i=0;i<nvar;i++) yp(i,kount)=y[i];
	kount++;
      }
      return;
    }
    if (abs(hnext) <= hmin) 
      throw std::runtime_error("Step size too small in odeint");
    h=hnext;
  }
  throw std::runtime_error("Too many steps in routine odeint");
}
#undef MAXSTP
#undef TINY
