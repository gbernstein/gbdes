// Code to fit a PixelMap to list of detections using Marquardt-Levenberg class.

#include "Match.h"
#include "Marquardt.h"

using namespace astrometry;

class MapMarq {
public:
  MapMarq(const list<Detection*>& dets_, PixelMap* pm_): dets(dets_), pm(pm_), sigma(1.) {}
  void operator()(const DVector& a, double& chisq, 
		  DVector& beta, tmv::SymMatrix<double>& alpha);
  void setSigma(double s) {sigma = s;}
private:
  const list<Detection*>& dets;
  PixelMap* pm;
  double sigma;
};

void
astrometry::mapFit(list<Detection*> testPoints, PixelMap* pm, double sigma) {
  MapMarq mm(testPoints, pm);
  mm.setSigma(sigma);
  Marquardt<MapMarq> marq(mm);
  const double RELATIVE_TOLERANCE=0.001;
  marq.setRelTolerance(RELATIVE_TOLERANCE);
  DVector p = pm->getParams();
  marq.fit(p);
  pm->setParams(p);
}

void
MapMarq::operator()(const DVector& a, 
		    double& chisq, 
		    DVector& beta, 
		    tmv::SymMatrix<double>& alpha) {
  chisq = 0;
  beta.setZero();
  alpha.setZero();
  int np = a.size();
  Assert(pm->nParams()==np);
  Assert(beta.size()==np);
  Assert(alpha.nrows()==np);
  DMatrix derivs(2,np);
  pm->setParams(a);

  for (list<Detection*>::const_iterator i = dets.begin();
       i != dets.end();
       ++i) {
    const Detection* d = *i;
    if (d->isClipped) continue;
    double xmod, ymod;
    pm->toWorldDerivs(d->xpix, d->ypix, xmod, ymod, derivs);
    xmod = d->xw - xmod;
    ymod = d->yw - ymod;

    chisq += xmod*xmod + ymod*ymod;
    for (int i=0; i<np; i++) {
      beta[i] += xmod*derivs(0,i) + ymod*derivs(1,i);
      for (int j=0; j<=i; j++) 
	alpha(i,j)+=derivs(0,i)*derivs(0,j) + derivs(1,i)*derivs(1,j);
    }
  }
  if (sigma != 1.) {
    double weight = 1./(sigma*sigma);
    chisq *= weight;
    beta *= weight;
    alpha *= weight;
  }
}
