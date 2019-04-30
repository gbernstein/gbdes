// Statistics accumulator class
#include "Accum.h"
#include "FitSubroutines.h"
using astrometry::WCS_UNIT;
using astrometry::RESIDUAL_UNIT;

// Methods of the statistics-accumulator class

template <class S>
Accum<S>::Accum(): sumxw(0.), sumyw(0.), 
		   sumx(0.), sumy(0.), sumw(0.),
		   sumxx(0.), sumyy(0.),
		   sum_m(0.), sum_mw(0.), sum_mm(0.), sum_mmw(0.),
		   chisq(0.), sumdof(0.), n(0), ntot(0), nclipped(0) {}

// Specialization for Astro
template <>
void 
Accum<Astro>::add(const typename Astro::Detection* d,
		  double magMean,
		  double wtot,
		  double dof) {
  ++ntot;
  if (d->isClipped) {
    ++nclipped;
    return;
  }

  // Get a rough sigma to use in weighting the centroids
  double sigma = d->getSigma();
  
  if (sigma <=0. || d->itsMatch->getDOF() <= 0) {
    // Residual statistics are meaningless if there
    // is no valid error nor multi-exposure fit
    return;
  }
  double wt = pow(sigma, -2.);
  
  auto dxy = d->residWorld(); // Returned in RESIDUAL_UNIT
  double dx = dxy[0];
  double dy = dxy[1];
  sumx += dx;
  sumy += dy;
  sumxw += dx*wt;
  sumyw += dy*wt;
  sumxx += dx*dx;
  sumyy += dy*dy;
  sumw  += wt;

  chisq += d->trueChisq();
  sumdof += d->expectedTrueChisq;
  ++n;
}

// Specialization for Photo
template<>
void 
Accum<Photo>::add(const Photo::Detection* d,
		  double magMean,
		  double wtot,
		  double dof) {
  ++ntot;
  if (d->isClipped) {
    ++nclipped;
    return;
  }
  if (wtot <=0. || dof<=0. ) return;  // Can't do anything more with this
  double dm = d->magOut - magMean;
  sum_m += dm;
  sum_mw += dm * d->wt;
  sum_mm += dm * dm;
  sum_mmw += dm * dm *d->wt;
  sumw += d->wt;
  double wtFrac = 1. - d->wt/wtot;
  if (wtFrac > 0.001) {
    chisq += dm * dm * d->wt / wtFrac;
  }
  sumx += d->args.xExposure;
  sumy += d->args.yExposure;
  sumdof += dof;
  ++n;
}

template <class S>
double 
Accum<S>::rms() const {
  if (n==0)
    return 0.;
  else if (S::isAstro)
    return sqrt( (sumxx+sumyy)/(2.*n));
  else
    return sqrt( sum_mm / n);
}
template <class S>
double
Accum<S>::reducedChisq() const {
  if (sumdof==0)
    return 0.;
  else if (S::isAstro)
    return chisq / sumdof;
  else
    return chisq / n;
}
template <class S>
string 
Accum<S>::summary() const {
  ostringstream oss;
  oss << setw(5) << n 
      << fixed << setprecision(0) << noshowpoint
      << " " << setw(6) << sumdof
      << fixed << setprecision(1) 
      << "  " << setw(5) << (n*100.)/ntot
      << fixed << setprecision(1) 
      << "  " << setw(4) << (nclipped*100.)/ntot;
  if (S::isAstro) {
    oss << fixed << setprecision(1)
	<< " " << setw(6) << rms();
  } else {
    oss << fixed << setprecision(1)
	<< " " << setw(6) << rms()*1000.;
  }
  oss << fixed << setprecision(2) 
      << "  " << setw(6) << reducedChisq();
  return oss.str();
}
template <class S>
string
Accum<S>::header() {
  return "  N     DOF   %use   %clip  RMS   ChiRed";
}
// Instantiate both cases
template
class Accum<Astro>;
template
class Accum<Photo>;
