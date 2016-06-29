// Statistics accumulator class
#include "Accum.h"
#include "FitSubroutines.h"

// Methods of the statistics-accumulator class

template <class S>
Accum<S>::Accum(): sumxw(0.), sumyw(0.), 
		   sumx(0.), sumy(0.), sumwx(0.), sumwy(0.), 
		   sumxx(0.), sumyy(0.), sumxxw(0.), sumyyw(0.),
		   sum_m(0.), sum_mw(0.), sum_mm(0.), sum_mmw(0.),
		   chisq(0.), sumdof(0.), n(0) {}

template <class S>
void 
Accum<S>::add(const typename S::Detection* d, double xoff, double yoff, double dof) {
  double dx = (d->xw-xoff);
  double dy = (d->yw-yoff);
  sumx += dx;
  sumy += dy;
  sumxw += dx*d->wtx;
  sumyw += dy*d->wty;
  sumxxw += dx*dx*d->wtx;
  sumyyw += dy*dy*d->wty;
  sumxx += dx*dx;
  sumyy += dy*dy;
  sumwx += d->wtx;
  sumwy += d->wty;
  chisq += (dx * dx * d->wtx + dy*dy*d->wty) / dof;
  sumdof += dof;
  ++n;
}

// Specialization for Photo
template<>
void 
Accum<Photo>::add(const Photo::Detection* d, double xoff, double yoff, double dof) {
  // mag comes in as xoff here
  double dm = d->magOut - xoff;
  sum_m += dm;
  sum_mw += dm * d->wt;
  sum_mm += dm * dm;
  sum_mmw += dm * dm *d->wt;
  sumw += d->wt;
  chisq += dm * dm * d->wt / dof;
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
    return chisq / (2.* n);
  else
    return chisq / n;
}
template <class S>
string 
Accum<S>::summary() const {
  ostringstream oss;
  oss << setw(5) << n 
      << fixed << setprecision(0) << noshowpoint
      << " " << setw(6) << sumdof;
  if (S::isAstro) {
    oss << fixed << setprecision(1)
	<< " " << setw(6) << rms()*1000.*DEGREE/ARCSEC;
  } else {
    oss << fixed << setprecision(1)
	<< " " << setw(6) << rms()*1000.*DEGREE;
  }
  oss << fixed << setprecision(2) 
      << " " << setw(6) << reducedChisq();
  return oss.str();
}
template <class S>
string
Accum<S>::header() {
  return "  N     DOF    RMS  ChiRed";
}
// Instantiate both cases
template
class Accum<Astro>;
template
class Accum<Photo>;
