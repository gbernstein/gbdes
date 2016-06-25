// Subroutines used by PhotoFit.cpp

#include "PhotoSubs.h"

using namespace photometry;

// Methods of the statistics-accumulator class

Accum::Accum(): sum_m(0.), sum_mw(0.),
		sum_mm(0.), sum_mmw(0.),
		sum_w(0.),
		chisq(0.), 
		sum_x(0.), sum_y(0.),
		n(0) {}
void 
Accum::add(const Detection* d, double magoff, double dof) {
  double dm = d->magOut - magoff;
  sum_m += dm;
  sum_mw += dm * d->wt;
  sum_mm += dm * dm;
  sum_mmw += dm * dm *d->wt;
  sum_w += d->wt;
  chisq += dm * dm * d->wt / dof;
  sum_x += d->args.xExposure;
  sum_y += d->args.yExposure;
  ++n;
}
double 
Accum::rms() const {
  return n > 0 ? sqrt( sum_mm/n ) : 0.;
}
double 
Accum::reducedChisq() const {
  return n>0 ? chisq / n : 0.;
}
string 
Accum::summary() const {
  ostringstream oss;
  double dm = 0., sigm=0.;
  double x=0., y=0.;
  if (n>0) {
    dm = sum_mw / sum_w;
    sigm = 1./sqrt(sum_w);
    x = sum_x / n;
    y = sum_y / n;
  }
  oss << setw(5) << n 
      << fixed << setprecision(1)
    /** << " " << setw(5) << dm*1000.
	<< " " << setw(5) << sigm*1000.**/
      << " " << setw(5) << rms()*1000.
      << setprecision(2) 
      << " " << setw(5) << reducedChisq()
      << setprecision(5) << showpoint << showpos
      << " " << setw(9) << x 
      << " " << setw(9) << y ;
  return oss.str();
}
