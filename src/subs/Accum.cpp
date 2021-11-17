// Statistics accumulator class
#include "Accum.h"
#include "FitSubroutines.h"
using astrometry::WCS_UNIT;
using astrometry::RESIDUAL_UNIT;
using photometry::MMAG;

// Methods of the statistics-accumulator class

template <class S>
Accum<S>::Accum(): sumxw(0.), sumyw(0.), sumw(0.),
		   sumx(0.), sumy(0.),
		   sumxx(0.), sumyy(0.),
		   chisq(0.), sumdof(0.), n(0), ntot(0), nclipped(0) {}

// Specialization for Astro
template <>
void 
Accum<Astro>::add(const typename Astro::Detection& d) {
  ++ntot;
  if (d.isClipped) {
    ++nclipped;
    return;
  }

  // Get a rough sigma to use in weighting the centroids
  double sigma = d.getSigma();
  
  if (sigma <=0. || d.itsMatch->getDOF() <= 0) {
    // Residual statistics are meaningless if there
    // is no valid error nor multi-exposure fit
    return;
  }
  double wt = pow(sigma, -2.);
  
  auto dxy = d.residWorld(); // Returned in RESIDUAL_UNIT
  double dx = dxy[0];
  double dy = dxy[1];
  sumx += dx;
  sumy += dy;
  sumxw += dx*wt;
  sumyw += dy*wt;
  sumxx += dx*dx;
  sumyy += dy*dy;

  chisq += d.trueChisq();
  sumdof += d.expectedTrueChisq;
  ++n;
}

// Specialization for Photo
template<>
void 
Accum<Photo>::add(const Photo::Detection & d) {
  ++ntot;
  if (d.isClipped) {
    ++nclipped;
    return;
  }

  // Get a rough sigma to use in weighting the centroids
  double sigma = d.getSigma();
  
  if (sigma <=0. || d.itsMatch->getDOF() <= 0) {
    // Residual statistics are meaningless if there
    // is no valid error nor multi-exposure fit
    return;
  }
  double wt = pow(sigma, -2.);
  
  double dm = d.residMag();
  sumx += dm;
  sumxx += dm * dm;
  sumw += wt;

  chisq += d.trueChisq();
  sumdof += d.expectedTrueChisq;
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
    return sqrt( sumxx / n);
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
    oss << fixed << setprecision(3)
	<< " " << setw(6) << rms();
  } else {
    oss << fixed << setprecision(3)
	<< " " << setw(6) << rms()/MMAG;
  }
  oss << fixed << setprecision(5) 
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
