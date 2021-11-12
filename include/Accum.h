// Code for statistics accumulator.
// Templated to either Astro or Photo, with
// appropriate specializations in its compilation.
// The Photo case will use the "x" variable as magnitude.

#ifndef ACCUM_H
#define ACCUM_H
#include "Std.h"

template <class S>
class Accum {
public:
  double sumx;
  double sumy;
  double sumxx;
  double sumyy;
  double sumxw;
  double sumyw;
  double sumw;
  int n; // Number of points with meaningful residuals
  int ntot; // total number of points
  int nclipped; // Points that were clipped
  double chisq;
  double sumdof;
  Accum();
  // Add to statistics this detection, with Match mean having
  // value of (xoff, yoff) with inverse variance wtot.
  // dof is fraction of a degree of freedom to assign for each xoff, yoff.
  void add(const typename S::Detection & d);
  double rms() const;
  double reducedChisq() const;
  string summary() const;
  static string header();
};


#endif
