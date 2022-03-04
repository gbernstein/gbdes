// Adapt Numerical Recipes Runge-Kutta solver to my C++ world.

#ifndef ODEINT_H
#define ODEINT_H

#include "LinearAlgebra.h"

namespace ode {

  class ODEFunction {
  public:
    virtual linalg::DVector operator()(double x, const linalg::DVector& y) const =0;
    virtual int order() const =0;
    virtual ~ODEFunction() {};
  };

  extern
  void rkck(const linalg::DVector& y,
	    const linalg::DVector& dydx,
	    double x,
	    double h,
	    linalg::DVector& yout,
	    linalg::DVector& yerr,
	    const ODEFunction& derivs);

  extern
  void rkqs(linalg::DVector& y, 
	    linalg::DVector& dydx,
	    double& x,
	    double htry,
	    double eps,
	    const linalg::DVector yscal,
	    double& hdid,
	    double& hnext,
	    const ODEFunction& derivs);

  // Numerical Recipes driver for adaptive-step Runge-Kutta
  extern
  void odeint(linalg::DVector& ystart,	// Initial values of y's
	      linalg::DVector& xp,	// Vector to be filled with saved x values
	      linalg::DMatrix& yp,	// Matrix to fill with y's at xp's
	      double dxsave,	// Minimum dx to save in xp & yp
	      int& kount,	// Number of steps saved in xp & yp
	      double x1,	// starting x
	      double x2,	// ending x
	      double eps,	// desired accuracy
	      double h1,	// suggested step in x
	      double hmin,	// minimum step in x (can be zero)
	      int& nok,		// number of good steps taken
	      int& nbad,	// number of bad steps taken
	      const ODEFunction& derivs);	// The function defining ODE

} // namespace
#endif
