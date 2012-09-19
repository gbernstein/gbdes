// $Id: odeint.h,v 1.1.1.1 2009/10/30 21:20:52 garyb Exp $
// Adapt Numerical Recipes Runge-Kutta solver to my C++ world.

#ifndef ODEINT_H
#define ODEINT_H

#include "UseTMV.h"

namespace ode {

  class ODEFunction {
  public:
    virtual DVector operator()(double x, const DVector& y) const =0;
    virtual int order() const =0;
    virtual ~ODEFunction() {};
  };

  extern
  void rkck(const DVector& y,
	    const DVector& dydx,
	    double x,
	    double h,
	    DVector& yout,
	    DVector& yerr,
	    const ODEFunction& derivs);

  extern
  void rkqs(DVector& y, 
	    DVector& dydx,
	    double& x,
	    double htry,
	    double eps,
	    const DVector yscal,
	    double& hdid,
	    double& hnext,
	    const ODEFunction& derivs);

  // Numerical Recipes driver for adaptive-step Runge-Kutta
  extern
  void odeint(DVector& ystart,	// Initial values of y's
	      DVector& xp,	// Vector to be filled with saved x values
	      DMatrix& yp,	// Matrix to fill with y's at xp's
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
