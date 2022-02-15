// Code for a PixelMap that represents a small astrometric shift defined by a
// piecewise linearly interpolation between a series of variable values on
// linearly spaced nodes.
// Values of 0 are assumed at the boundary nodes.

// The assumption is that the deviations being entered here are small, and in particular always
// preserve one-to-one nature of the mapping.  This means that the difference between successive
// deviations must be < nodeStep/scaling.  Also requiring positive nodeStep for simplicity.

#ifndef PIECEWISE_MAP_H
#define PIECEWISE_MAP_H

#include "PixelMap.h"

namespace astrometry {
  // This class is a PixelMap with a piecewise linear shift as
  // a function of either X axis, Y axis, or radius from some origin.  Value at
  // each of the regularly spaced nodes is a free parameter.
  // A fixed value of 0 is placed in the bounding nodes
  // to turn off the correction outside the requested range.

  class PiecewiseMap: public PixelMap {
  public:
    enum Coordinate {X,Y,R}; // Which axis is the Table being applied to

    // Build a map with zeros at all nodes
    PiecewiseMap(string name_,
	     double argStart_,
	     double argEnd, 
	     double argStep_,
	     PiecewiseMap::Coordinate axis_,
	     double xCenter_ = 0.,
	     double yCenter_ = 0.);

    // Build a map with starting values
    // Value at first and last nodes is forced to zero.
    PiecewiseMap(string name_,
	     double argStart_,
	     double argStep_,
	     const DVector& values,
	     PiecewiseMap::Coordinate axis_,
	     double xCenter_ = 0.,
	     double yCenter_ = 0.);

    ~PiecewiseMap() {}

    // Implement the PixelMap interface:
    virtual PiecewiseMap* duplicate() const;

    void toWorld(double xpix, double ypix,
		 double& xworld, double& yworld,
		 double color=astrometry::NODATA) const;
    void toPix( double xworld, double yworld,
		double &xpix, double &ypix,
		double color=astrometry::NODATA) const;
    Matrix22 dWorlddPix(double xpix, double ypix,
			double color=astrometry::NODATA) const;    
    // Use base class matrix inversion for dPixdWorld

    virtual void toPixDerivs( double xworld, double yworld,
			      double &xpix, double &ypix,
			      DMatrix& derivs,
			      double color=astrometry::NODATA) const;
    virtual void toWorldDerivs(double xpix, double ypix,
			       double& xworld, double& yworld,
			       DMatrix& derivs,
			       double color=astrometry::NODATA) const;

    virtual void setParams(const DVector& p);
    virtual DVector getParams() const;
    virtual int nParams() const {return v.size()-2;}

    static string type() {return "Piecewise";}
    virtual string getType() const {return type();}

#ifdef USE_YAML
    static PixelMap* create(const YAML::Node& node,
			    bool& defaulted,
			    string name_="");
    virtual void write(YAML::Emitter& os) const;
#endif
    
  private:
    string name;

    Coordinate axis;
    // Lookups into the 1d table:
    double lookup(double arg) const;
    void valueSlope(double arg, double& value, double& slope) const;
    // Get the indexing info for an argument:
    void getIndex(double arg, int& index, double& frac) const;
    // Find location where arg + f(arg) = val:
    double inverse(double val) const;
    double argStart;
    double argStep;
    // For radial functions, define center:
    double xCenter;
    double yCenter;
    DVector v;
    DVector dvda;

    void setup(); // Recalculate slopes table after changing values
  };

} // end namespace astrometry
#endif //  PIECEWISE_MAP_H
