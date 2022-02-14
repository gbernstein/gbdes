// Code for a PixelMap that represents a small astrometric shift defined by a
// piecewise linearly interpolation between a series of variable values on
// linearly spaced nodes.
// Values of 0 are assumed at the boundary nodes.

// The assumption is that the deviations being entered here are small, and in particular always
// preserve one-to-one nature of the mapping.  This means that the difference between successive
// deviations must be < nodeStep/scaling.  Also requiring positive nodeStep for simplicity.

#ifndef PHOTOPIECEWISE_H
#define PHOTOPIECEWISE_H

#include "PhotoMap.h"

namespace photometry {
  // This class is a PhotoMap with a piecewise linear photometric correction as
  // a function of either X axis, Y axis, or radius from some origin.  Value at
  // each of the regularly spaced nodes is a free parameter.
  // A fixed value of 0 is placed in the bounding nodes
  // to turn off the correction outside the requested range.

  class PhotoPiecewise: public PhotoMap {
  public:
    enum Coordinate {X,Y,R}; // Which axis is the Table being applied to

    // Build a map with zeros at all nodes
    PhotoPiecewise(string name_,
	     double argStart_,
	     double argEnd,
	     double argStep_,
	     PhotoPiecewise::Coordinate axis_,
	     double xCenter_ = 0.,
	     double yCenter_ = 0.);

    // Build a map with starting values
    // Value at first and last nodes is forced to zero.
    PhotoPiecewise(string name_,
	     double argStart_,
	     double argStep_,
	     const DVector& values,
	     PhotoPiecewise::Coordinate axis_,
	     double xCenter_ = 0.,
	     double yCenter_ = 0.);

    virtual ~PhotoPiecewise() {}

    // Implement the PhotoMap interface:
    virtual PhotoPiecewise* duplicate() const;

    virtual double forward(double magIn, const PhotoArguments& args) const {
      return magIn + lookup(getarg(args));
    }
    virtual double inverse(double magOut, const PhotoArguments& args) const {
      return magOut - lookup(getarg(args));
    }
    virtual double forwardDerivs(double magIn, const PhotoArguments& args,
				 DVector& derivs) const;

    virtual void setParams(const DVector& p);
    virtual DVector getParams() const;
    virtual int nParams() const {return v.size()-2;}

    static string type() {return "Piecewise";}
    virtual string getType() const {return type();}

    static PhotoMap* create(const YAML::Node& node, string name_="");
    virtual void write(YAML::Emitter& os) const;

  private:
    string name;

    Coordinate axis;
    double argStart;
    double argStep;
    // For radial functions, define center:
    double xCenter;
    double yCenter;
    DVector v;
    // Lookups into the 1d table:
    double lookup(double arg) const;
    double getarg(const PhotoArguments& args) const; // Get the X/Y/R value we want
  };

} // end namespace photometry
#endif //  PHOTOPIECEWISE_H
