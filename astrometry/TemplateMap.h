// Code for a PixelMap that represents a small astrometric shift defined by a
// template.  The template(s) are interpolated between values at a given grid of nodes.
// An overall scaling of the astrometric shift is a free parameter of the map.
// Currently will implement linear interpolation between equally-spaced nodes.
// The edge values are applied to any point outside the bounds of the grid.

// The assumption is that the deviations being entered here are small, and in particular always
// preserve one-to-one nature of the mapping.  This means that the difference between successive
// deviations must be < nodeStep/scaling.  Also requiring positive nodeStep for simplicity.

#ifndef TEMPLATE_MAP_H
#define TEMPLATE_MAP_H

#include "PixelMap.h"

namespace astrometry {
  class TemplateMap1d: public PixelMap {
  public:
    enum Coordinate {X,Y}; // Which axis is the Table being applied to

    TemplateMap1d(Coordinate c, 
		  double nodeStart_, double nodeStep_, 
		  const vector<double>& deviations_,
		  string name="");
    ~TemplateMap1d() {}

    // Implement the PixelMap interface:
    virtual TemplateMap1d* duplicate() const;
    static string mapType() {return "Template1d";}
    virtual string getType() const {return mapType();}
    static PixelMap* create(std::istream& is, string name_="");
    virtual void write(std::ostream& os, int precision=PixelMap::DEFAULT_PRECISION) const;

    void toWorld(double xpix, double ypix,
		 double& xworld, double& yworld) const;
    void toPix( double xworld, double yworld,
		double &xpix, double &ypix) const;
    Matrix22 dPixdWorld(double xworld, double yworld) const;    
    Matrix22 dWorlddPix(double xpix, double ypix) const;    

    virtual void toPixDerivs( double xworld, double yworld,
			      double &xpix, double &ypix,
			      DMatrix& derivs) const;
    virtual void toWorldDerivs(double xpix, double ypix,
			       double& xworld, double& yworld,
			       DMatrix& derivs) const;

    virtual void setParams(const DVector& p);
    virtual DVector getParams() const {return DVector(1,scaling);}
    virtual int nParams() const {return 1;}

  private:
    bool applyToX;	//If false, applying to Y coord
    double nodeStart;
    double nodeStep;
    vector<double> deviations;
    double scaling;
    double lookup(double arg) const;  // Do a lookup into the table.
    double slope(double arg) const;  // Just looking up slope at the arg
    double reverse_lookup(double val) const;	// Invert the map.
  };

  // ??? Do Template2d on a grid.
  // ??? Radial grid
  // ??? **Some way of caching information for each point???

} // end namespace astrometry
#endif //  TEMPLATE_MAP_H
