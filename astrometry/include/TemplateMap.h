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
#include "Lookup1d.h"

namespace astrometry {
  // This class is a PixelMap with a small tabulated deviation that is function of,
  // and in direction of, either X axis, Y axis, or radius from some origin.
  // A single parameter, the scaling applied to the deviation, is initialized to 1.
  // To handle some DECam stuff, it can be two tables that apply above/below some x 
  // value but have a common scaling factor.

  // The LookupTables themselves are kept in a cache whenever one is found in an input file.
  // They are reused (but with independent scaling factor) if we make a new version.
  // They all stick around until the program ends.

  class TemplateMap: public PixelMap {
  public:

    // Build a template with single lookup table:
    TemplateMap(string tableName,
		string filename_,
		string name_);

    // Build a template that has two lookup tables for x above/below the split
    TemplateMap(double xSplit_,
		string lowName, string highName,
		string filename_,
		string name_);

    ~TemplateMap() {} // Don't delete tables, they belong to cache

    // Implement the PixelMap interface:
    virtual TemplateMap* duplicate() const;

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
    virtual DVector getParams() const {return DVector(1,scaling);}
    virtual int nParams() const {return 1;}

    static string type() {return "Template";}
    virtual string getType() const {return type();}

#ifdef USE_YAML
    static PixelMap* create(const YAML::Node& node,
			    bool& defaulted,
			    string name_="");
    virtual void write(YAML::Emitter& os) const;
#endif
    
  private:
    // If and where to split into two lookup tables at a specific X
    bool hasSplit;
    double xSplit;
    // Tables for below and above xSplit
    const lookup::Lookup1d* tableLow;
    const lookup::Lookup1d* tableHigh;

    // Origin of the tabulated data:
    string tableLowName;
    string tableHighName;
    string filename;

    // The parameter of this map is multiplicative factor on the table:
    double scaling;
  };

} // end namespace astrometry
#endif //  TEMPLATE_MAP_H
