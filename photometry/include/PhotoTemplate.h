// Code for PhotoMaps that allow scaled versions of fixed functions of x, y, or radius from
// a chosen pixel position.
//   It will be assumed that the templates will be read from files.  For compactness of the
// serialization, we will just serialize the filename (and the scaling factor that we apply
// to the template).
//   Nodes of tables assumed to be equally spaced, so files will just give first argument,
// the spacing between args of nodes, and then a list of values at nodes.  Linear interpolation done
// between the nodes.  Values outside tables take the endpoint values.
//
// The format of the files holding templates is:
// First line has "X" or "Y" or "R" to say whether it's function of x, y, or radius.  If "R", 
//   then two additional floats give the center position for radii.
// Then comes a line with the nodeStart and nodeStep (i.e. first and delta of arguments)
// Successive lines are 1 node value per line.
// blank lines and # are skipped.

#ifndef PHOTOTEMPLATE_H
#define PHOTOTEMPLATE_H

#include "PhotoMap.h"
#include "Lookup1d.h"

namespace photometry {

  class PhotoTemplate: public PhotoMap {
  public:
    // Build a template with single lookup table:
    PhotoTemplate(string tableName,
		  string filename_,
		  string name_);

    // Build a template that has two lookup tables for x above/below the split
    PhotoTemplate(double xSplit_,
		  string lowName, string highName,
		  string filename_,
		  string name_);

    virtual ~PhotoTemplate() {} // Don't delete tables, they belong to cache
    virtual PhotoMap* duplicate() const {return new PhotoTemplate(*this);}

    static string type() {return "Template";}
    virtual string getType() const {return type();}

    virtual double forward(double magIn, const PhotoArguments& args) const {
      return magIn + scaling*value(args);
    }
    virtual double inverse(double magOut, const PhotoArguments& args) const {
      return magOut - scaling*value(args);
    }

    virtual void setParams(const DVector& p) {
      Assert(p.size()==1);
      scaling = p[0];
    }
    virtual int nParams() const {return 1.;}
    virtual DVector getParams() const {return DVector(1,scaling);}
    virtual double forwardDerivs(double magIn, const PhotoArguments& args,
				 DVector& derivs) const {
      derivs = DVector(1,value(args));
      return magIn + scaling * derivs[0];
    }

    virtual void write(YAML::Emitter& os) const;
    static PhotoMap* create(const YAML::Node& node, string name);

  private:
    // If and where to split into two lookup tables at a specific X
    bool hasSplit;
    double xSplit;
    // Tables for below and above xSplit
    const lookup::Lookup1d* tableLow;
    const lookup::Lookup1d* tableHigh;

    double value(const PhotoArguments& args) const {
      if (hasSplit && args.xDevice>xSplit)
	return tableHigh->value(args.xDevice, args.yDevice);
      else
	return tableLow->value(args.xDevice, args.yDevice);
    }
    // Origin of the tabulated data:
    string tableLowName;
    string tableHighName;
    string filename;

    // The parameter of this map is multiplicative factor on the table:
    double scaling;
  };

} // namespace photometry

#endif
