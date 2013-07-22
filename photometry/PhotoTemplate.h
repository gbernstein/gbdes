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

namespace photometry {

  // The class that can look up values of a 1d function; only for internal use of these classes
  class LookupTable {
  public:
    LookupTable() {}
    LookupTable(istream& is) {read(is);}
    void read(istream& is);
    double lookup(double arg) const;
  private:
    double nodeStart;
    double nodeStep;
    vector<double> values;
  };

  class PhotoTemplate1d: public PhotoMap {
  public:
    PhotoTemplate1d(string sourceFile, string name="");
    virtual ~PhotoTemplate1d() {}
    virtual PhotoMap* duplicate() const {return new PhotoTemplate1d(*this);}

    static string mapType() {return "Template1d";}
    virtual string getType() const {return mapType();}

    virtual double forward(double magIn, const PhotoArguments& args) const {
      return magIn + scaling*lookup( applyToX? args.xDevice : args.yDevice);
    }
    virtual double inverse(double magOut, const PhotoArguments& args) const {
      return magOut - scaling*lookup( applyToX? args.xDevice : args.yDevice);
    }

    virtual void setParams(const DVector& p) {
      Assert(p.size()==1);
      scaling = p[0];
    }
    virtual int nParams() const {return 1.;}
    virtual DVector getParams() const {return DVector(1,scaling);}
    virtual double forwardDerivs(double magIn, const PhotoArguments& args,
				 DVector& derivs) const {
      derivs = DVector(1,lookup( applyToX? args.xDevice : args.yDevice));
      return magIn + scaling * derivs[0];
    }

    static PhotoMap* create(std::istream& is, string name_="");
    virtual void write(std::ostream& os, int precision=DEFAULT_PRECISION) const;

  private:
    bool applyToX;	// Select function of X or of Y pixel coordinate
    string filename;	// File where template was read from
    double scaling;	// overall multiplicative factor applied to template
    LookupTable table;
    double lookup(double arg) const {return table.lookup(arg);} // table interpolator
    void read();	// read table from filename
  };

  class PhotoRings: public PhotoMap {
  public:
    PhotoRings();
    PhotoRings(string sourceFileLeft, string sourceFileRight, string name="");
    virtual ~PhotoRings() {}
    virtual PhotoMap* duplicate() const {return new PhotoRings(*this);}

    static string mapType() {return "Rings";}
    virtual string getType() const {return mapType();}

    virtual double forward(double magIn, const PhotoArguments& args) const {
      return magIn + scaling*lookup(args.xDevice,args.yDevice);
    }
    virtual double inverse(double magOut, const PhotoArguments& args) const {
      return magOut - scaling*lookup(args.xDevice,args.yDevice);
    }

    virtual void setParams(const DVector& p) {
      Assert(p.size()==1);
      scaling = p[0];
    }
    virtual int nParams() const {return 1.;}
    virtual DVector getParams() const {return DVector(1,scaling);}
    virtual double forwardDerivs(double magIn, const PhotoArguments& args,
				 DVector& derivs) const {
      derivs = DVector(1,lookup(args.xDevice, args.yDevice));
      return magIn + scaling * derivs[0];
    }

    static PhotoMap* create(std::istream& is, string name_="");
    virtual void write(std::ostream& os, int precision=DEFAULT_PRECISION) const;

  private:
    const double xSplit;	// *** I am going to hard-wire a split into two tables at x=1024.
    string fileLeft;	// File from which template was read
    string fileRight;	
    LookupTable tableLeft;
    LookupTable tableRight;
    double scaling;	// overall multiplicative factor applied to template
    double xCenterLeft;
    double yCenterLeft;
    double xCenterRight;
    double yCenterRight;
    double lookup(double xpix, double ypix) const; // table interpolator
    void read();	// read tables from 2 files
  };


} // namespace photometry

#endif
