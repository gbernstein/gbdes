// PhotoMap is a class representing a remapping of apparent magnitude.
// Each PhotoMap operates on a magnitude and a PhotoArguments instance, where the latter
// contains all the information we think the map might depend upon.
// PhotoMaps can have parameters.
//   For purposes of serialization, each derived class should define static methods
//   string Derived::mapType() giving a unique name to this class.
//   PhotoMap* create(istream& is, string name)  which produces member of the class from stream input
//   ... and implement the virtual method void write(ostream& os).
//   ... and also every PixelMap has a string name() bestowed on creation.  If name is not given, 
//   a name will be created from running counter.
#ifndef PHOTOMAP_H
#define PHOTOMAP_H

#include <string>
#include <list>
#include <iostream>
#include "Std.h"
#include "UseTMV.h"

//??? need setPixelStep back for maps???

namespace photometry {

  class PhotometryError: public std::runtime_error {
  public:
    PhotometryError(const string &m=""): 
      std::runtime_error("Photometry Error: " +m) {}
  };

  class PhotoArguments {
  public:
    double xDevice;
    double yDevice;
    double xExposure;
    double yExposure;
    double color;
    double airmass; // ??
  }

  class PhotoMap {
  public:
    PhotoMap(string name_="");
    virtual ~PhotoMap() {}
    // Return pointer to deep copy of self
    virtual PhotoMap* duplicate() const =0;

    // forward and inverse magnitude transformations:
    virtual double forward(double magIn, const PhotoArguments& args) const=0;
    // d(magOut)/d(magIn); base implementation is a finite-difference estimate:
    virtual double derivative(double magIn, const PhotoArguments& args) const;

    // Get parameters of the mapping.  Note default behavior of derived classes
    // is to ignore anything to do with map parameters.
    virtual void setParams(const DVector& p) {}
    // Note here that a zero-length vector is returned as parameters by default
    virtual DVector getParams() const {return DVector(0);}
    virtual int nParams() const {return 0;}

    // These calls include derivs of world w.r.t. parameters of map.  
    // Default implementation does not calculate derivs as default is no parameters
    virtual double forwardDerivs(double magIn, const PhotoArguments& args,
				 DVector& derivs) const {return forward(magIn, args);}

    // For serialization:
    virtual void write(std::ostream& os) const =0;
    string getName() const {return name;}
    virtual string getType() const =0;
  private:
    string name;
    static int anonymousCounter;
  };

  inline std::ostream& operator<<(std::ostream& os, const PhotoMap& pm) {pm.write(os); return os;}

  class IdentityMap: public PhotoMap {
  public:
    IdentityMap(): PhotoMap("Identity") {}
    virtual PhotoMap* duplicate() const {return new IdentityMap();}
    static string mapType() {return "Identity";}
    virtual string getType() const {return mapType();}
    static PhotoMap* create(std::istream& is, string name_="") {return new IdentityMap;}
    virtual double forward(double magIn, const PhotoArguments& args) const {return magIn;}
    virtual double derivative(double magIn, const PhotoArguments& args) const {return 1.;}

    virtual void write(std::ostream& os) const {} // Nothing to write
  };

  // Class adds a polynomial function of coordinates to the magnitude. 
  // Optionally this polynomial is multiplying the object color.

  class PolyMap: public PhotoMap {
  public:
    // Enum decides whether class is using device coordinates or exposure coords as polynomial arguments:
    enum ArgumentType {Device, Exposure};
    PolyMap(const poly2d::Poly2d& p,
	    ArgumentType argType,
	    string name="",
	    bool isColorTerm=false,
	    double tol_=0.0001);
    // Constructor with 1 order has terms with sum of x and y powers up to this order.
    // Constructor with 2 orders has all terms w/ powers of x up to orderx, y to ordery
    PolyMap(int orderx, int ordery,
	    ArgumentType argType,
	    string name="",
	    bool isColorTerm=false,
	    double tol_=0.0001);
    PolyMap(int order, 
	    ArgumentType argType,
	    string name="",
	    bool isColorTerm=false,
	    double tol_=0.0001);
    // Note that default tolerance is set to be 1 mas if world units are degrees.
    virtual PhotoMap* duplicate() const {return new PolyMap(*this);}
      
    ~PolyMap() {}

    virtual double forward(double magIn, const PhotoArguments& args) const;
    virtual double derivative(double magIn, const PhotoArguments& args) const;
    virtual double forwardDerivs(double magIn, const PhotoArguments& args,
				 DVector& derivs) const;


    void setParams(const DVector& p);
    DVector getParams() const;
    int nParams() const {return poly.nCoeffs();}

    // Access routines for this derived class:
    poly2d::Poly2d getPoly() const {return poly;}

    // Set tolerance in world coords for soln of inverse
    void setWorldTolerance(double wt) {worldTolerance=wt;}

    static string mapType() {return "Poly";}
    virtual string getType() const {return mapType();}
    static PhotoMap* create(std::istream& is, string name="");
    void write(std::ostream& os) const;

  private:
    poly2d::Poly2d poly;
    double worldTolerance;
    bool useExposureCoords;
    bool useColor;
  };

} // namespace photometry
#endif //PHOTOMAP_H
