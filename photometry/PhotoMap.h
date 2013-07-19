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
#include "Poly2d.h"

namespace photometry {

  class PhotometryError: public std::runtime_error {
  public:
    PhotometryError(const string &m=""): 
      std::runtime_error("Photometry Error: " +m) {}
  };

  class PhotoArguments {
  public:
    PhotoArguments(): xDevice(NODATA), yDevice(NODATA), xExposure(NODATA), yExposure(NODATA),
		      color(NODATA) {}
    const static double NODATA;
    double xDevice;
    double yDevice;
    double xExposure;
    double yExposure;
    double color;
  };


  class PhotoMap {
  public:
    PhotoMap(string name_="");
    virtual ~PhotoMap() {}
    // Return pointer to deep copy of self
    virtual PhotoMap* duplicate() const =0;

    // Enum decides whether class is using device coordinates or exposure coords as position arguments:
    enum ArgumentType {Device, Exposure};

    // forward transformations:
    virtual double forward(double magIn, const PhotoArguments& args) const=0;
    // There will be a default inverse which assumes that transform is an additive quantity
    // independent of the input mag:
    virtual double inverse(double magOut, const PhotoArguments& args) const;
    // d(magOut)/d(magIn); base implementation is unity
    virtual double derivative(double magIn, const PhotoArguments& args) const {return 1.;}

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

    // Does the transform require object colors?  Default is no.
    virtual bool needsColor() const {return false;}

    // For serialization:
    static const int DEFAULT_PRECISION=8;
    virtual void write(std::ostream& os, int precision=DEFAULT_PRECISION) const =0;
    string getName() const {return name;}
    virtual string getType() const =0;
  protected:
    // This can be used by derived classes to do inversion if they are nonlinear in magnitude:
    // Notice 1 mmag default tolerance.
    void NewtonInverse(double magOut,
		       double& magIn,
		       const PhotoArguments& args,
		       double magTolerance=0.001) const;
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
    virtual void write(std::ostream& os, int precision) const {} // Nothing to write
  };

  // Map that shifts magnitudes by a scalar.  Scalar is a free parameter.
  class ConstantMap: public PhotoMap {
  public:
    ConstantMap(double c_=0., string name=""): PhotoMap(name), c(c_) {}
    virtual PhotoMap* duplicate() const {return new ConstantMap(*this);}
    static string mapType() {return "Constant";}
    virtual string getType() const {return mapType();}
    static PhotoMap* create(std::istream& is, string name_="");
    virtual void write(std::ostream& os, int precision) const;
    virtual double forward(double magIn, const PhotoArguments& args) const {return magIn+c;}
    virtual double inverse(double magIn, const PhotoArguments& args) const {return magIn+c;}
    virtual double forwardDerivs(double magIn, const PhotoArguments& args,
				 DVector& derivs) const {
      derivs = DVector(1,1.);
      return magIn + c;
    }

    virtual int nParams() const {return 1;}
    virtual DVector getParams() const {return DVector(1,c);}
    virtual void setParams(const DVector& p) {
      Assert(p.size()==1);
      c = p[0];
    }
  private:
    double c;
  };


  // A ColorTerm will take any other kind of PhotoMap and multiply its shift by the color given
  // in the PhotoArguments.
  // The ColorTerm will assume ownership of the PhotoMap that it wraps and will delete it.  
  class ColorTerm: public PhotoMap {
  public:
    ColorTerm(PhotoMap* pm_, string name=""): PhotoMap(name),
					      pm(pm_) {
      /* ??? possibly check not wrapping another ColorMap */
    }
    virtual ~ColorTerm() {delete pm;}
    // Access the wrapped PhotoMap:
    const PhotoMap* map() const {return pm;}

    virtual PhotoMap* duplicate() const {return new ColorTerm(pm->duplicate(), getName());}
    static string mapType() {return "Color";}
    virtual string getType() const {return mapType();}
    /* will not have Create() call since color maps will be built by PhotoMapCollection */
    virtual void write(std::ostream& os, int precision) const {
      throw PhotometryError("ColorTerm " + getName() + " should not be getting serialized");
    }
    virtual bool needsColor() const {return true;}

    double forward(double magIn, const PhotoArguments& args) const {
      return magIn + args.color * (pm->forward(magIn,args)-magIn);
    }

    double derivative(double magIn, const PhotoArguments& args) const {
	return 1. + args.color * (pm->derivative(magIn,args)-1.);
    }

    double forwardDerivs(double magIn, const PhotoArguments& args, 
			 DVector& derivs) const {
      double magOut = pm->forwardDerivs(magIn, args, derivs);
      derivs *= args.color;
      return magIn + args.color*(magOut-magIn);
    }

    int nParams() const {return pm->nParams();}
    void setParams(const DVector& p) {pm->setParams(p);}
    DVector getParams() const {return pm->getParams();}
    
  private:
    PhotoMap* pm;
    /*Hide - do not want confused ownership*/
    ColorTerm(const ColorTerm& rhs);
    void operator=(const ColorTerm& rhs);
  };

  // Class adds a polynomial function of coordinates to the magnitude. 
  // Optionally this polynomial is multiplying the object color.

  class PolyMap: public PhotoMap {
  public:
    PolyMap(const poly2d::Poly2d& p,
	    ArgumentType argType,
	    string name="");
    // Constructor with 1 order has terms with sum of x and y powers up to this order.
    // Constructor with 2 orders has all terms w/ powers of x up to orderx, y to ordery
    PolyMap(int orderx, int ordery,
	    ArgumentType argType,
	    string name="");
    PolyMap(int order, 
	    ArgumentType argType,
	    string name="");
    virtual PhotoMap* duplicate() const {return new PolyMap(*this);}
      
    virtual ~PolyMap() {}

    virtual double forward(double magIn, const PhotoArguments& args) const;
    virtual double forwardDerivs(double magIn, const PhotoArguments& args,
				 DVector& derivs) const;


    void setParams(const DVector& p);
    DVector getParams() const {return poly.getC();}
    int nParams() const {return poly.nCoeffs();}

    // Access routines for this derived class:
    poly2d::Poly2d getPoly() const {return poly;}

    static string mapType() {return "Poly";}
    virtual string getType() const {return mapType();}
    static PhotoMap* create(std::istream& is, string name="");
    void write(std::ostream& os, int precision=PhotoMap::DEFAULT_PRECISION) const;

  private:
    poly2d::Poly2d poly;
    bool useExposureCoords;
    void setArgumentType(ArgumentType argType);
  };

} // namespace photometry
#endif //PHOTOMAP_H
