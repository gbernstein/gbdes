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
#include "LinearAlgebra.h"
#include "Poly2d.h"
#include "Bounds.h"
#include "yaml-cpp/yaml.h"

namespace photometry {

  typedef linalg::Matrix<double> DMatrix;
  typedef linalg::Vector<double> DVector;
  const double NODATA=-888.;
  
  class PhotometryError: public std::runtime_error {
  public:
    PhotometryError(const string &m=""): 
      std::runtime_error("Photometry Error: " +m) {}
  };

  class PhotoArguments {
  public:
    PhotoArguments(): xDevice(NODATA), yDevice(NODATA), xExposure(NODATA),
      yExposure(NODATA), color(NODATA) {}
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
    virtual void write(YAML::Emitter& os) const =0;
    string getName() const {return name;}
    virtual string getType() const {return "None";}
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

  inline YAML::Emitter& operator<<(YAML::Emitter& os, const PhotoMap& pm) {pm.write(os); return os;}

  class IdentityMap: public PhotoMap {
  public:
    IdentityMap(string name_="Identity"): PhotoMap(name_) {}
    virtual PhotoMap* duplicate() const {return new IdentityMap();}
    static string type() {return "Identity";}
    virtual string getType() const {return type();}
    virtual double forward(double magIn, const PhotoArguments& args) const {return magIn;}

    static PhotoMap* create(const YAML::Node& node, string name_="") {
      return name_.empty() ? new IdentityMap() : new IdentityMap(name_);}
    virtual void write(YAML::Emitter& os) const {
      os << YAML::BeginMap << YAML::Key << "Type" << YAML::Value << type()
	 << YAML::EndMap;
    }
  };

  // Map that shifts magnitudes by a scalar.  Scalar is a free parameter.
  class ConstantMap: public PhotoMap {
  public:
    ConstantMap(double c_=0., string name=""): PhotoMap(name), c(c_) {}
    virtual PhotoMap* duplicate() const {return new ConstantMap(*this);}
    static string type() {return "Constant";}
    virtual string getType() const {return type();}
    virtual double forward(double magIn, const PhotoArguments& args) const {return magIn+c;}
    virtual double inverse(double magOut, const PhotoArguments& args) const {return magOut-c;}
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
    static PhotoMap* create(const YAML::Node& node, string name_="");
    virtual void write(YAML::Emitter& os) const;
  private:
    double c;
  };


  // A ColorTerm will take any other kind of PhotoMap and multiply its shift
  // by the color given in the PhotoArguments.
  // The ColorTerm will assume ownership of the PhotoMap that it wraps and will delete it.  
  class ColorTerm: public PhotoMap {
  public:
    ColorTerm(PhotoMap* pm_, double referenceColor=0., string name=""):
      PhotoMap(name),
      pm(pm_),
      reference(referenceColor) {
      /* ??? possibly check not wrapping another ColorMap */
    }
    virtual ~ColorTerm() {delete pm;}
    // Access the wrapped PhotoMap:
    const PhotoMap* map() const {return pm;}
    const double refColor() const {return reference;}
    
    virtual PhotoMap* duplicate() const {return new ColorTerm(pm->duplicate(),
							      reference,
							      getName());}
    static string type() {return "Color";}
    virtual string getType() const {return type();}
    virtual bool needsColor() const {return true;}

    double forward(double magIn, const PhotoArguments& args) const {
      return magIn + (args.color-reference) * (pm->forward(magIn,args)-magIn);
    }

    double derivative(double magIn, const PhotoArguments& args) const {
	return 1. + (args.color-reference) * (pm->derivative(magIn,args)-1.);
    }

    double forwardDerivs(double magIn, const PhotoArguments& args, 
			 DVector& derivs) const {
      double magOut = pm->forwardDerivs(magIn, args, derivs);
      derivs *= (args.color-reference);
      return magIn + (args.color-reference)*(magOut-magIn);
    }

    int nParams() const {return pm->nParams();}
    void setParams(const DVector& p) {pm->setParams(p);}
    DVector getParams() const {return pm->getParams();}
    
    // Will not have Create() call since color maps will be built by PhotoMapCollection.
    // But for some convenience, provide a separate serialization.
    virtual void write(YAML::Emitter& os) const;

  private:
    PhotoMap* pm;
    double reference;	// Color value yielding zero correction
    /*Hide - do not want confused ownership*/
    ColorTerm(const ColorTerm& rhs) = delete;
    void operator=(const ColorTerm& rhs) = delete;
  };

  // Class adds a polynomial function of coordinates to the magnitude. 
  // Optionally this polynomial is multiplying the object color.

  class PolyMap: public PhotoMap {
  public:
    PolyMap(const poly2d::Poly2d& p,
	    ArgumentType argType,
	    string name="",
	    Bounds<double> domain=Bounds<double>(-1.,1.,-1.,1.));
    // Constructor with 1 order has terms with sum of x and y powers up to this order.
    // Constructor with 2 orders has all terms w/ powers of x up to orderx, y to ordery
    PolyMap(int orderx, int ordery,
	    ArgumentType argType,
	    string name="",
	    Bounds<double> domain=Bounds<double>(-1.,1.,-1.,1.));
    PolyMap(int order, 
	    ArgumentType argType,
	    string name="",
	    Bounds<double> domain=Bounds<double>(-1.,1.,-1.,1.));
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

    // Set up domain scaling
    void setDomain(Bounds<double> domain);
    // Return rectangular bounds that rescale to (-1,1)
    Bounds<double> getDomain() const;

    static string type() {return "Poly";}
    virtual string getType() const {return type();}
    static PhotoMap* create(const YAML::Node& node, string name="");
    void write(YAML::Emitter& os) const;
  private:
    poly2d::Poly2d poly;
    bool usePixelCoords; // True to use device/pixel coords, false for exposure
    void setArgumentType(ArgumentType argType);
    // scale input data so that x->xscale*(x-shift), nominally (-1,1.)
    double xshift, xscale;
    double yshift, yscale;
    void rescale(double& xpix, double& ypix) const {
      // Apply domain scaling
      xpix = xscale * (xpix-xshift);
      ypix = yscale * (ypix-yshift);
    }
  };

} // namespace photometry
#endif //PHOTOMAP_H
