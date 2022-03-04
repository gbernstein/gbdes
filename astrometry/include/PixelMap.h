// PixelMap is a base class for maps from 2d "pixel" space to 2d "world" space.
// Maps should be bijective over useful range as inverse should be defined.
// Designed to be fit to data by adjusting parameters of the map.
// pixelStep is the increment in pixel space that should be used for taking numerical derivatives.
//   For purposes of serialization, each derived class should define static method
//   string Derived::type() giving a unique name to this class.
//   PixelMap* create(istream& is, string name)  which produces member of the class from stream input
//   ... and implement the virtual method write(YAML::Emitter& os) for serialization
//       and create(YAML::Node& node) for deserialization
//   ... and also every PixelMap has a string name() bestowed on creation.  If name is not given, 
//   a name will be created from running counter.
#ifndef PIXELMAP_H
#define PIXELMAP_H

#include <string>
#include <iostream>
#ifdef USE_YAML
#include "yaml-cpp/yaml.h"
#endif
#include "Std.h"
#include "LinearAlgebra.h"
#include "Astrometry.h"
// Note that Astrometry.h has typdefs for the 2- and 3-element small vectors & matrices

namespace astrometry {

  class PixelMap {
  public:
    PixelMap(string name_="");
    virtual ~PixelMap() {}
    // Return pointer to deep copy of self
    virtual PixelMap* duplicate() const =0;

    // pixel coords to world coords map.
    // Derived classes that use the color are responsible for checking
    // that they did not receive the default NODATA value, and declarations
    // must include the default value too.
    virtual void toWorld(double xpix, double ypix,
			 double& xworld, double& yworld,
			 double color=astrometry::NODATA) const=0;
    // inverse map:
    virtual void toPix( double xworld, double yworld,
			double &xpix, double &ypix,
			double color=astrometry::NODATA) const=0;
    // Default implementations take finite derivatives:
    virtual Matrix22 dPixdWorld(double xworld, double yworld, 
				double color=astrometry::NODATA) const;
    virtual Matrix22 dWorlddPix(double xpix, double ypix, 
				double color=astrometry::NODATA) const;
    // Solid angle of a pixel (sr)
    double pixelArea(double xpix, double ypix, 
		     double color=astrometry::NODATA) const;
    // Step size for derivatives in "pixel" space:
    virtual double getPixelStep() const {return pixStep;}
    virtual void setPixelStep(double ps) {pixStep=ps;}

    // Get parameters of the mapping.  Note default behavior of derived classes
    // is to ignore anything to do with map parameters.
    virtual void setParams(const DVector& p) {}
    // Note here that a zero-length vector is returned as parameters by default
    virtual DVector getParams() const {return DVector(0);}
    virtual int nParams() const {return 0;}
    // Does the transform require object colors?  Default is no.
    virtual bool needsColor() const {return false;}

    // These calls include derivs of world w.r.t. parameters of map.  
    // Default implementation does not calculate derivs as default is no parameters
    virtual void toWorldDerivs(double xpix, double ypix,
			       double& xworld, double& yworld,
			       DMatrix& derivs,
			       double color=astrometry::NODATA) const {
      toWorld(xpix,ypix,xworld,yworld,color);}
    // In toPixDerivs, the derivatives are of world coords w.r.t. map parameters, even
    // though the call is going from world to pix coords:
    virtual void toPixDerivs( double xworld, double yworld,
			      double &xpix, double &ypix,
			      DMatrix& derivs,
			      double color=astrometry::NODATA) const {
      toPix(xworld,yworld,xpix,ypix,color);}


    // For serialization:  the precision is going to be a suggestion to derived
    // classes about number of digits to use in outputs.
    static const int DEFAULT_PRECISION=8;
#ifdef USE_YAML
    virtual void write(YAML::Emitter& os) const =0;
#endif
    string getName() const {return name;}
    virtual string getType() const {return "None";}
    
  private:
    string name;
    double pixStep;
    static int anonymousCounter;
  protected:
    // Invert the pix->world map using Newton iteration.
    // Likely that many derived PixelMap classes will want this.
    void NewtonInverse(double xworld, double yworld, 
		       double& xpix, double& pix,
		       double worldTolerance,
		       double color) const;
    // A helper function to check whether color has been defaulted
    void checkColor(double color) const {
      if (color==astrometry::NODATA) 
	throw AstrometryError("Defaulted color for color-dependent PixelMap");
    }
  };

#ifdef USE_YAML
  inline YAML::Emitter& operator<<(YAML::Emitter& os, const PixelMap& pm) {
    pm.write(os);
    return os;
  }
#endif
  
  class IdentityMap: public PixelMap {
  public:
    IdentityMap(string name_="Identity"): PixelMap(name_) {}
    virtual PixelMap* duplicate() const {return new IdentityMap(*this);}
    static string type() {return "Identity";}
    void toWorld(double xpix, double ypix,
		 double& xworld, double& yworld,
		 double color=astrometry::NODATA) const;
    void toPix( double xworld, double yworld,
		double &xpix, double &ypix,
		double color=astrometry::NODATA) const;
    Matrix22 dPixdWorld(double xworld, double yworld, 
			double color=astrometry::NODATA) const;    
    Matrix22 dWorlddPix(double xpix, double ypix, 
			double color=astrometry::NODATA) const;    
#ifdef USE_YAML
    static PixelMap* create(const YAML::Node& node,
			    bool& defaulted,
			    string name_="") {
      defaulted = false;
      return name_.empty() ? new IdentityMap() : new IdentityMap(name_);}
    virtual void write(YAML::Emitter& os) const {
      os << YAML::BeginMap << YAML::Key << "Type" << YAML::Value << type()
	 << YAML::EndMap;
    }
#endif
  };

  class ConstantMap: public PixelMap {
    // Class representing a constant shift of position
  public:
    ConstantMap(string name_, double dx=0., double dy=0.): PixelMap(name_),
							   dxy(2) {
      dxy[0] = dx; dxy[1] = dy;
    }
    virtual PixelMap* duplicate() const {return new ConstantMap(*this);}
    static string type() {return "Constant";}
    virtual string getType() const {return type();}
    int nParams() const {return 2;}
    void setParams(const DVector& p) {Assert(p.size()==2); dxy = p;}
    DVector getParams() const {return dxy;}

    void toWorld(double xpix, double ypix,
		 double& xworld, double& yworld,
		 double color=astrometry::NODATA) const {
      xworld = xpix + dxy[0];
      yworld = ypix + dxy[1];
    }
    void toPix( double xworld, double yworld,
		double &xpix, double &ypix,
		double color=astrometry::NODATA) const {
      xpix = xworld - dxy[0];
      ypix = yworld - dxy[1];
    }
    Matrix22 dPixdWorld(double xworld, double yworld, 
			double color=astrometry::NODATA) const {
      return Matrix22().setToIdentity();
    }

    Matrix22 dWorlddPix(double xpix, double ypix, 
			double color=astrometry::NODATA) const {
      return Matrix22().setToIdentity();
    }
    void toPixDerivs( double xworld, double yworld,
		      double &xpix, double &ypix,
		      DMatrix& derivs,
		      double color=astrometry::NODATA) const {
      toPix(xworld,yworld,xpix,ypix,color);
      Assert(derivs.ncols()==2 && derivs.nrows()==2);
      derivs.setToIdentity();
    }
    void toWorldDerivs(double xpix, double ypix,
		       double& xworld, double& yworld,
		       DMatrix& derivs,
		       double color=astrometry::NODATA) const {
      toWorld(xpix,ypix,xworld,yworld,color);
      Assert(derivs.ncols()==2 && derivs.nrows()==2);
      derivs.setToIdentity();
    }
  
#ifdef USE_YAML
    static PixelMap* create(const YAML::Node& node,
			    bool& defaulted,
			    string name_="");
    virtual void write(YAML::Emitter& os) const {
      vector<double> vv(2);
      vv[0] = dxy[0]; vv[1] = dxy[1];
      os << YAML::BeginMap << YAML::Key << "Type" << YAML::Value << type()
	 << YAML::Key << "Parameters" << YAML::Flow << YAML::Value << vv
	 << YAML::EndMap;
    }
#endif    
      
  private:
    DVector dxy;

  };

  // PixelMap that takes coordinates in one projection
  // into another.  No adjustable parameters
  class ReprojectionMap: public PixelMap {
  public:
    // On creation, give reference to a set of coordinates in the
    // source (pix) and destination (world) systems.  
    // scaleFactor will be applied to the lat/lon in both systems
    // to put them into the native radian units of this class.
    // For instance enter AstronomicalConstants::DEGREE to
    // have coordinates in degrees.
    // setPixelScale is 1 arcsecond for any numerical derivatives of this projection,
    // although derivatives are done analytically.
    ReprojectionMap(const SphericalCoords& pixCoords,
		    const SphericalCoords& worldCoords,
		    double scale_=1.,
		    string name_=""):
      PixelMap(name_),
      pix(pixCoords.duplicate()), world(worldCoords.duplicate()), 
      scaleFactor(scale_)  {
      setPixelStep(ARCSEC/scale_);
    }
    ReprojectionMap(const ReprojectionMap& rhs);
    ReprojectionMap& operator=(const ReprojectionMap& rhs);
    ~ReprojectionMap() {delete pix; delete world;}
    virtual PixelMap* duplicate() const;

    void toWorld(double xpix, double ypix,
		 double& xworld, double& yworld,
		 double color=astrometry::NODATA) const;
    void toPix( double xworld, double yworld,
		double &xpix, double &ypix,
		double color=astrometry::NODATA) const;
    Matrix22 dPixdWorld(double xworld, double yworld, 
			double color=astrometry::NODATA) const;    
    Matrix22 dWorlddPix(double xpix, double ypix, 
			double color=astrometry::NODATA) const;    

    static string type() {return "Reprojection";}
#ifdef USE_YAML
    static PixelMap* create(const YAML::Node& node,
			    bool& defaulted,
			    string name="");
    virtual void write(YAML::Emitter& os) const;
#endif
    virtual string getType() const {return type();}
  private:
    SphericalCoords* pix;
    SphericalCoords* world;
    double scaleFactor;
  };

  // A ColorTerm will take any other kind of PixelMap and multiply its shift by the color given
  // in the arguments.
  // The ColorTerm will assume ownership of the PixelMap that it wraps and will delete it.  
  class ColorTerm: public PixelMap {
  public:
    ColorTerm(PixelMap* pm_, double referenceColor=0., string name=""):
      PixelMap(name),
      pm(pm_),
      reference(referenceColor) {}
    virtual ~ColorTerm() {delete pm;}
    // Access the wrapped PhotoMap:
    const PixelMap* map() const {return pm;}
    const double refColor() const {return reference;}
    
    virtual PixelMap* duplicate() const {return new ColorTerm(pm->duplicate(),
							      reference,
							      getName());}
    static string type() {return "Color";}
    virtual string getType() const {return type();}
    virtual bool needsColor() const {return true;}
    // Will not have Create() call since color maps will be built by PixelMapCollection.
    // But for some convenience, provide a separate serialization.
#ifdef USE_YAML
    virtual void write(YAML::Emitter& os) const;
#endif
    virtual void setParams(const DVector& p) {pm->setParams(p);}
    virtual DVector getParams() const {return pm->getParams();}
    virtual int nParams() const {return pm->nParams();}
    virtual double getPixelStep() const {return pm->getPixelStep();}
    virtual void setPixelStep(double ps) {pm->setPixelStep(ps);}

    // Implement the interface with color arguments:
    virtual void toWorld(double xpix, double ypix, 
			 double& xworld, double& yworld,
			 double color=astrometry::NODATA) const {
      // Here is the definition of the color term action:
      checkColor(color);  // Exception if we enter with default
      double c = color-reference; 
      pm->toWorld(xpix,ypix,xworld,yworld);
      xworld = xpix + c*(xworld-xpix);
      yworld = ypix + c*(yworld-ypix);
    }
    virtual void toPix( double xworld, double yworld,
			double &xpix, double &ypix,
			double color=astrometry::NODATA) const;
    virtual Matrix22 dPixdWorld(double xworld, double yworld, 
				double color=astrometry::NODATA) const;
    virtual Matrix22 dWorlddPix(double xpix, double ypix, 
				double color=astrometry::NODATA) const;
    virtual void toWorldDerivs(double xpix, double ypix,
			       double& xworld, double& yworld,
			       DMatrix& derivs,
			       double color=astrometry::NODATA) const;
    virtual void toPixDerivs( double xworld, double yworld,
			      double &xpix, double &ypix,
			      DMatrix& derivs,
			      double color=astrometry::NODATA) const;
    
private:
    PixelMap* pm;
    double reference;	// Color value yielding zero correction
    /*Hide - do not want confused ownership of pm*/
    ColorTerm(const ColorTerm& rhs) = delete;
    void operator=(const ColorTerm& rhs) = delete;
  
  };
} // namespace astrometry
#endif //PIXMAP_H
