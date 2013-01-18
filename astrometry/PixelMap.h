// PixelMap is a base class for maps from 2d "pixel" space to 2d "world" space.
// Maps should be bijective over useful range as inverse should be defined.
// Designed to be fit to data by adjusting parameters of the map.
// pixelStep is the increment in pixel space that should be used for taking numerical derivatives.
//   For purposes of serialization, each derived class should define static methods
//   string Derived::mapType() giving a unique name to this class.
//   PixelMap* create(istream& is, string name)  which produces member of the class from stream input
//   ... and implement the virtual method void write(ostream& os).
//   ... and also every PixelMap has a string name() bestowed on creation.  If name is not given, 
//   a name will be created from running counter.
#ifndef PIXELMAP_H
#define PIXELMAP_H

#include <string>
#include <list>
#include <iostream>
#include "Std.h"
#include "UseTMV.h"
#include "Astrometry.h"

// Note that Astrometry.h has typdefs for the 2- and 3-element small vectors & matrices

namespace astrometry {

  class PixelMap {
  public:
    PixelMap(string name_="");
    virtual ~PixelMap() {}
    // Return pointer to deep copy of self
    virtual PixelMap* duplicate() const =0;

    // pixel coords to world coords map
    virtual void toWorld(double xpix, double ypix,
			 double& xworld, double& yworld) const=0;
    // inverse map:
    virtual void toPix( double xworld, double yworld,
			double &xpix, double &ypix) const=0;
    virtual Matrix22 dPixdWorld(double xworld, double yworld) const;
    virtual Matrix22 dWorlddPix(double xpix, double ypix) const;
    // Solid angle of a pixel (sr)
    double pixelArea(double xpix, double ypix) const;

    // Get parameters of the mapping.  Note default behavior of derived classes
    // is to ignore anything to do with map parameters.
    virtual void setParams(const DVector& p) {}
    // Note here that a zero-length vector is returned as parameters by default
    virtual DVector getParams() const {return DVector(0);}
    virtual int nParams() const {return 0;}
    // These calls include derivs of world w.r.t. parameters of map.  
    // Default implementation does not calculate derivs as default is no parameters
    virtual void toWorldDerivs(double xpix, double ypix,
			       double& xworld, double& yworld,
			       DMatrix& derivs) const {
      toWorld(xpix,ypix,xworld,yworld);}
    // In toPixDerivs, the derivatives are of world coords w.r.t. map parameters, even
    // thought the call is going from world to pix coords:
    virtual void toPixDerivs( double xworld, double yworld,
			      double &xpix, double &ypix,
			      DMatrix& derivs) const {
      toPix(xworld,yworld,xpix,ypix);}
    // Step size for derivatives in "pixel" space:
    virtual double getPixelStep() const {return pixStep;}
    virtual void setPixelStep(double ps) {pixStep=ps;}

    // For serialization:
    virtual void write(std::ostream& os) const =0;
    string getName() const {return name;}
  private:
    string name;
    double pixStep;
    static int anonymousCounter;
  protected:
    // Invert the pix->world map using Newton iteration.
    // Likely that many derived PixelMap classes will want this.
    void NewtonInverse(double xworld, double yworld, 
		       double& xpix, double& pix,
		       double worldTolerance) const;
  };

  inline std::ostream& operator<<(std::ostream& os, const PixelMap& pm) {pm.write(os); return os;}

  class CompoundPixelMap: public PixelMap {
  private:
    // PixelMap whose pix - to - world map follows series of transformations in order of list
    // Keeps pointers to elements of compound map; destruction is up to user.
    // Note that a deep copy using duplicate() will make cleanup difficult!
    // This class is deprecated; prefer PixelMapCollection::SubMap

    list<PixelMap*> pmlist;
  public:
    CompoundPixelMap(PixelMap* pm1, string name_=""): PixelMap(name_) {
      pmlist.push_back(pm1);
    }
    ~CompoundPixelMap() {}
    virtual PixelMap* duplicate() const;

    static string mapType() {return "Compound";}
    // * will not implement a static read from stream since compound must be built in code that is aware of 
    // * the elements to be compounded.
    void append(PixelMap* pmnew) {pmlist.push_back(pmnew);}
    void prepend(PixelMap* pmnew) {pmlist.push_front(pmnew);}
    // pixel coords to world coords map
    virtual void toWorld(double xpix, double ypix,
			 double& xworld, double& yworld) const;
    // inverse map:
    virtual void toPix( double xworld, double yworld,
			double &xpix, double &ypix) const;

    // Here the derivatives are of world coords w.r.t. map parameters
    virtual void toPixDerivs( double xworld, double yworld,
			      double &xpix, double &ypix,
			      DMatrix& derivs) const;
    // And derivs of world w.r.t. parameters of map.  Resizing of derivs
    // done if needed.
    virtual void toWorldDerivs(double xpix, double ypix,
			       double& xworld, double& yworld,
			       DMatrix& derivs) const;

    // Get parameters of the mapping.  Note default behavior of derived classes
    // is to ignore anything to do with map parameters.
    virtual void setParams(const DVector& p);
    virtual DVector getParams() const;
    virtual int nParams() const;

    virtual double getPixelStep() const {return pmlist.front()->getPixelStep();}
    virtual void setPixelStep(double p) {pmlist.front()->setPixelStep(p);}
    // ??? should set all steps by propagating derivatives??

    virtual void write(std::ostream& os) const;
  };

  class IdentityMap: public PixelMap {
  public:
    IdentityMap(): PixelMap("Identity") {}
    virtual PixelMap* duplicate() const {return new IdentityMap();}
    static string mapType() {return "Identity";}
    static PixelMap* create(std::istream& is, string name_="") {return new IdentityMap;}
    void toWorld(double xpix, double ypix,
		 double& xworld, double& yworld) const;
    void toPix( double xworld, double yworld,
		double &xpix, double &ypix) const;
    Matrix22 dPixdWorld(double xworld, double yworld) const;    
    Matrix22 dWorlddPix(double xpix, double ypix) const;    
    virtual void write(std::ostream& os) const {} // Nothing to write
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
    // setPixelScale to 1 arcsecond for any numerical derivatives of this projection,
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

    static string mapType() {return "Reprojection";}
    static PixelMap* create(std::istream& is, string name="");
    virtual void write(std::ostream& os) const;

    void toWorld(double xpix, double ypix,
		 double& xworld, double& yworld) const;
    void toPix( double xworld, double yworld,
		double &xpix, double &ypix) const;
    Matrix22 dPixdWorld(double xworld, double yworld) const;    
    Matrix22 dWorlddPix(double xpix, double ypix) const;    

  private:
    SphericalCoords* pix;
    SphericalCoords* world;
    double scaleFactor;
  };

} // namespace astrometry
#endif //PIXMAP_H
