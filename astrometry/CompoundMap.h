// DEPRECATED compound map class - prefer SubMap
#ifndef COMPOUND_MAP_H
#define COMPOUND_MAP_H

#include <list>
#include "PixelMap.h"

namespace astrometry {
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

    const list<PixelMap*>& getList() const {return pmlist;}

    static string mapType() {return "Compound";}
    virtual string getType() const {return mapType();}
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

    virtual void write(std::ostream& os) const {
      throw AstrometryError("write() not implemented for CompoundPixelMap, name is " + getName());
    }
  };
} // end namespace astrometry

#endif
