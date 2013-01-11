// $Id: PixelMapCollection.h,v 1.2 2012/02/03 19:20:11 garyb Exp $
// MapCollection is a storage area for many PixelMap components that are going to be 
// chained together in different combinations to consitute the astrometric maps
// for different regions of a dataset.
//    Use MapCollection to store all the components and produce all the CompoundMaps
// that represent the chains; this class will then do the bookkeeping and destruction of
// them all.
//    SubMap is a derived class of PixelMap that wraps another PixelMap but maps its 
// parameters into a larger parameter vector containing many maps' parameters.

#ifndef MAPCOLLECTION_H
#define MAPCOLLECTION_H

#include <list>
using std::list;
#include <string>
#include "Std.h"
#include "UseTMV.h"
#include "PixelMap.h"
#include "NameIndex.h"

namespace astrometry {

  class SubMap;

  typedef int PixelMapKey;
  class PixelMapChain: public list<PixelMapKey> {
  public:
    // will use front(), size(), empty(), iterators, begin(), end()
    void append(PixelMapKey m) {push_back(m);}
    void append(PixelMapChain& mc) {insert(end(), mc.begin(), mc.end());}
  };

  class PixelMapCollection {
  public:
    PixelMapCollection(): parameterCount(0) {}
    ~PixelMapCollection();
    PixelMapKey add(PixelMap* pm, string name=""); // Assign a name if one is not given.
    SubMap* issue(PixelMapChain& chain);
    void setParams(const DVector& p);
    DVector getParams() const;
    int nParams() const {return parameterCount;}
    int nMaps() const {return members.size();}
    PixelMapKey indexOf(string name) const {return names.indexOf(name);}
    string nameOf(PixelMapKey m) const {return names.nameOf(m);}
    PixelMap* element(PixelMapKey m) {return members[m];}
    const PixelMap* element(PixelMapKey m) const {return members[m];}
    int startIndex(PixelMapKey k) const {return startIndices[k];}
    int nParams(PixelMapKey k) const {return members[k]->nParams();}
  private:
    // hide:
    PixelMapCollection(const PixelMapCollection& rhs) {}
    void operator=(const PixelMapCollection& rhs) {}
    vector<PixelMap*> members;
    vector<int> startIndices;
    int parameterCount;
    list<CompoundPixelMap*> createdCompoundMaps;
    list<SubMap*> createdSubMaps;
    NameIndex names;
  };

  class SubMap: public PixelMap {
  public:
    SubMap(PixelMap* pm): map(pm) {}
    PixelMap* map;
    // These two arrays are all that we are adding to PixelMap:
    vector<int> startIndices;
    vector<int> nSubParams;
    PixelMapChain chain;
    // All the PixelMap methods are forwarded to the referred PixelMap:
    void setParams(const DVector& p) {map->setParams(p);}
    DVector getParams() const {return map->getParams();}
    int nParams() const {return map->nParams();}
    void toPixDerivs( double xworld, double yworld,
		      double &xpix, double &ypix,
		      DMatrix& derivs) const {
      map->toPixDerivs(xworld,yworld,xpix,ypix,derivs);
    }
    void toWorldDerivs(double xpix, double ypix,
		       double& xworld, double& yworld,
		       DMatrix& derivs) const {
      map->toWorldDerivs(xpix,ypix,xworld,yworld,derivs);
    }
    void toWorld(double xpix, double ypix,
		 double& xworld, double& yworld) const {
      map->toWorld(xpix,ypix,xworld,yworld);
    }
    void toPix( double xworld, double yworld,
		double &xpix, double &ypix) const {
      map->toPix(xworld,yworld,xpix,ypix);
    }
    Matrix22 dPixdWorld(double xworld, double yworld) const {
      return map->dPixdWorld(xworld,yworld);
    }
    Matrix22 dWorlddPix(double xpix, double ypix) const {
      return map->dWorlddPix(xpix,ypix);
    }
    // Step size for derivatives in "pixel" space:
    double getPixelStep() const {return map->getPixelStep();}
    void setPixelStep(double ps) {map->setPixelStep(ps);}
  };

} // namespace astrometry
#endif
