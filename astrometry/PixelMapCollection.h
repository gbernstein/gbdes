// PixelMapCollection is a storage area for many PixelMap components that are going to be 
// chained together in different combinations to consitute the astrometric maps
// for different regions of a dataset.  It also stores and creates Wcs's from these maps.
//    Use PixelMapCollection to store all the components and produce all the CompoundMaps
// that represent the chains; this class will then do the bookkeeping and destruction of
// them all.
//    PixelMapCollection also can serialize itself to a stream or reconstruct maps
// from a stream.
//    SubMap is a derived class of PixelMap that wraps another PixelMap but augments
// it with information on where the PixelMap's parameters fit in a vector that
// is the union of all PixelMap's parameters.
//    Parameters of any component map can be frozen or un-frozen for fitting purposes.
// 
// PixelMaps issued by PixelMapCollection are all SubMaps.  They contain reference
// to the parent PixelMapCollection and are aware of what elements have been frozen.

#ifndef MAPCOLLECTION_H
#define MAPCOLLECTION_H

#include <list>
#include <map>
#include <set>
using std::list;
#include <string>
#include "Std.h"
#include "UseTMV.h"
#include "PixelMap.h"
#include "Wcs.h"

namespace astrometry {

  class SubMap: public PixelMap {
    // A PixelMap that is the compounded action of one or more PixelMaps potentially
    // taken from a PixelMapCollection.  If so, the
    // startIndices and nParams index each compounded map's parameters into a master vector.
    friend class PixelMapCollection;
  private:
    // These two arrays are all that we are adding to PixelMap:
    vector<int> vStartIndices;
    vector<int> vNSubParams;
    int totalFreeParameters;	// Cache this total
    void countFreeParameters(); // update countFreeParameters from vNSubParams
  public:
    // ?? names ??
    SubMap();
    SubMap(list<PixelMap*> pixelMaps); 
    virtual ~SubMap();
    vector<PixelMap*> vMaps;
    int nMaps() const {return vMaps.size();}
    int startIndex(int iMap) const {return vStartIndices[iMap];}
    int nSubParams(int iMap) const {return vNSubParams[iMap];}

    // Implement the full PixelMap interface; parameter getting/setting
    // works with a single contiguous vector of parameters of all components.
    void setParams(const DVector& p);
    DVector getParams() const;
    int nParams() const {return totalFreeParameters;}
    void toWorld(double xpix, double ypix,
		 double& xworld, double& yworld) const;
    void toPix( double xworld, double yworld,
		double &xpix, double &ypix) const;
    void toPixDerivs( double xworld, double yworld,
		      double &xpix, double &ypix,
		      DMatrix& derivs) const;
    void toWorldDerivs(double xpix, double ypix,
		       double& xworld, double& yworld,
		       DMatrix& derivs) const;
    Matrix22 dPixdWorld(double xworld, double yworld) const;
    Matrix22 dWorlddPix(double xpix, double ypix) const;
    // Step size for derivatives in "pixel" space - applied to first map in the chain, if any
    double getPixelStep() const;
    void setPixelStep(double ps);
    void write(ostream& os) const;
  };

  // Convenience class to hold sequences of PixelMaps (by names)
  class PixelMapChain: public list<string> {
  public:
    // will use front(), size(), empty(), iterators, begin(), end()
    void append(string m) {push_back(m);}
    void append(PixelMapChain& mc) {insert(end(), mc.begin(), mc.end());}
  };

  // ?? Implement the PixelMap and/or Wcs interfaces by selecting a member to use??
  // get rid of Keys and Chains and just use names??
  class PixelMapCollection {
  public:
    // add some or all maps from another collection
    // add wcs from another collection or external
    // create the last wcs or pixmap specified in some file??
    PixelMapCollection();
    ~PixelMapCollection();

    // Include a new map or Wcs into this collection (will also absorb any parts
    // of compound maps).  PixelMapCollection will assume ownership.
    // Duplicate names are left as is when adding.
    void addMap(PixelMap* pm); 
    void addWcs(Wcs* pm); 

    // Get pointer to a new map constructed as the specified chain.
    SubMap* issueMap(PixelMapChain& chain, string name="");
    // Return pointer to map with specified name.  Returns 0 if name is not in collection
    SubMap* issueMap(string name);

    // Get pointer to a new Wcs constructed as the specified chain and 
    // projected with the given Orientation or projection.  Name is given to
    // both the PixelMap and the Wcs that uses it.
    Wcs* issueWcs(PixelMapChain& chain, const Orientation& nativeOrient, string name="");
    Wcs* issueWcs(PixelMapChain& chain, const SphericalCoords& nativeCoords, string name="");
    // Issue a WCS that is projection of existing PixelMap into specified coord system
    Wcs* issueWcs(string mapName, const Orientation& nativeOrient, string name="");
    Wcs* issueWcs(string mapName, const SphericalCoords& nativeCoords, string name="");
    // Issue a previously specified WCS
    Wcs* issueWcs(string name);
    
    // Create and read serializations of all the maps or just one
    void read(istream& is, string namePrefix=""); // prefix added to names of all elements
    void write(ostream& os) const;
    void writeMap(ostream& os, string name) const;
    void writeWcs(ostream& os, string name) const;

    // Fix or free parameters associated with one or more PixelMap(s).
    // If the chosen map is compound, all components are frozen.
    // This will change parameter index assignments in all issued SubMaps and Wcs's
    void setFixed(list<string> nameList, bool isFixed);
    void setFixed(string name, bool isFixed) {setFixed(list<string>(1,name), isFixed);}
    bool getFixed(string name) const;

    // Set/get the master parameter vector for all PixelMaps
    void setParams(const DVector& p);
    DVector getParams() const;
    int nParams() const {return parameterCount;}

    // Some access to members of the collection of maps:
    int nMaps() const {return mapElements.size();}
    bool mapExists(string name) const {return mapElements.count(name);}
    vector<string> allMapNames() const;

    // And the WCS specified for this collection:
    int nWcs() const {return wcsElements.size();}
    bool WcsExists(string name) const {return wcsElements.count(name);}
    vector<string> allWcsNames() const;

    // This routine adds a new type of PixelMap to the parsing dictionary
    template <class MapType>
    static void registerMapType();

  private:
    // hide:
    PixelMapCollection(const PixelMapCollection& rhs);
    void operator=(const PixelMapCollection& rhs);

    // Static lists of what kind of PixelMaps the parser will be able to identify
    static vector<string> mapTypeNames;
    typedef PixelMap* (*Creator)(istream& is, string name);
    static vector<Creator> mapCreators;
    static bool creatorsInitialized;
    // Call to give parser an initial inventory of map types to create:
    static void PixelMapTypeInitialize();

    // Structure for every PixelMap that we know about:
    struct MapElement {
      MapElement(): realization(0), atom(0) {}
      list<string> subordinateMaps;  // If it's compound, what it will be made from
      SubMap* realization;	     // pointer to its SubMap, if it's been built
      PixelMap* atom;		     // Pointer to the PixelMap itself, if atomic
      int startIndex;		     // Location in the union parameter array
      int nParams;		     // Number of parameters (only atomic is nonzero)
      bool isFixed;		     // True if its parameters are currently fixed.
    };

    map<string, MapElement> mapElements;  // all known PixelMaps, indexed by name.

    // Structure for every WCS that we know about:
    struct WcsElement {
      WcsElement(): nativeCoords(0), realization(0) {}
      WcsElement(const WcsElement& rhs): mapName(rhs.mapName), nativeCoords(0),
					 realization(rhs.realization) {
	if (rhs.nativeCoords) nativeCoords = rhs.nativeCoords->duplicate();
      }
      string mapName;	// The PixelMap it should use.
      // Projection it will use.  We own it.
      SphericalCoords* nativeCoords;
      Wcs* realization;	// Pointer to realization if we've got it.
    private:
      void operator=(const MapElement& rhs); // hide assignment
    };

    map<string, WcsElement> wcsElements; // all known Wcs's, indexed by name.
    typedef map<string, MapElement>::iterator MapIter;
    typedef map<string, WcsElement>::iterator WcsIter;
    typedef map<string, MapElement>::const_iterator ConstMapIter;
    typedef map<string, WcsElement>::const_iterator ConstWcsIter;

    int parameterCount; // Total parameters currently free to vary.

    // Return list of pointers to all map elements needed to fully specify the target.
    // Includes self, after checking for circular dependence on any of the ancestors.
    set<string> dependencies(string target,
			     const set<string>& ancestors = set<string>()) const;
    // Reassign all parameter indices and update all issued SubMaps (call whenever
    // setFixed or adding new elements with parameters)
    void rebuildParameterVector();

  };


} // namespace astrometry
#endif
