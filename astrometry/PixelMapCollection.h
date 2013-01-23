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
    //    Calls to get/setParams work with a vector that is the concatenated
    // parameter vectors of the chain elements, in order that they are chained.
    // If, however, the vNSubParams entry for any referenced PixelMap is set to zero, it
    // will be treated as a map with no free parameters, and skipped in the chain.
    friend class PixelMapCollection;
    // Since SubMaps can be issued by a PixelMapCollection, it is a friend class that
    // can fix/free the parameters of an element by accessing vNSubParams, and also
    // adjusts the vStartIndices to reflect each element's position in a master parameter
    // vector.  
  private:
    // The constituent maps and their indices into a master parameter vector,
    // plus an ordinal label in the full collection for each map component:
    vector<PixelMap*> vMaps;
    vector<int> vStartIndices;
    vector<int> vNSubParams;
    vector<int> vMapNumbers;  

    bool ownMaps;		// true if we must delete maps on destruction
    int totalFreeParameters;	// Cache this total
    void countFreeParameters(); // update countFreeParameters from vNSubParams

    // Hide to avoid inadvertent ownership confusion:
    SubMap(const SubMap& rhs);
    void operator=(const SubMap& rhs);

  public:
    // If shareMaps = true, SubMap will make and own duplicates of the input maps.
    // Otherwise is will just copy pointers and assume external ownership.
    SubMap(const list<PixelMap*>& pixelMaps, string name="", bool shareMaps=false); 
    virtual ~SubMap();
    // A duplicate of SubMap will have the same sharing policies as its parent.
    // Also, any information on start indices / nSubParams is lost in duplicate.
    virtual PixelMap* duplicate() const;

    const PixelMap* getMap(int i) const {return vMaps[i];}
    int nMaps() const {return vMaps.size();}
    int startIndex(int iMap) const {return vStartIndices[iMap];}
    int nSubParams(int iMap) const {return vNSubParams[iMap];}
    int mapNumber(int iMap) const {return vMapNumbers[iMap];}

    // Implement the full PixelMap interface; parameter getting/setting
    // works with a single contiguous vector of parameters of all components.
    void setParams(const DVector& p);
    DVector getParams() const;
    int nParams() const {return totalFreeParameters;} // ??? need nSubParams too ???
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
    void write(ostream& os) const {
      throw AstrometryError("SubMap " + getName() + " should not be getting serialized");
    }
  };

  // ?? Implement the PixelMap and/or Wcs interfaces by selecting a member to use??

  class PixelMapCollection {
  public:
    // add some or all maps from another collection
    // add wcs from another collection or external
    PixelMapCollection();
    ~PixelMapCollection();

    // The following routines have the PixelMapCollection learn the names and behavior
    // of individual PixelMaps, Wcs's, or an entire other collection.  This Collection
    // will make copies of everything needed to issue new SubMaps or Wcs's that reproduce
    // the behavior of the contributed items.  You just need to request something issued
    // under the same name as the input object(s).

    // If called with duplicateNamesAreExceptions=true, then incoming names that match 
    // any name already in the collection will 
    // throw an exception.  If the flag is false, any new PixelMap or Wcs with same name as
    // an existing will be assumed to have the same behavior as the existing one, and will
    // therefore just be ignored.

    void learnMap(const PixelMap& pm, bool duplicateNamesAreExceptions=false); 
    // Note that if you pass a Wcs to learnMap, it will be treated as a map followed by
    // a reprojection, and will not be available to be issued as a Wcs.
    void learnWcs(const Wcs& pm, bool duplicateNamesAreExceptions=false);
    void learn(PixelMapCollection& source, bool duplicateNamesAreExceptions=false);

    // Define a new pixelMap that is compounding of a list of other PixelMaps.  Order
    // of the list is from pixel to world coordinates.
    void defineChain(string chainName, const list<string>& elements);
    // Define a WCS system to be the given PixelMap followed by projection to the
    // sky described by the nativeCoords system.
    void defineWcs(string wcsName, const SphericalCoords& nativeCoords, string mapName,
		   double wScale = DEGREE);

    // Return pointer to a SubMap realizing the named coordinate transformation
    SubMap* issueMap(string mapName);
    // Return a pointer to a Wcs built from a SubMap and realizing the named coord system
    Wcs* issueWcs(string wcsName);
    
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
    bool wcsExists(string name) const {return wcsElements.count(name);}
    vector<string> allWcsNames() const;

    // This routine adds a new type of PixelMap to the parsing dictionary
    template <class MapType>
    static void registerMapType();

  private:
    // hide:
    PixelMapCollection(const PixelMapCollection& rhs);
    void operator=(const PixelMapCollection& rhs);

    int parameterCount; // Total parameters currently free to vary.

    // ***Main data of the class are these two containers of structures listing
    // all the PixelMaps and Wcs's curated by this class:

    // Structure for every PixelMap that we know about:
    struct MapElement {
      MapElement(): realization(0), atom(0), isFixed(false), number(-1) {}
      list<string> subordinateMaps;  // If it's compound, what it will be made from
      SubMap* realization;	     // pointer to its SubMap, if it's been built
      PixelMap* atom;		     // Pointer to the PixelMap itself, if atomic
      int startIndex;		     // Location in the union parameter array
      int nParams;		     // Number of parameters (only atomic is nonzero)
      int number;		     // sequential index among all maps with free parameters
      bool isFixed;		     // True if its parameters are currently fixed.
    };

    map<string, MapElement> mapElements;  // all known PixelMaps, indexed by name.

    // Structure for every WCS that we know about:
    struct WcsElement {
      WcsElement(): nativeCoords(0), realization(0), wScale(DEGREE) {}
      WcsElement(const WcsElement& rhs): mapName(rhs.mapName), nativeCoords(0),
					 realization(rhs.realization),
					 wScale(rhs.wScale)    {
	if (rhs.nativeCoords) nativeCoords = rhs.nativeCoords->duplicate();
      }
      ~WcsElement() {if (nativeCoords) delete nativeCoords;}
      string mapName;	// The PixelMap it should use.
      // Projection it will use.  We own it.
      SphericalCoords* nativeCoords;
      Wcs* realization;	// Pointer to realization if we've got it.
      double wScale;	// Scaling to apply to PixelMap world coords to get radians
    private:
      void operator=(const MapElement& rhs); // hide assignment
    };

    map<string, WcsElement> wcsElements; // all known Wcs's, indexed by name.

    typedef map<string, MapElement>::iterator MapIter;
    typedef map<string, WcsElement>::iterator WcsIter;
    typedef map<string, MapElement>::const_iterator ConstMapIter;
    typedef map<string, WcsElement>::const_iterator ConstWcsIter;

    // **** Useful utilities: *****

    // Reassign all parameter indices and update all issued SubMaps (call whenever
    // setFixed or adding new elements with parameters)
    void rebuildParameterVector();

    // See if the chain of dependence of a map has any cycles.  Throws exception if so.
    void checkCircularDependence(string mapName,
				 const set<string>& ancestors = set<string>()) const;

    // Return list of pointers to all map elements needed to fully specify the target.
    // Includes self.  Assumes no dependence cycles.
    set<string> dependencies(string mapName) const;

    // Produce a list giving the atomic transformation sequence needed to implement a
    // specified map. Assumes no dependence cycles.
    list<string> orderAtoms(string mapName) const;

    // **** Static structures for serialization: ****
    // Static lists of what kind of PixelMaps the parser will be able to identify
    static vector<string> mapTypeNames;
    typedef PixelMap* (*Creator)(istream& is, string name);
    static vector<Creator> mapCreators;
    static bool creatorsInitialized;
    // Call to give parser an initial inventory of map types to create:
    static void PixelMapTypeInitialize();

  };


} // namespace astrometry
#endif
