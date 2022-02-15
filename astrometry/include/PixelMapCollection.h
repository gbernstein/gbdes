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
#include <vector>

#include "Std.h"
#include "PixelMap.h"
#include "Wcs.h"
// Note that these classes won't be very useful without YAML for serialization, but
// we do want to be able to compile the repo without it.  No PixelMapCollection at all.
#ifdef USE_YAML
#include "yaml-cpp/yaml.h"
#endif

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

    bool anyColor;	// True if any of the maps require color

    bool ownMaps;		// true if we must delete maps on destruction
    int totalFreeParameters;	// Cache this total
    void countFreeParameters(); // update countFreeParameters from vNSubParams

    // Hide to avoid inadvertent ownership confusion:
    SubMap(const SubMap& rhs);
    void operator=(const SubMap& rhs);

  public:
    // If shareMaps = false, SubMap will make and own duplicates of the input maps.
    // Otherwise it will just copy pointers and assume external ownership.
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
    int nParams() const {return totalFreeParameters;}

    void toWorld(double xpix, double ypix,
		 double& xworld, double& yworld,
		 double color=astrometry::NODATA) const;
    void toPix( double xworld, double yworld,
		double &xpix, double &ypix,
		double color=astrometry::NODATA) const;
    void toPixDerivs( double xworld, double yworld,
		      double &xpix, double &ypix,
		      DMatrix& derivs,
		      double color=astrometry::NODATA) const;
    void toWorldDerivs(double xpix, double ypix,
		       double& xworld, double& yworld,
		       DMatrix& derivs,
		       double color=astrometry::NODATA) const;
    Matrix22 dPixdWorld(double xworld, double yworld, 
			double color=astrometry::NODATA) const;
    Matrix22 dWorlddPix(double xpix, double ypix, 
			double color=astrometry::NODATA) const;

    // Step size for derivatives in "pixel" space - applied to first map in the chain, if any
    virtual double getPixelStep() const;
    virtual void setPixelStep(double ps);
    static string type() {return "Composite";}
    virtual string getType() const {return type();}

    virtual bool needsColor() const {return anyColor;}
    
#ifdef USE_YAML
    void write(YAML::Emitter& os) const;
#endif
  };
  
#ifdef USE_YAML
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

    // One can shut off the rebuilding of indices since it can be slow.
    void learnMap(const PixelMap& pm, bool duplicateNamesAreExceptions=false,
		  bool rebuildIndices=true); 
    // Note that if you pass a Wcs to learnMap, it will be treated as a map followed by
    // a reprojection, and will not be available to be issued as a Wcs.
    void learnWcs(const Wcs& pm, bool duplicateNamesAreExceptions=false,
		  bool rebuildIndices=true); 
    void learn(PixelMapCollection& source, bool duplicateNamesAreExceptions=false);

    // Mark a WCS as invalid / unused
    void invalidate(string wcsName);

    // Purge all map elements that are only used by invalidated WCS's
    void purgeInvalid();
      
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

    // These two are like the above, except that they return objects that own themselves and
    // all subcomponents, so this PixelMapCollection can be destroyed while cloned map/wcs remain
    // useful.  All parameters of the cloned objects are free though, and disconnected from the
    // PixelMapCollection's master parameter vector.
    PixelMap* cloneMap(string mapName) const;
    Wcs* cloneWcs(string wcsName) const;
    
    // Reassign all parameter indices and update all issued SubMaps (call whenever
    // setFixed or adding new elements with parameters)
    void rebuildParameterVector();

    // Create and read serializations of all the maps or just one
    void write(ostream& os, string comment="") const;
    void writeMap(ostream& os, string name, string comment="") const;
    void writeWcs(ostream& os, string name, string comment="") const;
    // The read returns false if the stream does not begin with the magicWord for
    // serialized PixelMapCollections. True if success; throws exception for format errors.
    bool read(istream& is, string namePrefix=""); // prefix added to names of all elements
    
    // Fix or free parameters associated with one or more PixelMap(s).
    // If the chosen map is compound, all components are frozen.
    // This will change parameter index assignments in all issued SubMaps and Wcs's
    // Anything with zero free parameters is always considered fixed
    void setFixed(set<string> nameList, bool isFixed=true);
    void setFixed(string name, bool isFixed=true);
    bool getFixed(string name) const;

    // Get whether a map has any parameters that are still at defaults,
    // i.e. parameters have not yet been set or fit to any data
    bool getDefaulted(string name) const;

    // Tell whether a map is atomic
    bool isAtomic(string name) const;

    // Set/get the master parameter vector for all PixelMaps
    void setParams(const DVector& p);
    DVector getParams() const;
    int nParams() const {return parameterCount;}

    // Set parameters of a member map by copying from
    // an input PixelMap of the same name.  Exception if
    // map does not exist or is fixed.
    void copyParamsFrom(const PixelMap& pm);
    
    // Some access to members of the collection of maps:
    int nMaps() const {return mapElements.size();}
    // Number of atomic maps in the collection:
    int nAtomicMaps() const {return atomCount;}
    // Number of maps with free parameters:
    int nFreeMaps() const {return freeCount;}
    bool mapExists(string name) const {return mapElements.count(name);}
    vector<string> allMapNames() const;

    // And the WCS specified for this collection:
    int nWcs() const {return wcsElements.size();}
    bool wcsExists(string name) const {return wcsElements.count(name);}
    vector<string> allWcsNames() const;

    // Return whether the map given by first name depends on
    // a map of the second name.
    bool dependsOn(const string mapName, const string targetName) const;

    // Return list of names of all map elements needed to fully specify the target.
    // Includes self.  Assumes no dependence cycles.
    set<string> dependencies(string mapName) const;


    // This is a routine useful for debugging: return the name of the atomic
    // map that a certain parameter in the vector belongs to.
    string atomHavingParameter(int parameterIndex) const;
    // And a convenience to get the indices of a map's params without
    // having to issue the map.  Returns 0's if map is fixed or compound.
    void parameterIndicesOf(string mapname, int& startIndex, int& nParams) const;
    

    // This routine adds a new type of PixelMap to the parsing dictionary
    template <class MapType>
    static void registerMapType() {
      mapCreators[MapType::type()] = MapType::create;
    }

    // This is a keyword that must appear in a YAML deserialization to
    // indicate a valid PixelMapCollection
    static const string magicKey;
    
  private:
    // hide:
    PixelMapCollection(const PixelMapCollection& rhs);
    void operator=(const PixelMapCollection& rhs);

    int parameterCount; // Total parameters currently free to vary.
    int atomCount;	// Number of atomic map components
    int freeCount;	// Number of atomic components with free parameters

    // ***Main data of the class are these two containers of structures listing
    // all the PixelMaps and Wcs's curated by this class:

    // Structure for every PixelMap that we know about:
    struct MapElement {
      MapElement(): realization(0), atom(0), isFixed(false), isDefaulted(false), number(-1) {}
      list<string> subordinateMaps;  // If it's compound, what it will be made from
      SubMap* realization;	     // pointer to its SubMap, if it's been built
      PixelMap* atom;		     // Pointer to the PixelMap itself, if atomic
      int startIndex;		     // Location in the union parameter array
      int nParams;		     // Number of parameters (only atomic is nonzero)
      int number;		     // sequential index among all maps with free parameters
      bool isFixed;		     // True if atom has parameters currently fixed.
      bool isDefaulted;		     // True if atom has default parameters
    };

    map<string, MapElement> mapElements;  // all known PixelMaps, indexed by name.

    // Structure for every WCS that we know about:
    struct WcsElement {
      WcsElement(): nativeCoords(0), realization(0), wScale(DEGREE) {}
      WcsElement(const WcsElement& rhs): mapName(rhs.mapName), nativeCoords(0),
					 realization(rhs.realization),
					 wScale(rhs.wScale),
					 isValid(true) {
	if (rhs.nativeCoords) nativeCoords = rhs.nativeCoords->duplicate();
      }
      ~WcsElement() {if (nativeCoords) delete nativeCoords;}
      string mapName;	// The PixelMap it should use.
      // Projection it will use.  We own it.
      SphericalCoords* nativeCoords;
      Wcs* realization;	// Pointer to realization if we've got it.
      double wScale;	// Scaling to apply to PixelMap world coords to get radians
      bool isValid;
    private:
      void operator=(const MapElement& rhs); // hide assignment
    };

    map<string, WcsElement> wcsElements; // all known Wcs's, indexed by name.

    // Remove knowledge of a particular map or WCS and its realization if it exists.
    void removeMap(string mapName);
    void removeWcs(string wcsName);  // Does NOT remove the map on which it's based.

    // **** Useful utilities: *****

    // See if the chain of dependence of a map has any cycles.  Throws exception if so.
    void checkCircularDependence(string mapName,
				 const set<string>& ancestors = set<string>()) const;

    // Check that all referenced names exist, and that there are no circular dependences.
    void checkCompleteness() const;

    // Produce a list giving the atomic transformation sequence needed to implement a
    // specified map. Assumes no dependence cycles.
    list<string> orderAtoms(string mapName) const;

    // **** Static structures / methods for serialization: ****
    // This routine will write a complete YAML key/value to emitter for this PhotoMap
    void writeSingleMap(YAML::Emitter& os, const MapElement& mel, string name) const;
    void writeSingleWcs(YAML::Emitter& os, const WcsElement& wel, string name) const;
    PixelMap* createAtomFromNode(const YAML::Node& node,
				 bool& defaulted,
				 string name) const;
    // Signature of routines to create maps from yaml data; defaulted if
    // parameters were not given explicitly.
    typedef PixelMap* (*Creator)(const YAML::Node& node, bool& defaulted,
				 string name);
    // Static map of what kind of PixelMaps the parser will be able to identify
    static map<string, Creator> mapCreators;
    static bool creatorsInitialized;
    // Call to give parser an initial inventory of map types to create:
    static void PixelMapTypeInitialize();
  };

#endif  // end USE_YAML

} // namespace astrometry
#endif
