// PhotoMapCollection is a storage area for many PhotoMap components that are going to be 
// chained together in different combinations to constitute the photometric maps.
// Also enables global fitting of the parameters of these maps, including the use of
// priors relating the zeropoints of different exposures.
//    Use PhotoMapCollection to store all the components and produce all the SubMaps
// that represent the chains; this class will then do the bookkeeping and destruction of
// them all.
//    PhotoMapCollection also can serialize itself to a stream or reconstruct maps
// from a stream.
//    SubMap is a derived class of PhotoMap that wraps another PhotoMap but augments
// it with information on where the PhotoMap's parameters fit in a vector that
// is the union of all PhotoMap's parameters.
//    Parameters of any component map can be frozen or un-frozen for fitting purposes.
// 
// PhotoMaps issued by PhotoMapCollection are all SubMaps.  They contain reference
// to the parent PhotoMapCollection and are aware of what elements have been frozen.

#ifndef PHOTOMAPCOLLECTION_H
#define PHOTOMAPCOLLECTION_H

#include <list>
#include <map>
#include <set>
using std::list;
#include <string>
#include "Std.h"
#include "PhotoMap.h"

#include "yaml-cpp/yaml.h"

namespace photometry {

  class SubMap: public PhotoMap {
    // A PhotoMap that is the compounded action of one or more PhotoMaps potentially
    // taken from a PhotoMapCollection.  If so, the
    // startIndices and nParams index each compounded map's parameters into a master vector.
    //    Calls to get/setParams work with a vector that is the concatenated
    // parameter vectors of the chain elements, in order that they are chained.
    // If, however, the vNSubParams entry for any referenced PhotoMap is set to zero, it
    // will be treated as a map with no free parameters, and skipped in the chain.
    friend class PhotoMapCollection;
    // Since SubMaps can be issued by a PhotoMapCollection, it is a friend class that
    // can fix/free the parameters of an element by accessing vNSubParams, and also
    // adjusts the vStartIndices to reflect each element's position in a master parameter
    // vector.  
  private:
    // The constituent maps and their indices into a master parameter vector,
    // plus an ordinal label in the full collection for each map component:
    vector<PhotoMap*> vMaps;
    vector<int> vStartIndices;
    vector<int> vNSubParams;
    vector<int> vMapNumbers;  

    bool ownMaps;		// true if we must delete maps on destruction
    int totalFreeParameters;	// Cache this total
    void countFreeParameters(); // update countFreeParameters from vNSubParams

    bool anyColor;		// True if any of the components needs color info.
    // Hide to avoid inadvertent ownership confusion:
    SubMap(const SubMap& rhs);
    void operator=(const SubMap& rhs);

  public:
    // If shareMaps = false, SubMap will make and own duplicates of the input maps.
    // Otherwise is will just copy pointers and assume external ownership.
    SubMap(const list<PhotoMap*>& photoMaps, string name="", bool shareMaps=false); 
    virtual ~SubMap();
    // A duplicate of SubMap will have the same sharing policies as its parent.
    // Also, any information on start indices / nSubParams is lost in duplicate.
    virtual PhotoMap* duplicate() const;

    const PhotoMap* getMap(int i) const {return vMaps[i];}
    int nMaps() const {return vMaps.size();}
    int startIndex(int iMap) const {return vStartIndices[iMap];}
    int nSubParams(int iMap) const {return vNSubParams[iMap];}
    int mapNumber(int iMap) const {return vMapNumbers[iMap];}

    // Implement the full PhotoMap interface; parameter getting/setting
    // works with a single contiguous vector of parameters of all components.
    void setParams(const DVector& p);
    DVector getParams() const;
    int nParams() const {return totalFreeParameters;} 

    virtual double forward(double magIn, const PhotoArguments& args) const;
    virtual double inverse(double magOut, const PhotoArguments& args) const;
    virtual double derivative(double magIn, const PhotoArguments& args) const;
    virtual double forwardDerivs(double magIn, const PhotoArguments& args,
				 DVector& derivs) const;
    virtual bool needsColor() const {return anyColor;}

    static string type() {return "Composite";}
    virtual string getType() const {return type();}

    void write(YAML::Emitter& os) const; // Provided for future convenience,
    // serialization should happen through PhotoMapCollection
  };

  class PhotoMapCollection {
  public:
    // add some or all maps from another collection
    // add wcs from another collection or external
    PhotoMapCollection();
    ~PhotoMapCollection();

    // The following routines have the PhotoMapCollection learn the names and behavior
    // of individual PhotoMaps, Wcs's, or an entire other collection.  This Collection
    // will make copies of everything needed to issue new SubMaps or Wcs's that reproduce
    // the behavior of the contributed items.  You just need to request something issued
    // under the same name as the input object(s).

    // If called with duplicateNamesAreExceptions=true, then incoming names that match 
    // any name already in the collection will 
    // throw an exception.  If the flag is false, any new PhotoMap or Wcs with same name as
    // an existing will be assumed to have the same behavior as the existing one, and will
    // therefore just be ignored.

    // One can shut off the rebuilding of indices since it can be slow.
    void learnMap(const PhotoMap& pm, bool duplicateNamesAreExceptions=false,
		  bool rebuildIndices=true); 
    void learn(PhotoMapCollection& source, bool duplicateNamesAreExceptions=false);

    // Mark a map as invalid / unused
    void invalidate(string mapName);

    // Purge all map elements that are only used by invalidated maps
    void purgeInvalid();

    // Define a new PhotoMap that is compounding of a list of other PhotoMaps.  Maps
    // are applied to mag in the order they occur in list
    void defineChain(string chainName, const list<string>& elements);

    // Return pointer to a SubMap realizing the named coordinate transformation
    SubMap* issueMap(string mapName);

    // Like the above, except that they return objects that own themselves and
    // all subcomponents, so this PhotoMapCollection can be destroyed while cloned map remains
    // useful.  All parameters of the cloned objects are free though, and disconnected from the
    // PhotoMapCollection's master parameter vector.
    PhotoMap* cloneMap(string mapName) const;
    
    // Reassign all parameter indices and update all issued SubMaps (call whenever
    // setFixed or adding new elements with parameters)
    void rebuildParameterVector();

    // Create and read serializations of all the maps or just one
    void write(ostream& os, string comment="") const;
    void writeMap(ostream& os, string name, string comment="") const;

    // The read returns false if the stream does not appear to be a
    // serialized PhotoMapCollection. True if success; throws exception for format errors.
    bool read(istream& is, string namePrefix=""); // prefix added to names of all elements

    // Fix or free parameters associated with one or more PhotoMap(s).
    // If the chosen map is compound, all components are frozen.
    // This will change parameter index assignments in all issued SubMaps and Wcs's
    void setFixed(set<string> nameList, bool isFixed=true);
    void setFixed(string name, bool isFixed=true);
    bool getFixed(string name) const;

    // Tell whether a map is atomic
    bool isAtomic(string name) const;

    // Set/get the master parameter vector for all PhotoMaps
    void setParams(const DVector& p);
    DVector getParams() const;
    int nParams() const {return parameterCount;}

    // Set parameters of a member map by copying from
    // an input PhotoMap of the same name.  Exception if
    // map does not exist or is fixed.
    void copyParamsFrom(const PhotoMap& pm);

    // Some access to members of the collection of maps:
    int nMaps() const {return mapElements.size();}
    // Number of atomic maps in the collection:
    int nAtomicMaps() const {return atomCount;}
    // Number of maps with free parameters:
    int nFreeMaps() const {return freeCount;}
    bool mapExists(string name) const {return mapElements.count(name);}
    vector<string> allMapNames() const;

    // Return list of pointers to all map elements needed to fully specify the target.
    // Includes self.  Assumes no dependence cycles.
    set<string> dependencies(string mapName) const;

    // Return whether the map given by first name depends on
    // a map of the second name.
    bool dependsOn(const string mapName, const string targetName) const;

    // Get whether a map has any parameters that are still at defaults,
    // i.e. parameters have not yet been set or fit to any data.
    // (For PhotoMaps, we are not yet tracking defaults ???)
    bool getDefaulted(string name) const {return false;}

    // This is a routine useful for debugging: return the name of the atomic
    // map that a certain parameter in the vector belongs to.
    string atomHavingParameter(int parameterIndex) const;
    // And a convenience to get the indices of a map's params without
    // having to issue the map.  Returns 0's if map is fixed or compound.
    void parameterIndicesOf(string mapname, int& startIndex, int& nParams) const;

    // This routine adds a new type of PhotoMap to the parsing dictionary
    template <class MapType>
    static void registerMapType() {
      mapCreators[MapType::type()] = MapType::create;
    }

    // This is a keyword that must appear in a YAML deserialization to
    // indicate a valid PixelMapCollection
    static const string magicKey;
    
  private:
    // hide:
    PhotoMapCollection(const PhotoMapCollection& rhs);
    void operator=(const PhotoMapCollection& rhs);

    int parameterCount; // Total parameters currently free to vary.
    int atomCount;	// Number of atomic map components
    int freeCount;	// Number of atompic components with free parameters

    // ***Main data of the class is in this container of structures listing
    // all the PhotoMaps curated by this class:

    // Structure for every PhotoMap that we know about:
    struct MapElement {
      MapElement(): realization(nullptr), atom(nullptr), isFixed(false), 
	nParams(0), number(-1), isValid(true) {}
      list<string> subordinateMaps;  // If it's compound, what it will be made from
      SubMap* realization;	     // pointer to its SubMap, if it's been built
      PhotoMap* atom;		     // Pointer to the PhotoMap itself, if atomic
      int startIndex;		     // Location in the union parameter array
      int nParams;		     // Number of parameters (only atomic is nonzero)
      int number;		     // sequential index among all maps with free parameters
      bool isFixed;		     // True if its parameters are currently fixed.
      bool isValid;		     // False if no info available to fit it
    };

    map<string, MapElement> mapElements;  // all known PhotoMaps, indexed by name.

    // Remove knowledge of a particular map or WCS and its realization if it exists.
    void removeMap(string mapName);

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
    PhotoMap* createAtomFromNode(const YAML::Node& node, string name) const;
    typedef PhotoMap* (*Creator)(const YAML::Node& node, string name);
    // Static map of what kind of PhotoMaps the parser will be able to identify
    static map<string,Creator> mapCreators;
    static bool creatorsInitialized;
    // Call to give parser an initial inventory of map types to create:
    static void PhotoMapTypeInitialize();

  };


} // namespace photometry
#endif  // PHOTOMAPCOLLECTION_H
