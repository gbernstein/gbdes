// A very restricted form of lookup table that is used by TemplateMap and PhotoTemplate
// in the astrometry and photometry packages, respectively.
//
// The lookup table argument is either x, y, or radius component of an (x,y) position.
//
// The value is a 1d correction (to either magnitude or the x/y/r coordinate).
//
// Lookup table nodes are assumed to be linearly spaced.
//
// Table information is read from YAML files. The root node of the YAML file is
// a map with a string name for each distinct table.  All tables in a file are
// read and stored in a static internal cache using the ingestFile() method.
// There can be multiple caches.  Then one requests pointers to read-only tables
// by their name (and cache source) with static find() method.  Constructor
// is private.  No provision for clearing caches at this point.
//
// Duplicate table names within a cache will throw exception.
//
// Methods of the class, once created, are documented below.

#ifndef LOOKUP1D_H
#define LOOKUP1D_H

#include <string>
#include <list>
#include <set>
#include <map>
#include <iostream>
#include <vector>
#include "Std.h"

#ifdef USE_YAML  // We won't be able to construct any instance of this without YAML
#include "yaml-cpp/yaml.h"
#endif

namespace lookup {
  class LookupError: public std::runtime_error {
  public:
    LookupError(const std::string &m=""): std::runtime_error("Lookup1d Error: " +m) {}
  };

  class Lookup1d {
  public:
    enum Coordinate {X,Y,R}; // Which axis is the Table being applied to
    Coordinate getAxis() const {return axis;}
    double getArg(double x, double y) const; // Return appropriate x/y/r value
    void getCenter(double& cx, double& cy) const {cx=xCenter; cy=yCenter;}

    // Read value (or value and derivative) of the table
    void valueSlope(double arg, double& val, double& slope) const;
    double value(double arg) const {
      double val, s;
      valueSlope(arg,val,s);
      return val;
    }
    void valueSlope(double x, double y, double& val, double& slope) const {
      return valueSlope(getArg(x,y),val,slope);
    }
    double value(double x, double y) const {
      return value(getArg(x,y));
    }
    // Find location where arg + scale*f(arg) = val:
    double inverse(double val, double scale) const;

#ifdef USE_YAML
    // Serialization routine, really for documentation purposes
    // as we have not given ways to make tables here:
    void write(YAML::Emitter& os) const;
#endif
    
    // Cache management routines:  first, get all tables from some file.
    // Seperate named caches can be kept.
    static void ingestFile(const string& filename, const string& cachename="default");

    // If an environment variable matching this path is defined, it will
    // be used as colon-separated list of paths in which to search for table files.
    static string CAL_PATH;

    // Request pointer to named Lookup1d.  If it does not exist in cache, return 0
    static const Lookup1d* find(const string& tablename, const string& cachename="default");
    
    typedef std::map<string, const Lookup1d*> Cache;

  private:
#ifdef USE_YAML
    Lookup1d(const YAML::Node& node);
#endif
    // Hide copy constructor & assignment
    Lookup1d(const Lookup1d& rhs);
    void operator=(const Lookup1d& rhs);

    Coordinate axis;
    double argStart;
    double argStep;
    vector<double> v;
    vector<double> dvda;
    // For radial functions, define center:
    double xCenter;
    double yCenter;
    

    // The caches
    static std::map<string, Cache> caches;
    // Names of files already ingested
    static std::map<string,std::set<string> > cachedFilenames;
  };

} // end namespace lookup

#endif  // End include gaurd


