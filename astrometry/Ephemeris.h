// $Id: Ephemeris.h,v 1.2 2009/11/02 15:54:48 garyb Exp $
// Classes to interface to the JPL DE405 Planetary Ephemeris.
// This is built atop David Hoffman's C interface software to the DE405.
// See earlier orbit-code distribution.
#ifndef EPHEMERIS_H
#define EPHEMERIS_H

#include "Astrometry.h"
#include <string>
#include "JPL.h"	//C definitions from Hoffman code, only used here.
#include <cstdio>
#include <deque>

namespace ephem {
  using namespace astrometry;

  enum Planet {
    MERCURY=0,
    VENUS=1,
    EARTH=2,		/*actually E-M barycenter for JPL tables*/
    MARS=3,
    JUPITER=4,
    SATURN=5,
    URANUS=6,
    NEPTUNE=7,
    PLUTO=8,
    MOON=9,		/* Moon relative to geocenter */
    SUN=10,		/* All coords (including this) relative to SSBARY*/
    NUTATIONS=11,
    LIBRATIONS=12
  };

  class StateCoords {
  public:
    StateCoords(): x(), v() {}
    StateCoords(const Vector3& x_, const Vector3& v_): x(x_), v(v_) {}
    CartesianICRS x;
    CartesianICRS v;
    bool operator==(const StateCoords& rhs) const {
      return (x==rhs.x && v==rhs.v);
    }
    StateCoords& operator-=(const StateCoords& rhs) {
      x-=rhs.x; v-=rhs.v; return *this;
    }
    StateCoords& operator+=(const StateCoords& rhs) {
      x+=rhs.x; v+=rhs.v; return *this;
    }
  };

  // class represents the JPL ephemeris
  class JPL {
  public:
    JPL(string ephemFile="");
    ~JPL();
    CartesianICRS position(Planet target, const UT& t);
    StateCoords state(Planet target,  const UT& t);
    // This function makes (if needed) and returns pointer to 
    // a static realization, i.e. it always gives a system-wide version.
    static JPL* theEphemeris(string filename="");
  private:
    HeaderOneType  H1;
    HeaderTwoType  H2;
    RecOneType   R1;
    std::FILE    *Ephemeris_File;
    double       Coeff_Array[ARRAY_SIZE];
    double       T_beg , T_end , T_span;
    bool         Debug;
    void	 ReadCoefficients(double Time);
    // The interpolation function; v=0 if no velocities desired.
    void         interpolate(Planet target, const UT& t, 
			     Vector3* x, Vector3* v=0);
  };


  class Observatory {
  protected:
    int code;
    string name;
  public:
    Observatory(): code(-1), name("unnamed") {}
    Observatory(int code_, const string& name_): code(code_), name(name_) {}
    virtual ~Observatory() {}
    int getCode() const {return code;}
    string getName() const {return name;}
    virtual CartesianICRS position(const UT& t) const =0;
    virtual StateCoords   state(const UT& t)    const =0;

    static  Observatory* build(const string& s);
    virtual void write(std::ostream& os) const =0;
    friend  std::ostream& operator<<(std::ostream& os, const Observatory& rhs);
  };

  std::ostream& operator<<(std::ostream& os, const Observatory& rhs); 

  class GroundObservatory: public Observatory {
  private:
    double lat;		//Observatory position in radians
    double lon;
    double altitude;    // and altitude in m
    mutable JPL* ephemeris;  
  public:
    GroundObservatory(int code_, string name_, double lat_, double lon_, 
		      double alt_, JPL* eph=JPL::theEphemeris()):
      Observatory(code_, name_), lat(lat_), lon(lon_), altitude(alt_),
      ephemeris(eph) {}
    GroundObservatory(const string& spec);
    void useEphemeris(JPL* eph) {ephemeris=eph;}
    CartesianICRS position(const UT& t) const;
    StateCoords   state(const UT& t)    const;
    void write(std::ostream& os) const;
  };

  // "Observatory" at geocenter
  class GeocenterObservatory: public Observatory {
  private:
    mutable JPL* ephemeris;
  public:
    GeocenterObservatory(JPL* eph=JPL::theEphemeris());
    CartesianICRS position(const UT& t) const {
      return ephemeris->position(EARTH,t);
    }
    StateCoords   state(const UT& t)    const {
      return ephemeris->state(EARTH,t);
    }
    void write(std::ostream& os) const;
  };
  // "Observatory" at Solar System Barycenter - origin of the coord system.
  class BarycenterObservatory: public Observatory {
  public:
    BarycenterObservatory();
    CartesianICRS position(const UT& t) const {return CartesianICRS();}
    StateCoords   state(const UT& t)    const {return StateCoords();}
    void write(std::ostream& os) const;
  };

  // class OrbitingObservatory;
  // class HSTObservatory;

  // Class which holds the full catalog of observatories.
  // Destroying this also destroys all the observatories.
  class ObservatoryCatalog {
  private:
    vector<Observatory*> vobs;
    ObservatoryCatalog(const ObservatoryCatalog& rhs) {}	//no copying
    void operator=(const ObservatoryCatalog& rhs) {}
  public:
    ObservatoryCatalog(string filename="");
    ~ObservatoryCatalog();
    const Observatory* find(int obscode) const;
    // This function makes (if needed) and returns reference to 
    // a static realization, i.e. it always gives a system-wide version.
    static const ObservatoryCatalog& theCatalog(const string filename="");
  };

  

  // Class for integrated orbits of KBOs
  class NumericalOrbit {
  protected:
    StateCoords xv0;	//Initial state & time
    UT          t0;
    // Calculate time of light emission for reception by observer at obstime
    double timeDelay(const CartesianICRS& observer,
		     const UT& obstime) const;
  public:
    NumericalOrbit(const StateCoords& xv0_, const UT& t0_):
      xv0(xv0_), t0(t0_) {}
    virtual ~NumericalOrbit() {}
    virtual CartesianICRS position(const UT& t) const =0;
    virtual StateCoords   state(const UT& t) const =0;
    virtual CartesianICRS position(const UT& t, DMatrix& derivs) const=0;
    virtual StateCoords   state(const UT& t, DMatrix& derivs) const=0;
    virtual void flush() const {}	//flush saved states
    CartesianICRS astrometricPosition(const CartesianICRS& observer,
				      const UT& obstime) const;
    CartesianICRS astrometricPosition(const CartesianICRS& observer,
				      const UT& obstime,
				      DMatrix& derivs) const;
    // ?? Apparent motion as well??
  };

  // This one is a real integrator
  class NBodyOrbit: public NumericalOrbit {
  private:
    mutable     JPL* ephem;
    double	tStep;	//Time step to use for integrations
    // Keep all the computed time steps for this objects
    mutable     std::deque<StateCoords> statedeque;
    mutable     double tFirst;	//time of first/last list element.
    mutable     double tLast;
    CartesianICRS acceleration(CartesianICRS x,
				     const UT& t) const;
    void integrateTo(const UT& t) const;
    StateCoords interpolate(const UT& t) const;
    void makeDerivs(const UT& t, DMatrix& derivs,
		    StateCoords exact, bool doVelocity) const;
  public:
    NBodyOrbit(const StateCoords& xv0_, const UT& t0_,
	       double tStep_=10.*DAY,
	       JPL* eph=JPL::theEphemeris());
    CartesianICRS position(const UT& t) const;
    StateCoords   state(const UT& t) const;
    CartesianICRS position(const UT& t, DMatrix& derivs) const;
    StateCoords   state(const UT& t, DMatrix& derivs) const;
    void flush() const;
  };

  // No-gravity orbit "integrator"
  class InertialOrbit: public NumericalOrbit {
  private:
  void makeDerivs(const UT& t, DMatrix& derivs,
		  bool doVelocity) const;
  public:
    InertialOrbit(const StateCoords& xv0_, const UT& t0_):
      NumericalOrbit(xv0_,t0_) {}
    CartesianICRS position(const UT& t) const;
    StateCoords   state(const UT& t) const;
    CartesianICRS position(const UT& t, DMatrix& derivs) const;
    StateCoords   state(const UT& t, DMatrix& derivs) const;
  };


  class ABGState {
  private:
    double a, b, g, adot, bdot, gdot;
    ReferenceFrame frame;
    UT  UT0;
    mutable NumericalOrbit *itsOrbit;
    bool useInertialOrbit;
    void readit(std::istream& is, DMatrix* m);
    void dXYZdABG(DMatrix& d) const;
    void buildOrbit() const;
  public:
    ABGState();
    ABGState(double a_,double b_,double g_,
	     double adot_, double bdot_, double gdot_,
	     const UT& UT0_, const ReferenceFrame& frame_,
	     bool useInertial=false); 
    ABGState(const StateCoords& xv, const UT& epoch);
    ABGState(const ABGState& rhs);
    const ABGState& operator=(const ABGState& rhs);
    ~ABGState();
    double alpha() const {return a;}
    double beta() const {return b;}
    double gamma() const {return g;}
    double alphadot() const {return adot;}
    double betadot() const {return bdot;}
    double gammadot() const {return gdot;}
    double& alpha() {return a;}
    double& beta() {return b;}
    double& gamma() {return g;}
    double& alphadot() {return adot;}
    double& betadot() {return bdot;}
    double& gammadot() {return gdot;}
    DVector getVector() const;
    UT t0() const {return UT0;}
    const ReferenceFrame& refFrame() const {return frame;}
    StateCoords xv() const;
    StateCoords xv(DMatrix& derivs) const;
    CartesianICRS position(const UT& t) const {
      buildOrbit(); return itsOrbit->position(t);
    }
    StateCoords   state(const UT& t) const {
      buildOrbit(); return itsOrbit->state(t);
    }

    void setFrom(double a_,double b_,double g_,
		 double adot_, double bdot_, double gdot_,
		 const UT& UT0_, const ReferenceFrame& frame_,
		 bool useInertial=false); 
    void setFrom(const DVector& p,
		 const UT& UT0_, const ReferenceFrame& frame_); 
    void setFrom(const StateCoords& s,
		 const UT& UT0_, const ReferenceFrame& frame_); 
    void setFrom(const DVector& p);
    // read from old-style file, with covar matrix
    void readOldFormat(std::istream& is);
    void readOldFormat(std::istream& is, DMatrix& covar);
    // Note that neither of below write/read the reference frame info.
    void write(std::ostream& os) const;
    void read(std::istream& is);
    // Read/write ABG along with its reference frame and uncertainties
    void read(std::istream& is, 
	      DMatrix& covar);
    void write(std::ostream& os, 
	       DMatrix& covar) const;
    void setInertial(bool useInertial) {useInertialOrbit=useInertial;}

    // Get ICRS position at time, deriv matrix w.r.t. the 6 params
    SphericalICRS observation(const CartesianICRS& observer,
			      const UT time) const;
    SphericalICRS observation(const CartesianICRS& observer,
			      const UT time,
			      double& range) const;
    SphericalICRS observation(const CartesianICRS& observer,
			      const UT time,
			      DMatrix& derivs) const;
  };

  class OrbitalElements {
  private:
    double a, e, i;
    double longitudeOfAscendingNode;
    double argumentOfPerifocus;
    double meanAnomaly;
    UT     epoch;
    void fromXV(const StateCoords& xv, const UT& UT0);
    void toXV(StateCoords& xv, const UT& UT0) const;
    bool heliocentric;
  public:
    OrbitalElements() {}
    explicit OrbitalElements(const ABGState& agb, bool helio=false);
    explicit OrbitalElements(const StateCoords& xv, const UT& UT0,
			     bool helio=false);
    StateCoords getXV(const UT& UT0) const;
    UT getEpoch() const {return epoch;}
    void write(std::ostream& os) const;
    void writeHeader(std::ostream& os) const;
    bool read(std::istream& is);
    double getA() const {return a;}
    double getE() const {return e;}
    double getI() const {return i;}
    double getLAN() const {return longitudeOfAscendingNode;}
    double getAOP() const {return argumentOfPerifocus;}
    double getM() const {return meanAnomaly;}
    bool read(std::istream& is, bool helio=false);
  };

  void
  aei_derivs(const ephem::StateCoords& xv,
	     DMatrix& daei_dxv);

} // namespace ephem
#endif

