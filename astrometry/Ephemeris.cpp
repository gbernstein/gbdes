/* 	$Id: Ephemeris.cpp,v 1.4 2011/02/14 20:33:27 garyb Exp $	 */
#ifndef lint
static char vcid[] = "$Id: Ephemeris.cpp,v 1.4 2011/02/14 20:33:27 garyb Exp $";
#endif /* lint */

#include "Ephemeris.h"
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <fstream>
#include <sstream>
#include "StringStuff.h"

using namespace std;	//for all the JPL code

namespace ephem {

  using namespace astrometry;

  JPL* JPL::theEphemeris(string filename) {
    static bool initialized=false;
    static JPL* globalRealization;
    if (!initialized) {
      globalRealization = new JPL(filename);
      initialized = true;
    }
    return globalRealization;
  }

  // Default file & name of environment variable which can hold filename
  // for the binary JPL ephemeris
  const string EPHEM_ENVIRON="ORBIT_EPHEMERIS";
  const string DEFAULT_EPHEM_FILE="binEphem.405";
  const string OBSERVATORY_ENVIRON="ORBIT_OBSERVATORIES";
  const string DEFAULT_OBSERVATORY_FILE="observatories.dat";
  const int    OBSCODE_ORBITAL=2000;	//obscode above here are orbiting Earth
  const int    OBSCODE_GEOCENTER=500;	//obscode for geocenter
  const int    OBSCODE_SSBARY=1999;	//obscode for solar system barycenter
  const int    OBSCODE_HST=250;	        //obscode for HST

  // JPL-reading routines are from Hoffman's ephem_read.c

  //    This function is used by the functions below to read an array of     **/
  //    Tchebeychev coefficients from a binary ephemeris data file.          **/

  void 
  JPL::ReadCoefficients( double Time ) {
    double  T_delta = 0.0;
    long     Offset  =  0 ;		/*** ??? change to long 8/9/99 ***/

    /*--------------------------------------------------------------------------*/
    /*  Find ephemeris data that record contains input time. Note that one, and */
    /*  only one, of the following conditional statements will be true (if both */
    /*  were false, this function would not have been called).                  */
    /*--------------------------------------------------------------------------*/
    if ( Time < T_beg )                    /* Compute backwards location offset */
      {
	T_delta = T_beg - Time;
	Offset  = (int) -ceil(T_delta/T_span);	/***Needed negative sign here???*/
      }

    if ( Time > T_end )                    /* Compute forewards location offset */
      {
	T_delta = Time - T_end;
	Offset  = (int) ceil(T_delta/T_span);
      }

    /*--------------------------------------------------------------------------*/
    /*  Retrieve ephemeris data from new record.                                */
    /*--------------------------------------------------------------------------*/

    fseek(Ephemeris_File,(Offset-1)*ARRAY_SIZE*sizeof(double),SEEK_CUR);
    fread(&Coeff_Array,sizeof(double),ARRAY_SIZE,Ephemeris_File);
  
    T_beg  = Coeff_Array[0];
    T_end  = Coeff_Array[1];
    T_span = T_end - T_beg;

    if (Time < T_beg || Time > T_end) {
      fprintf(stderr,"JD %.2f is out of range of ephemeris file\n",Time);
      exit(1);
    }
    /*--------------------------------------------------------------------------*/
    /*  Debug print (optional)                                                  */
    /*--------------------------------------------------------------------------*/
    if ( Debug ) 
      {
	printf("\n  In: ReadCoefficients \n");
	printf("\n      ARRAY_SIZE = %4d",ARRAY_SIZE);
	printf("\n      Offset  = %3d",Offset);
	printf("\n      T_delta = %7.3f",T_delta);
	printf("\n      T_Beg   = %7.3f",T_beg);
	printf("\n      T_End   = %7.3f",T_end);
	printf("\n      T_Span  = %7.3f\n\n",T_span);
      }
  }

  // Constructor: opens file, reads header data & loads first record.
  JPL::JPL(string filename): Debug(false) {
    int headerID;
    /*--------------------------------------------------------------------------*/
    /*  Open ephemeris file.                                                    */
    /*--------------------------------------------------------------------------*/

    // use requested filename, or environment-specified file, or the default filename,
    // in that order
    if (filename.empty()) {
      if (getenv(EPHEM_ENVIRON.c_str())!=NULL)
	filename=getenv(EPHEM_ENVIRON.c_str());
      else
	filename=DEFAULT_EPHEM_FILE;
    }
  
    Ephemeris_File = fopen(filename.c_str(),"r");
    
    /*--------------------------------------------------------------------------*/
    /*  Read header & first coefficient array, then return status code.         */
    /*--------------------------------------------------------------------------*/
    
    if ( Ephemeris_File == NULL ) /*........................No need to continue */
      throw AstrometryError("Cannot open JPL ephemeris file <"
			    + filename + ">");
    
    /*.................Read first three header records from ephemeris file */
         
    fread(&H1,sizeof(double),ARRAY_SIZE,Ephemeris_File);
    fread(&H2,sizeof(double),ARRAY_SIZE,Ephemeris_File);
    fread(&Coeff_Array,sizeof(double),ARRAY_SIZE,Ephemeris_File);
    
    /*...............................Store header data in global variables */
    
    R1 = H1.data;
    
    /*..........................................Set current time variables */
    
    T_beg  = Coeff_Array[0];
    T_end  = Coeff_Array[1];
    T_span = T_end - T_beg;
    
    /*..............................Convert header ephemeris ID to integer */
    
    headerID = (int) R1.DENUM;
    
    /*..............................................Debug Print (optional) */
    
    if ( Debug ) 
      {
	printf("\n  In: Initialize_Ephemeris \n");
	printf("\n      ARRAY_SIZE = %4d",ARRAY_SIZE);
	printf("\n      headerID   = %3d",headerID);
	printf("\n      T_Beg      = %7.3f",T_beg);
	printf("\n      T_End      = %7.3f",T_end);
	printf("\n      T_Span     = %7.3f\n\n",T_span);
      }

    /*..................................................Return status code */
       
    if ( headerID != EPHEMERIS ) 
      throw AstrometryError("Opened wrong file <" 
			    + filename + "> for JPL ephemeris");
  }

  JPL::~JPL() {
    fclose(Ephemeris_File);
  }


  CartesianICRS
  JPL::position( Planet target, const UT& t) {
    Vector3 x;
    interpolate(target, t, &x);

    // Earth & Moon are special, need EMBary posn and E->M vector:
    if (target==EARTH) {
      Vector3 xEM;
      interpolate(MOON, t, &xEM);
      xEM /= 1. + R1.EMRAT;
      x -= xEM;
    }
    if (target==MOON) {
      Vector3 xEMBary;
      interpolate(EARTH, t, &xEMBary);
      x *= R1.EMRAT/(1. + R1.EMRAT);
      x += xEMBary;
    }
    x *= 1./R1.AU;
    return CartesianICRS(x);
  }

  StateCoords
  JPL::state( Planet target, const UT& t) {
    Vector3 x, v;
    interpolate(target, t, &x, &v);

    // Earth & Moon are special, need EMBary posn and E->M vector:
    if (target==EARTH) {
      Vector3 xEM, vEM;
      interpolate(MOON, t, &xEM, &vEM);
      xEM /= 1. + R1.EMRAT;
      x -= xEM;
      vEM /= 1. + R1.EMRAT;
      v -= vEM;
    }
    if (target==MOON) {
      Vector3 xEMBary, vEMBary;
      interpolate(EARTH, t, &xEMBary, &vEMBary);
      x *= R1.EMRAT/(1. + R1.EMRAT);
      x += xEMBary;
      v *= R1.EMRAT/(1. + R1.EMRAT);
      v += vEMBary;
    }
    x *= 1./R1.AU;
    v *= 1./(SECOND*R1.AU);
    return StateCoords(x,v);
  }

  // interpolate from the JPL table.  v=0 if velocity isn't required.
  void
  JPL::interpolate( Planet target, const UT& t,
		    Vector3* x, Vector3* v) {

    double    A[50]   , B[50] , Cp[50] , P_Sum[3] , V_Sum[3] , Up[50] ,
      T_break , T_seg , T_sub  , Tc;
    int       i , j;
    long int  C , G , N , offset = 0;

    double Time=t.getTT();	//TT is within 0.002s of TDB, the ephemeris' time

    if ( target >= NUTATIONS )
      throw AstrometryError("JPL code cannot handle request for nutations or librations");

    /*--------------------------------------------------------------------------*/
    /* Initialize local coefficient array.                                      */
    /*--------------------------------------------------------------------------*/

    for ( i=0 ; i<50 ; i++ )
      {
        A[i] = 0.0;
        B[i] = 0.0;
      }

    /*--------------------------------------------------------------------------*/
    /* Determine if a new record needs to be input.                             */
    /*--------------------------------------------------------------------------*/
  
    if (Time < T_beg || Time > T_end)  ReadCoefficients(Time);

    /*--------------------------------------------------------------------------*/
    /* Read the coefficients from the binary record.                            */
    /*--------------------------------------------------------------------------*/
  
    C = R1.coeffPtr[target][0] - 1;               /*    Coeff array entry point */
    N = R1.coeffPtr[target][1];                   /*          Number of coeff's */
    G = R1.coeffPtr[target][2];                   /* Granules in current record */

    /*...................................................Debug print (optional) */

    if ( Debug )
      {
	printf("\n  In: Interpolate_State\n");
	printf("\n  target = %2d",target);
	printf("\n  C      = %4d (before)",C);
	printf("\n  N      = %4d",N);
	printf("\n  G      = %4d\n",G);
      }

    /*--------------------------------------------------------------------------*/
    /*  Compute the normalized time, then load the Tchebeyshev coefficients     */
    /*  into array A[]. If T_span is covered by a single granule this is easy.  */
    /*  If not, the granule that contains the interpolation time is found, and  */
    /*  an offset from the array entry point for the ephemeris body is used to  */
    /*  load the coefficients.                                                  */
    /*--------------------------------------------------------------------------*/

    if ( G == 1 ) {
	Tc = 2.0*(Time - T_beg) / T_span - 1.0;
	for (i=C ; i<(C+3*N) ; i++)  A[i-C] = Coeff_Array[i];
    } else if ( G > 1 ) {
	T_sub = T_span / ((double) G);          /* Compute subgranule interval */
       
	for ( j=G ; j>0 ; j-- ) {
	  T_break = T_beg + ((double) j-1) * T_sub;
	  if ( Time > T_break ) {
	    T_seg  = T_break;
	    offset = j-1;
	    break;
	  }
	}
            
	Tc = 2.0*(Time - T_seg) / T_sub - 1.0;
	C  = C + 3 * offset * N;
       
	for (i=C ; i<(C+3*N) ; i++) A[i-C] = Coeff_Array[i];
    } else              /* Something has gone terribly wrong */
      throw AstrometryError(" JPL granules must be >= 1: check header data.");

    /*...................................................Debug print (optional) */

    if ( Debug ) {
      printf("\n  C      = %4d (after)",C);
      printf("\n  offset = %4d",offset);
      printf("\n  Time   = %12.7f",Time);
      printf("\n  T_sub  = %12.7f",T_sub);
      printf("\n  T_seg  = %12.7f",T_seg);
      printf("\n  Tc     = %12.7f\n",Tc);
      printf("\n  Array Coefficients:\n");
      for ( i=0 ; i<3*N ; i++ ) {
	printf("\n  A[%2d] = % 22.15e",i,A[i]);
      }
      printf("\n\n");
    }

    /*..........................................................................*/

    /*--------------------------------------------------------------------------*/
    /* Compute the interpolated position & velocity                             */
    /*--------------------------------------------------------------------------*/
  
    for ( i=0 ; i<3 ; i++ )  {              /* Compute interpolating polynomials */
      Cp[0] = 1.0;           
      Cp[1] = Tc;
      Cp[2] = 2.0 * Tc*Tc - 1.0;
      
      Up[0] = 0.0;
      Up[1] = 1.0;
      Up[2] = 4.0 * Tc;

      for ( j=3 ; j<N ; j++ )
	Cp[j] = 2.0 * Tc * Cp[j-1] - Cp[j-2];
      if (v) for ( j=3 ; j<N ; j++ ) 
	Up[j] = 2.0 * Tc * Up[j-1] + 2.0 * Cp[j-1] - Up[j-2];

      P_Sum[i] = 0.0;           /* Compute interpolated position & velocity */
      V_Sum[i] = 0.0;

      for ( j=N-1 ; j>-1 ; j-- )  P_Sum[i] = P_Sum[i] + A[j+i*N] * Cp[j];
      if (v) for ( j=N-1 ; j>0  ; j-- )  V_Sum[i] = V_Sum[i] + A[j+i*N] * Up[j];

      (*x)[i] = P_Sum[i];
      if (v) (*v)[i] = V_Sum[i] * 2.0 * ((double) G) / (T_span * 86400.0);
    }

    return;
  }

  std::ostream& operator<<(std::ostream& os, const Observatory& rhs) {
    rhs.write(os); return os;
  }

  // From Skycalc:
  double 
  lst(double jd,  double longit)
  {
    /* returns the local MEAN sidereal time (radians) at julian date jd
       at west longitude long (radians).  Follows
       definitions in 1992 Astronomical Almanac, pp. B7 and L2.
       Expression for GMST at 0h ut referenced to Aoki et al, A&A 105,
       p.359, 1982.  On workstations, accuracy (numerical only!)
       is about a millisecond in the 1990s. */

    const double  J2000= 2451545.;        /* Julian date at standard epoch */
    const double   SEC_IN_DAY=86400.;
  
    double t, ut, jdmid, jdint, jdfrac, sid_g, sid;
    long jdin, sid_int;

    jdin = jd;         /* fossil code from earlier package which
			  split jd into integer and fractional parts ... */
    jdint = jdin;
    jdfrac = jd - jdint;
    if(jdfrac < 0.5) {
      jdmid = jdint - 0.5;
      ut = jdfrac + 0.5;
    } else {
      jdmid = jdint + 0.5;
      ut = jdfrac - 0.5;
    }
    t = (jdmid - J2000)/36525;
    sid_g = (24110.54841+8640184.812866*t+0.093104*t*t-6.2e-6*t*t*t)/SEC_IN_DAY;
    sid_int = sid_g;
    sid_g = sid_g - (double) sid_int;
    sid_g = sid_g + 1.0027379093 * ut - longit/TPI;
    sid_int = sid_g;
    sid_g = (sid_g - (double) sid_int) * TPI;
    if(sid_g < 0.) sid_g = sid_g + TPI;
    return(sid_g);
  }

  // From Skycalc:  give geocentric cartesian coordinates given
  // geodetic lat/lon (radians) and altitude above the reference
  // ellipsoid (meters).  Output in AU now.
  // See '92 Ast Almanac, p K11, or 4.22 of the supplement.

  // note that this gives proper equatorial vector if lon
  // is LMST measured positive eastward.

  CartesianICRS
  topo(double geolong, double geolat, double height) {
    const double  FLATTEN=0.003352813;   /* flattening of earth, 1/298.257 */
    const double  EQUAT_RAD=6378137.;    /* equatorial radius of earth, meters */

    double denom, C_geo, S_geo;

    denom = (1. - FLATTEN) * sin(geolat);
    denom = cos(geolat) * cos(geolat) + denom*denom;
    C_geo = 1. / sqrt(denom);
    S_geo = (1. - FLATTEN) * (1. - FLATTEN) * C_geo;
    C_geo = C_geo* EQUAT_RAD + height ;  /* deviation from almanac
					    notation -- include height here. */
    S_geo = S_geo* EQUAT_RAD + height ;

    CartesianICRS x;
    x[0] = C_geo * cos(geolat) * cos(geolong);
    x[1] = C_geo * cos(geolat) * sin(geolong);
    x[2] = S_geo * sin(geolat);

    /* convert to AU, keeping in mind that Horizons was km and this is m*/
    x *= 1./AU;
    return x;
  }

  // Give ICRS cartesian coordinates of ground-based observatory 
  // w.r.t. SS Barycenter.
  // ??? Note that at the present time there is no allowance for the
  // precession or nutation of Earth's pole.  Therefore positions
  // are incorrect by ~0.6 km /yr from 2000.
  CartesianICRS
  GroundObservatory::position(const UT& t) const {

    /* Get the LMST and calculate to ICRS vector */
    double obslmst=lst(t.getJD(), lon);
    CartesianICRS geoToObs=topo(obslmst, lat, altitude);

    /* Use Horizons to get the coordinates of geocenter */
    CartesianICRS geocenter=ephemeris->position(EARTH,t);

    geoToObs += geocenter;
    return geoToObs;
  }

  // ???
  StateCoords
  GroundObservatory::state(const UT& t) const {
    throw AstrometryError("GroundObservatory::state() not implemented");
  }

  GroundObservatory::GroundObservatory(const string& spec) {
    std::istringstream iss(spec);
    if (!(iss >> code)) 
      throw AstrometryError("Bad obscode for GroundObservatory");
    string buffer;
    // ??? Update this using AngularCoords
    if (!(iss >> buffer)) throw AstrometryError("Bad longitude for GroundObservatory: " + spec);
    lon = dmsdeg(buffer) * DEGREE;
    if (!(iss >> buffer)) throw AstrometryError("Bad latitude for GroundObservatory");
    lat = dmsdeg(buffer) * DEGREE;

    if (!(iss >> altitude)) throw AstrometryError("Bad altitude for GroundObservatory");
    name.clear();
    while (iss >> buffer)
      name += buffer + " ";
    ephemeris = JPL::theEphemeris();
  }

  void 
  GroundObservatory::write(std::ostream& os) const {
    // ??? use angularCoords:
    os << setw(3) << setfill('0') << code << setfill(' ')
       << " " << degdms(lon/DEGREE,1)
       << " " << degdms(lat/DEGREE,1)
       << " " << setprecision(1) << fixed << setw(7) << altitude
       << " " << name;
  }

  GeocenterObservatory::GeocenterObservatory(JPL* eph):
    Observatory(OBSCODE_GEOCENTER, "Geocenter"), ephemeris(eph) {}
  void 
  GeocenterObservatory::write(std::ostream& os) const {
      os << setw(3) << setfill('0') << code << setfill(' ') 
	 << " " << name;
  }
  BarycenterObservatory::BarycenterObservatory(): 
    Observatory(OBSCODE_SSBARY, "Solar System Barycenter") {}
  void 
  BarycenterObservatory::write(std::ostream& os) const {
      os << setw(3) << setfill('0') << code << setfill(' ') 
	 << " " << name;
  }

  Observatory*
  Observatory::build(const string& s) {
    int obscode;
    std::istringstream iss(s);
    if (!(iss >> obscode))
      throw AstrometryError("Format error in observatory specification");
    if (obscode==OBSCODE_GEOCENTER)
      return new GeocenterObservatory();
    else if (obscode==OBSCODE_SSBARY)
      return new BarycenterObservatory();
    else if (obscode > 0 && obscode<OBSCODE_ORBITAL)
      return new GroundObservatory(s);
    // ??? add new types of observatories here
    else 
      throw AstrometryError("Cannot interpret observatory spec:\n"+s);
  }


  const ObservatoryCatalog& 
  ObservatoryCatalog::theCatalog(string filename) {
    static bool initialized=false;
    static ObservatoryCatalog* globalRealization;
    if (!initialized) {
      globalRealization = new ObservatoryCatalog(filename);
      initialized = true;
    }
    return *globalRealization;
  }

  ObservatoryCatalog::ObservatoryCatalog(string filename) {
    // use requested filename, or environment-specified file, or the default filename,
    // in that order
    if (filename.empty()) {
      if (getenv(OBSERVATORY_ENVIRON.c_str())!=NULL)
	filename=getenv(OBSERVATORY_ENVIRON.c_str());
      else
	filename=DEFAULT_OBSERVATORY_FILE;
    }

    std::ifstream ifs(filename.c_str());
    if (!ifs) 
      throw AstrometryError("Could not open observatories file <"
			    + filename + ">");
    string linebuffer;
    while (getlineNoComment(ifs, linebuffer)) {
      vobs.push_back(Observatory::build(linebuffer));
    }
  }
    
  ObservatoryCatalog::~ObservatoryCatalog() {
    for (vector<Observatory*>::const_iterator i=vobs.begin();
	 i!=vobs.end();
	 ++i)
      delete *i;
  }
  const Observatory*
  ObservatoryCatalog::find(int obscode) const {
    static int lastcode=-999;
    static const Observatory* lastobs;
    if (obscode==lastcode) return lastobs;
    for (vector<Observatory*>::const_iterator i=vobs.begin();
	 i!=vobs.end();
	 ++i)
      if ((*i)->getCode()==obscode) {
	lastcode = obscode;
	lastobs=*i;
	return *i;
      }
    return 0;
  }

  /////////////////////////////////////////////////////////////
  // Orbit basis manipulation
  /////////////////////////////////////////////////////////////

  ABGState::ABGState(): a(0.),b(0.),g(0.02),
			adot(0.),bdot(0.),gdot(0.),
			frame(CartesianICRS(), Orientation()),
			UT0(2000,1,1.), itsOrbit(0),
                        useInertialOrbit(false) {}

  ABGState::ABGState(double a_,double b_,double g_,
		     double adot_, double bdot_, double gdot_,
		     const UT& UT0_, const ReferenceFrame& frame_,
		     bool useInertial):
    itsOrbit(0)  {
    setFrom(a_,b_,g_,adot_,bdot_,gdot_,UT0_,frame_, useInertial);
  }

  ABGState::ABGState(const StateCoords& xv, const UT& epoch): itsOrbit(0) {
    // Make the frame be barycentric, directed to position
    CartesianICRS origin;
    SphericalICRS pole(xv.x);
    ReferenceFrame ff(origin, pole);
    setFrom(xv, epoch, ff);
  }

  ABGState::~ABGState() {if (itsOrbit) delete itsOrbit;}

  ABGState::ABGState(const ABGState& rhs): itsOrbit(0) {
    setFrom(rhs.a, rhs.b, rhs.g, rhs.adot, rhs.bdot, rhs.gdot,
	    rhs.UT0, rhs.frame, rhs.useInertialOrbit);
  }

  const ABGState& 
  ABGState::operator=(const ABGState& rhs) {
    if (this==&rhs) return *this;
    setFrom(rhs.a, rhs.b, rhs.g, rhs.adot, rhs.bdot, rhs.gdot,
	    rhs.UT0, rhs.frame, rhs.useInertialOrbit);
    return *this;
  }

  void
  ABGState::setFrom(double a_,double b_,double g_,
		    double adot_, double bdot_, double gdot_,
		    const UT& UT0_, const ReferenceFrame& frame_,
		    bool useInertial) {
    if (itsOrbit) {
      delete itsOrbit;
      itsOrbit = 0;
    }
    a = a_; b=b_; g=g_; adot=adot_; bdot=bdot_; gdot=gdot_;
    UT0=UT0_; 
    frame=frame_;
    useInertialOrbit = useInertial;
  }
  void
  ABGState::setFrom(const DVector& p,
		    const UT& UT0_, const ReferenceFrame& frame_) {
    Assert(p.size()>=6);
    setFrom(p[0], p[1], p[2], p[3], p[4], p[5], UT0_, frame_,
	    useInertialOrbit);
  }
  void
  ABGState::setFrom(const DVector& p) {
    setFrom(p, UT0, frame);
  }

  DVector
  ABGState::getVector() const {
    DVector p(6);
    p[0] = a; p[1]=b; p[2]=g; p[3]=adot; p[4]=bdot; p[5]=gdot;
    return p;
  }

  void
  ABGState::setFrom(const StateCoords& s,
		    const UT& UT0_, const ReferenceFrame& frame_) {
    // Rotate & translate posn into ref. frame
    CartesianCustom xc(s.x ,frame_);
    // Rotate velocity into ref frame (no translation!)
    Vector3 v=frame_.orient.fromICRS(s.v.getVector());
    // Convert to ABG:
    v  *= 1./xc[2];
    setFrom(xc[0]/xc[2], xc[1]/xc[2], 1./xc[2], 
	    v[0], v[1], v[2], UT0_, frame_, useInertialOrbit);
  }

  void
  ABGState::readOldFormat(std::istream& is) {
    ABGState::readit(is, 0);
  }

  void
  ABGState::readOldFormat(std::istream& is, DMatrix& covar) {
    ABGState::readit(is, &covar);
  }

  void
  ABGState::readit(std::istream& is, DMatrix *m) {
    string buffer;
    // First read abg:
    if (!getlineNoComment(is, buffer)) return;
    {
      std::istringstream iss(buffer);
      iss >> a >> adot >> b >> bdot >> g >> gdot;
      if (!iss) {
	is.setstate(is.rdstate() | std::ios_base::failbit);
	return;
      }
    }
    // And the covariance matrix:
    if (m) Assert( m->nrows()==6 && m->ncols()==6);
    for (int i=0; i<6; i++) {
      if (!getlineNoComment(is, buffer)) return;
      if (!m) continue;
      std::istringstream iss(buffer);
      iss >> (*m)(i,0) >> (*m)(i,1) >> (*m)(i,2) 
	  >> (*m)(i,3) >> (*m)(i,4) >> (*m)(i,5);
      if (!iss) {
	is.setstate(is.rdstate() | std::ios_base::failbit);
	return;
      }
    }
    // And specify the reference frame & time
    double lat, lon, jd0;
    Vector3 x;
 
    if (!getlineNoComment(is, buffer)) return;
    {
      std::istringstream iss(buffer);
      iss >> lat >> lon >> x[0] >> x[1] >> x[2] >> jd0;
      if (!iss) {
	is.setstate(is.rdstate() | std::ios_base::failbit);
	return;
      }
    }

    x *= -1.;	//Convention in old files is earth->SSBary
    SphericalEcliptic refDirection(lon*DEGREE, lat*DEGREE);
    Orientation o(refDirection);
    o.alignToEcliptic();
    CartesianICRS refPt(o.toICRS(x));

    frame = ReferenceFrame(refPt,o);
    UT0.set(jd0);

    if (itsOrbit) {
      delete itsOrbit;
      itsOrbit = 0;
    }
  }

  void
  ABGState::write(std::ostream& os) const {
    os << std::setprecision(6) << std::fixed 
       << 1./gamma()
       << std::setprecision(3) << std::fixed
       << " " << alpha()/ARCSEC
       << " " << beta()/ARCSEC
       << " " << std::setprecision(5) << alphadot()*HOUR/ARCSEC
       << " " << betadot()*HOUR/ARCSEC
       << " " << gammadot()*HOUR/ARCSEC 
       << "   " << UT0 ;
  };

  void
  ABGState::read(std::istream& is) {
    double a,b,d,adot,bdot,gdot;
    is >> d >> a >> b >> adot >> bdot >> gdot >> UT0;
    setFrom(a*ARCSEC, b*ARCSEC, 1./d, 
	    adot*ARCSEC/HOUR, bdot*ARCSEC/HOUR, gdot*ARCSEC/HOUR,
	    UT0, frame, useInertialOrbit);
  };


  StateCoords
  ABGState::xv() const {
    StateCoords sc;
    {
      Vector3 x;
      x[0] = a/g; x[1] = b/g; x[2] = 1./g;
      CartesianCustom cc(x,frame);
      CartesianICRS xICRS(cc);
      sc.x = xICRS;
    }

    // Just apply rotation to velocities
    Vector3 v;
    v[0] = adot/g; v[1] = bdot/g; v[2] = gdot/g;
    sc.v = frame.orient.toICRS(v);
    return sc;
  }

  StateCoords
  ABGState::xv(DMatrix& derivs) const {
    dXYZdABG(derivs);
    return xv();
  }

  void
  ABGState::dXYZdABG(DMatrix& d) const {
    Assert(d.ncols()==6 && d.nrows()==6);
    d.setZero();
    // Here is the orientation's rotation matrix, which goes
    // from ICRS to our orient, so we'll need the inverse.
    // This could be done with matrix multiplication but
    // there are many zero elements, so I'll write it out:
    const Matrix33& m=frame.orient.m();
    double z = 1./g;
    /*
    d(0,0) = d(1,1) = d(3,3) = d(4,4) = d(5,5) = z;
    d(0,2) = -a*z*z;
    d(1,2) = -b*z*z;
    d(2,2) = -z*z;
    d(3,2) = -adot*z*z;
    d(4,2) = -bdot*z*z;
    d(5,2) = -gdot*z*z;
    SqDMatrix bigR(6,0.);
    for (int i=0; i<3; i++)
      for (int j=0; j<3; j++)
	bigR(i,j) = bigR(i+3,j+3)=m(j,i);
    d = bigR*d;
    return;
    */
    d(0,0) = z*m(0,0);
    d(1,0) = z*m(0,1);
    d(2,0) = z*m(0,2);
    d(0,1) = z*m(1,0);
    d(1,1) = z*m(1,1);
    d(2,1) = z*m(1,2);
    d(0,2) = -z*z*(a*m(0,0) + b*m(1,0) + m(2,0));
    d(1,2) = -z*z*(a*m(0,1) + b*m(1,1) + m(2,1));
    d(2,2) = -z*z*(a*m(0,2) + b*m(1,2) + m(2,2));
    d(3,2) = -z*z*(adot*m(0,0) + bdot*m(1,0) + gdot*m(2,0));
    d(4,2) = -z*z*(adot*m(0,1) + bdot*m(1,1) + gdot*m(2,1));
    d(5,2) = -z*z*(adot*m(0,2) + bdot*m(1,2) + gdot*m(2,2));
    d(3,3) = z*m(0,0);
    d(3,4) = z*m(1,0);
    d(3,5) = z*m(2,0);
    d(4,3) = z*m(0,1);
    d(4,4) = z*m(1,1);
    d(4,5) = z*m(2,1);
    d(5,3) = z*m(0,2);
    d(5,4) = z*m(1,2);
    d(5,5) = z*m(2,2);
  }

  void
  ABGState::buildOrbit() const {
    if (!itsOrbit) {
      if (useInertialOrbit)
	itsOrbit = new InertialOrbit(xv(), UT0);
      else
	itsOrbit = new NBodyOrbit(xv(), UT0);
    }
  }

  SphericalICRS
  ABGState::observation(const CartesianICRS& observer,
			const UT time) const {
    buildOrbit();
    CartesianICRS target=itsOrbit->astrometricPosition(observer, time);
    return SphericalICRS(target);
  }

  SphericalICRS
  ABGState::observation(const CartesianICRS& observer,
			const UT time,
			double& range) const {
    buildOrbit();
    CartesianICRS target=itsOrbit->astrometricPosition(observer, time);
    range = target.radius();
    return SphericalICRS(target);
  }

  SphericalICRS
  ABGState::observation(const CartesianICRS& observer,
			const UT time,
			DMatrix& derivs) const {
    buildOrbit();
    DMatrix dICRSdXYZ(3,6);

    CartesianICRS target=
      itsOrbit->astrometricPosition(observer, time, dICRSdXYZ);
    Matrix23 dLLdICRS;
    SphericalICRS posn(target,dLLdICRS);
    DMatrix dXdA(6,6);
    dXYZdABG(dXdA);
    derivs = dLLdICRS * dICRSdXYZ * dXdA;
    return posn;
  }
    
  // Write the ABG and its reference frame & uncertainties to ostreadm
  void 
  ABGState::write(std::ostream& os, 
		  DMatrix& covar) const {

    os << "# Reference frame: " << endl;
    os << frame << endl;
    os << "# Best fit: " << endl;
    write(os);
    os << endl;
    os << "#Standard deviations: " << endl;
    os << " " << std::setprecision(5) << std::setw(8)
       << sqrt(covar(2,2))/(gamma()*gamma())
       << std::setprecision(4) << std::fixed 
       << "     " << sqrt(covar(0,0))/ARCSEC
       << " " << sqrt(covar(1,1))/ARCSEC
       << std::setprecision(5)  
       << " " << sqrt(covar(3,3))*HOUR/ARCSEC
       << " " << sqrt(covar(4,4))*HOUR/ARCSEC
       << std::setprecision(5)  
       << " " << sqrt(covar(5,5))*HOUR/ARCSEC
       << endl;
    os << "# Correlation matrix: " << std::setprecision(3) << endl;
    for (int i=0; i<6; i++) {
      for (int j=0; j<6; j++) 
	cout << std::setw(6) << covar(i,j)/sqrt(covar(i,i)*covar(j,j)) << " ";
      cout << endl;
    }
  }
    
  // Read a full specification of ABG orbit & uncertainties from file
  void
  ABGState::read(std::istream& is, 
		 DMatrix& covar) {
    Assert( covar.nrows()==6 && covar.ncols()==6);
    // First two non-comment lines give the reference frame
    string buffer;
    getlineNoComment(is, buffer);
    string buffer2;
    getlineNoComment(is, buffer2);
    buffer += " " + buffer2;
    {
      std::istringstream iss(buffer);
      iss >> frame;
    }

    // Next comes the line with ABG data
    getlineNoComment(is, buffer);
    {
      std::istringstream iss(buffer);
      read(iss);
    }

    // Next is sigmas line and the correlation matrix
    DVector sigma(6);
    getlineNoComment(is, buffer);
    {
      std::istringstream iss(buffer);
      iss >> sigma[2];
      sigma[2] *= gamma()*gamma();
      iss >> sigma[0] >> sigma[1];
      sigma[0] *= ARCSEC;
      sigma[1] *= ARCSEC;
      iss >> sigma[3] >> sigma[4] >> sigma[5];
      sigma[3] *= ARCSEC/HOUR;
      sigma[4] *= ARCSEC/HOUR;
      sigma[5] *= ARCSEC/HOUR;
    }

    for (int i=0; i<6; i++) {
      getlineNoComment(is, buffer);
      {
	std::istringstream iss(buffer);
	for (int j=0; j<6; j++) {
	  iss >> covar(i,j);
	  covar(i,j) *= sigma[i]*sigma[j];
	}
      }
    }
  }

  ////////////////////////////////////////////////////////////////////
  // Routines for orbit determinations
  ////////////////////////////////////////////////////////////////////


  NBodyOrbit::NBodyOrbit(const StateCoords& xv0_, const UT& t0_,
			 double tStep_, JPL* eph):
    NumericalOrbit(xv0_,t0_), ephem(eph), tStep(tStep_) {flush();}

  void
  NBodyOrbit::flush() const {
    statedeque.clear();
  }

  // acceleration due to Newtonian gravity.
  // ??? no checks for close encounters right now.
  CartesianICRS
  NBodyOrbit::acceleration(CartesianICRS x,
			   const UT& t) const {
    static Planet bodies[]={SUN, JUPITER, SATURN, URANUS, NEPTUNE};
    static int nbodies = 5; 
    static double mass[]={1.000006, 9.548e-4, 2.859e-4, 4.37e-5, 5.15e-5};
    /* Note Sun mass includes the inner planets*/
    
    CartesianICRS xb, accel;

    for (int i=0; i<nbodies; i++) {
      xb = ephem->position(bodies[i], t);
      xb -= x;
      double r2=xb[0]*xb[0]+xb[1]*xb[1]+xb[2]*xb[2];
      double factor= mass[i] * pow(r2,-1.5);
      xb *= factor;
      accel += xb;
    }
    accel *= GM;
    return accel;
  }

  void 
  NBodyOrbit::integrateTo(const UT& t) const {
    if (statedeque.empty()) {
      // Set up time steps 0 and -1 to straddle t0
      StateCoords step0=xv0;
      StateCoords step1=xv0;
      CartesianICRS dv=acceleration(step0.x, t0);
      dv *= 0.5*tStep;
      step0.v += dv;
      statedeque.push_back(step0);
      tLast = 0.;

      // Step backwards in time
      step1.v -= dv;
      CartesianICRS dx=step1.v;
      dx *= -tStep;
      step1.x += dx;
      statedeque.push_front(step1);
      tFirst = -tStep;
    }

    int direction;
    Assert(tFirst<=tLast);
    double dTT = t-t0;

    while (dTT < tFirst) {
      // Add backwards steps
      double dT = -tStep;
      StateCoords xv=statedeque.front();
      UT ta=t0;
      ta += tFirst;
      CartesianICRS dv=acceleration(xv.x, ta);
      dv *= dT;
      xv.v += dv;
      CartesianICRS dx=xv.v;
      dx *= dT;
      xv.x += dx;
      statedeque.push_front(xv);
      tFirst += dT;
    }
      
    while (dTT>tLast) {
      double dT = +tStep;
      StateCoords xv=statedeque.back();
      CartesianICRS dx=xv.v;
      dx *= dT;
      xv.x += dx;

      tLast += dT;
      UT ta = t0;
      ta += tLast;
      CartesianICRS dv=acceleration(xv.x, ta);
      dv *= dT;
      xv.v += dv;

      statedeque.push_back(xv);
    }
     
    return;
  }
    
  StateCoords 
  NBodyOrbit::interpolate(const UT& t) const {
    integrateTo(t);	//insure that table extends this far

    Assert ( abs( (statedeque.size()-1)*tStep - (tLast-tFirst))
	     <= SECOND);

    double dTT = t - t0;
    int nearStep = static_cast<int> (rint((dTT-tFirst)/tStep));
    Assert(nearStep>=0 && nearStep<statedeque.size());

    // ??? add some static laststep to reduce searching
    UT ta = t0;
    ta += nearStep*tStep + tFirst;
    StateCoords xv=statedeque[nearStep];
    double dT = t - ta;

    CartesianICRS a=acceleration(xv.x, ta);
    // Recall that v applies to time + 0.5*tStep
    for (int i=0; i<3; i++) {
      xv.x[i] += dT*(xv.v[i] + a[i]*(dT/2-0.5*tStep));
      xv.v[i] += a[i]*(dT-0.5*tStep);
    }

    return xv;
  }
 
  CartesianICRS
  NBodyOrbit::position(const UT& t) const {
    return interpolate(t).x;
  }
  StateCoords
  NBodyOrbit::state(const UT& t) const {
    return interpolate(t);
  }

  void
  NBodyOrbit::makeDerivs(const UT& t, DMatrix& derivs,
			 StateCoords exact, bool doVelocity) const {
    if (doVelocity) {
      Assert(derivs.nrows()==6 && derivs.ncols()==6);
    } else {
      Assert(derivs.nrows()==3 && derivs.ncols()==6);
    }
    double dT = t-t0;
    // Find non-inertial component of motion
    CartesianICRS inertial=xv0.v;
    inertial *= dT;
    inertial += xv0.x;
    CartesianICRS accel=exact.x;
    accel -= inertial;
    
    double ir2 = pow(xv0.x.radius(), -2.);
    derivs(0,0) = 1. - 2.*xv0.x[0]*ir2*accel[0];
    derivs(0,1) =     -2.*xv0.x[1]*ir2*accel[0];
    derivs(0,2) =     -2.*xv0.x[2]*ir2*accel[0];
    derivs(0,3) = dT;
    derivs(0,4) = 0.;
    derivs(0,5) = 0.;

    derivs(1,0) =     -2.*xv0.x[0]*ir2*accel[1];
    derivs(1,1) = 1.  -2.*xv0.x[1]*ir2*accel[1];
    derivs(1,2) =     -2.*xv0.x[2]*ir2*accel[1];
    derivs(1,3) = 0.;
    derivs(1,4) = dT;
    derivs(1,5) = 0.;

    derivs(2,0) =     -2.*xv0.x[0]*ir2*accel[2];
    derivs(2,1) =     -2.*xv0.x[1]*ir2*accel[2];
    derivs(2,2) = 1.  -2.*xv0.x[2]*ir2*accel[2];
    derivs(2,3) = 0.;
    derivs(2,4) = 0.;
    derivs(2,5) = dT;

    if (!doVelocity) return;

    // Find total acceleration:
    accel=exact.v;
    accel -= xv0.v;

    derivs(3,0) =      -2.*xv0.x[0]*ir2*accel[0];
    derivs(3,1) =      -2.*xv0.x[1]*ir2*accel[0];
    derivs(3,2) =      -2.*xv0.x[2]*ir2*accel[0];
    derivs(3,3) = 1.;
    derivs(3,4) = 0.;
    derivs(3,5) = 0.;

    derivs(4,0) =      -2.*xv0.x[0]*ir2*accel[1];
    derivs(4,1) =      -2.*xv0.x[1]*ir2*accel[1];
    derivs(4,2) =      -2.*xv0.x[2]*ir2*accel[1];
    derivs(4,3) = 0.;
    derivs(4,4) = 1.;
    derivs(4,5) = 0.;

    derivs(5,0) =      -2.*xv0.x[0]*ir2*accel[2];
    derivs(5,1) =      -2.*xv0.x[1]*ir2*accel[2];
    derivs(5,2) =      -2.*xv0.x[2]*ir2*accel[2];
    derivs(5,3) = 0.;
    derivs(5,4) = 0.;
    derivs(5,5) = 1.;
  }

  CartesianICRS 
  NBodyOrbit::position(const UT& t, DMatrix& derivs) const {
    StateCoords exact=state(t);
    makeDerivs(t, derivs, exact, false);
    return exact.x;
  }
  StateCoords
  NBodyOrbit::state(const UT& t, DMatrix& derivs) const {
    StateCoords exact=state(t);
    makeDerivs(t, derivs, exact, true);
    return exact;
  }

  CartesianICRS
  InertialOrbit::position(const UT& t) const {
    CartesianICRS x=xv0.x;
    CartesianICRS dx=xv0.v;
    dx *= (t-t0);
    x+=dx;
    return x;
  }
  StateCoords
  InertialOrbit::state(const UT& t) const {
    StateCoords xv=xv0;
    xv.x = position(t);
    return xv;
  }

  void
  InertialOrbit::makeDerivs(const UT& t, DMatrix& derivs,
			    bool doVelocity) const {
    if (doVelocity) {
      Assert(derivs.nrows()==6 && derivs.ncols()==6);
    } else {
      Assert(derivs.nrows()==3 && derivs.ncols()==6);
    }

    derivs = 0.;
    double dT = t-t0;
    derivs(0,0) = 1.;
    derivs(0,3) = dT;

    derivs(1,1) = 1.;
    derivs(1,4) = dT;

    derivs(2,2) = 1.;
    derivs(2,5) = dT;

    if (!doVelocity) return;

    derivs(3,3) = 1.;
    derivs(4,4) = 1.;
    derivs(5,5) = 1.;
  }

  CartesianICRS 
  InertialOrbit::position(const UT& t, DMatrix& derivs) const {
    makeDerivs(t, derivs, false);
    return position(t);
  }
  StateCoords
  InertialOrbit::state(const UT& t, DMatrix& derivs) const {
    makeDerivs(t, derivs, true);
    return state(t);
  }


  // Get the apparent astrometric position of an orbiting target
  // given the position of the observer and a time of observation.
  // This takes light-travel time into account (to 1st order in v/c for
  // the derivatives) but currently does NOT account for differential
  // gravitational lensing.  The differential lensing is <1 mas unless
  // objects are <10 AU away or we look at solar elongations <90 degrees.

  double
  NumericalOrbit::timeDelay(const CartesianICRS& observer,
			    const UT& obstime) const {
    // First solve for the time of light departure
    const double DelayTolerance=0.001*SECOND;
    const int MaxIterations=10;

  
    double delay=0.;
    double olddelay=0.;
    UT depart=obstime;
    int iter=0;
    do {
      CartesianICRS x1=position(depart);
      x1 -= observer;
      delay = x1.radius() / SpeedOfLightAU;
      if ( abs(delay-olddelay)< DelayTolerance) break;
      olddelay = delay;
      depart = obstime;
      depart -= delay;
      iter++;
    } while (iter<MaxIterations);

    if (iter>=MaxIterations) 
      throw AstrometryError("Too many time-delay iterations"
			    " in astrometicPosition");
    return delay;
  }

  CartesianICRS
  NumericalOrbit::astrometricPosition(const CartesianICRS& observer,
				      const UT& obstime) const {
    UT depart=obstime;
    depart -= timeDelay(observer,obstime);
    CartesianICRS dx=position(depart);
    dx -= observer;
    return dx;
  }

  // Also return derivatives w.r.t. 6-d initial StateCoord
  CartesianICRS
  NumericalOrbit::astrometricPosition(const CartesianICRS& observer,
				      const UT& obstime,
				      DMatrix& derivs) const {
    UT depart=obstime;
    double delay= timeDelay(observer,obstime);
    depart -= delay;
    DMatrix dxv(6,6);
    StateCoords xv1=state(depart, dxv);
    xv1.x -= observer;
    // Derivatives of apparent position w.r.t. true state vector:
    DMatrix dposn(3,6,0.);
    for (int i=0; i<3; i++) dposn(i,i)=1.;
    // And time-delay terms:
    for (int i=0; i<3; i++) {
      for (int j=0; j<3; j++)
	dposn(i,j) -= xv1.v[i]*xv1.x[j]/(SpeedOfLightAU*xv1.x.radius());
    }
    derivs = dposn * dxv;
    return xv1.x;
  }


  ///////////////////////////////////////////////////////////////////////////////


    /* Orbiting observatory  calculation
    double pole, rapole, phase;
    int i;
    SPACECRAFT *s;
    for (i=0; obscode!=spacecraftlis[ti].code && i<nspacecraft; i++)  ;
    if (i>=nspacecraft) {
      fprintf(stderr,"Unknown spacecraft code %d\n",obscode);
      exit(1);
    }
    s = &(spacecraftlist[i]);
    pole = s->i;
    phase = -(jd - s->jd0) / s->precess;//
					// for precess > 0
    phase -= floor(phase);
    rapole = s->ra0 + phase*TPI;

    phase = (jd - s->jd0) / s->P;	// orbit is increasing RA 
					// for P > 0 
    phase -= floor(phase);
    phase *= TPI;

    *xobs = s->a * (cos(phase)*cos(rapole) 
		    - cos(pole)*sin(phase)*sin(rapole));
    *yobs = s->a * (cos(phase)*sin(rapole) 
		    + cos(pole)*sin(phase)*cos(rapole));
    *zobs = s->a * sin(phase)*sin(pole);
  }
  return;
  }
    */


/* Return the angle from zenith to the horizon for this observatory
double
zenith_horizon(int obscode) {
  if (obscode < OBSCODE_ORBITAL)
    return PI/2.;
  else {
    // Do spherical-Earth limb calculation for orbiting spacecraft
    double a;
    int i;
    if (nsites<=0 && nspacecraft<=0) read_observatories(NULL); 
    for (i=0; obscode!=spacecraftlist[i].code && i<nspacecraft; i++)  ;
    if (i>=nspacecraft) {
      fprintf(stderr,"Unknown spacecraft code %d\n",obscode);
      exit(1);
    }
    a = spacecraftlist[i].a;
    a /= EQUAT_RAD/(1000.*R1.AU);
    if (a < 1.) {
      fprintf(stderr,"Your spacecraft is underground.  Oops.\n");
      exit(1);
    }
    return PI - asin(1./a);
  }
  }*/
  /* And likewise some information for orbiting observatories: */
  static int nspacecraft=0;
  typedef struct {
    int code;
    double	i;		/*inclination of orbit (degrees)*/
    double	P;		/*sidereal orbital period, (days)*/
    double	precess;	/*period for orbital precession (days)*/
    double	jd0;		/*time of zero orbit phase*/
    double	ra0;		/*RA of orbit pole at jd0 (degrees)*/
    double	a;		/*semi-major axis of orbit (m)*/
    char	name[80];
  } SPACECRAFT;


  // Get orbital elements from phase space coordinates.
  // Both are barycentric ICRS system.
  // Taken essentially directly from Bharat Khushalani's code.
  void 
  OrbitalElements::fromXV(const StateCoords& xv, const UT& UT0) {

    double combinedMass; /* mass of Sun + mass of KBO */ 

    double R[4], V[4];
    double rMagnitude, velSquare, radVelDotProduct;
    double eccentricityVector[4], angularMomentum[4], ascendingNode[4];
    double semiLatusRectum;
    double hMagnitude, ascendingNodeMagnitude; /* magnitude of angular momentum */
    double ascEccDotProduct;
    double xBar, yBar;
    double cosE, sinE, E1, E2, eccentricAnomaly;
    /* E1 and E2 are used to decide the quadrant of Eccentric Anomaly */
    double meanMotion;

    if (heliocentric) combinedMass = 4*PI*PI*(1.0000137398);
    combinedMass = GM*SolarSystemMass ; 

    // Use total SS mass & barycentric coords as approximation for outer
    // solar-system bodies.

    CartesianEcliptic xx(xv.x);
    CartesianEcliptic vv(xv.v);
    R[1] = xx[0]; 
    R[2] = xx[1]; 
    R[3] = xx[2]; 
    V[1] = vv[0]; 
    V[2] = vv[1]; 
    V[3] = vv[2]; 

    rMagnitude = sqrt(R[1]*R[1] + R[2]*R[2] + R[3]*R[3]);
    velSquare = V[1]*V[1] + V[2]*V[2] + V[3]*V[3];
    radVelDotProduct = R[1]*V[1] + R[2]*V[2] + R[3]*V[3];
    for (int k=1;k<=3;k++)
      {
	eccentricityVector[k] = (velSquare/combinedMass - 1/rMagnitude)*R[k]
	  - (radVelDotProduct/combinedMass)*V[k];
      }

    /* Angular mom is cross product of rad and vel */
    angularMomentum[1] = R[2]*V[3] - R[3]*V[2];
    angularMomentum[2] = R[3]*V[1] - R[1]*V[3];
    angularMomentum[3] = R[1]*V[2] - R[2]*V[1];

    /* Ascending Node vector is cross product of k-unit vector and angular momentum
       vector. k = [0 0 1] and h = [hx, hy, hz] */
    ascendingNode[1] = -angularMomentum[2];
    ascendingNode[2] = angularMomentum[1];
    ascendingNode[3] = 0.0;

    a = 1/(2/rMagnitude - velSquare/combinedMass);
    e = sqrt(eccentricityVector[1]*eccentricityVector[1] + 
			eccentricityVector[2]*eccentricityVector[2] + 
			eccentricityVector[3]*eccentricityVector[3]);
    semiLatusRectum = (angularMomentum[1]*angularMomentum[1] + 		    
		       angularMomentum[2]*angularMomentum[2] +	
		       angularMomentum[3]*angularMomentum[3])/combinedMass;
    /* p = h-square by mu */
    hMagnitude = sqrt( angularMomentum[1]*angularMomentum[1] + 		    
		       angularMomentum[2]*angularMomentum[2] +	
		       angularMomentum[3]*angularMomentum[3] );
    i = acos(angularMomentum[3]/hMagnitude); /* in radians here */
    ascendingNodeMagnitude = sqrt(ascendingNode[1]*ascendingNode[1] +
				  ascendingNode[2]*ascendingNode[2] +
				  ascendingNode[3]*ascendingNode[3]);
    longitudeOfAscendingNode = acos(ascendingNode[1]/ascendingNodeMagnitude);
    /* Capital Omega in radians here */
    if (ascendingNode[2] < 0) longitudeOfAscendingNode = 
				2*PI - longitudeOfAscendingNode;
    /* ???could use atan2 here?? */
    ascEccDotProduct = ascendingNode[1]*eccentricityVector[1] +
      ascendingNode[2]*eccentricityVector[2] +
      ascendingNode[3]*eccentricityVector[3];
    argumentOfPerifocus = acos(ascEccDotProduct/
			       (ascendingNodeMagnitude*e)); 
    /* Small omega in radians here */
    if (eccentricityVector[3] < 0) argumentOfPerifocus = 
				     2*PI - argumentOfPerifocus;
    xBar = (semiLatusRectum - rMagnitude)/e;
    yBar = radVelDotProduct*sqrt(semiLatusRectum/combinedMass)/e;

    /* From here, we assume that the motion is elliptical */

    cosE = (xBar/a) + e;
    sinE = yBar/(a*sqrt(1-e*e));
    /* where a*sqrt(1-e*e) is semiminor */
    eccentricAnomaly = atan2(sinE,cosE);

    meanAnomaly = eccentricAnomaly - e*sinE; /* radians */
    epoch = UT0;
  }

  void 
  OrbitalElements::toXV(StateCoords& xv, const UT& UT0) const {

    double eccentricAnomaly, r0[3], v0[3], r1[3], v1[3], r2[3], v2[3];
    double c, s, t, dt;

    double mu = GM*SolarSystemMass;	/*use SSMASS, work in AU/yrs*/
    if (heliocentric) mu = 4*PI*PI*(1.0000137398);

    if (e >= 1. || a <=0.) {
      cerr << "elements_to_xv only for closed orbits now" << endl;
      exit(1);
    }

    /* get the eccentric Anomaly from the current mean anomaly */
    double currentMA;
    currentMA = (UT0 - epoch) * pow(a,-1.5) * sqrt(mu)
      + meanAnomaly;
    /* Put it near 0 */
    t = floor (currentMA / TPI);
    currentMA -= t*TPI;
    /* Use bisection to find solution (adapt Numerical Recipes)*/  
    {
#define JMAX 40
#define TOLERANCE SECOND
      double f, fmid, x1, x2, dx, xmid, rtb;
      int j;
      x1 = currentMA - e;
      x2 = currentMA + e;
      f   = x1 - e * sin(x1) - currentMA;
      fmid= x2 - e * sin(x2) - currentMA;
      if (f*fmid > 0.0) {
	cerr << "Error, eccentricAnomaly root not bracketed" << endl;
	cerr << "f " << f << " fmid " << fmid << endl;
	exit(1);
      }

      rtb = f < 0.0 ? (dx=x2-x1,x1) : (dx=x1-x2,x2);
      for (j=1;j<=JMAX;j++) {
	xmid=rtb+(dx *= 0.5);
	fmid= xmid - e * sin(xmid) - currentMA;
	if (fmid <= 0.0) rtb=xmid;
	if (fabs(dx) < TOLERANCE || fmid == 0.0) break;
      }
      if (j>=JMAX) {
	fprintf(stderr,"eccentricAnomaly took too long\n");
	exit(1);
      }
      currentMA = rtb;
#undef JMAX
#undef TOLERANCE   
    }

    /*Coordinates and velocity in the system aligned with orbit: */
    c = cos(currentMA);
    s = sin(currentMA);
    r0[0] = a * (c - e);
    r0[1] = a * s * sqrt(1-e*e);
    dt = sqrt(pow(a,3.)/mu) * ( 1 - e*c);
    v0[0] = -a * s / dt;
    v0[1] = a * c * sqrt(1-e*e) / dt;

    /* Rotate about z to put perihelion at arg of peri */
    c = cos(argumentOfPerifocus);  s = sin(argumentOfPerifocus);
    r1[0] = r0[0]*c - r0[1]*s;
    r1[1] = r0[0]*s + r0[1]*c;
    v1[0] = v0[0]*c - v0[1]*s;
    v1[1] = v0[0]*s + v0[1]*c;

    /* Rotate about x axis to incline orbit */
    c = cos(i);  s = sin(i);
    r2[0] = r1[0];
    r2[1] = r1[1]*c;
    r2[2] = r1[1]*s;
    v2[0] = v1[0];
    v2[1] = v1[1]*c;
    v2[2] = v1[1]*s;

    /* Rotate about z axis to bring node to longitude */
    c = cos(longitudeOfAscendingNode);  s = sin(longitudeOfAscendingNode);
    CartesianEcliptic xx;
    CartesianEcliptic vv;
    xx[0] = r2[0]*c - r2[1]*s;
    xx[1] = r2[0]*s + r2[1]*c;
    xx[2] = r2[2];
    vv[0] = v2[0]*c - v2[1]*s;
    vv[1] = v2[0]*s + v2[1]*c;
    vv[2] = v2[2];

    xv.x = CartesianICRS(xx);
    xv.v = CartesianICRS(vv);
    SphericalICRS posn(xv.x);
  }

  OrbitalElements::OrbitalElements(const ABGState& agb, bool helio):
    heliocentric(helio) {
    fromXV(agb.xv(), agb.t0());
  }
  OrbitalElements::OrbitalElements(const StateCoords& xv, const UT& UT0,
				   bool helio):
    heliocentric(helio) {
    fromXV(xv, UT0);
  }
  StateCoords 
  OrbitalElements::getXV(const UT& UT0) const {
    StateCoords xv;
    toXV(xv,UT0);
    return xv;
  }

  bool
  OrbitalElements::read(std::istream& is, bool helio) {
    heliocentric = helio;
    is >> epoch >> a >> e >> i 
       >> longitudeOfAscendingNode >> argumentOfPerifocus
       >> meanAnomaly;
    i *= DEGREE;
    longitudeOfAscendingNode *= DEGREE;
    argumentOfPerifocus *= DEGREE;
    meanAnomaly *= DEGREE;
    return (is);
  }

  void
  OrbitalElements::write(std::ostream& os) const {
    os << std::setprecision(6) << std::fixed << a
       << " " << e
       << " " << i/DEGREE
       << " " << longitudeOfAscendingNode/ DEGREE
       << " " << argumentOfPerifocus / DEGREE
       << " " << meanAnomaly / DEGREE
       << " " << epoch.getJD();
  }

  void
  OrbitalElements::writeHeader(std::ostream& os) const {
    os << "    a    "
       << "    e    "
       << "    i    "
       << "  AscNode  "
       << "  ArgPeri  "
       << "  meanAnomaly "
       << "     jd0 ";
  }
} // namespace ephem
