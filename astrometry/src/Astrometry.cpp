#include "Astrometry.h"

#include <string>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <vector>
using namespace std;
#include <cmath>

#include "LeapSeconds.h" //This has the TAIminusUTC function/table.
#include "StringStuff.h"

namespace astrometry {
  // function to give difference between TT (dynamical) time and UTC:
  double
  UT::TTminusUTC(double jd_) {
    const double TTminusTAI=32.184*TIMESEC;
    return TTminusTAI + TAIminusUTC(jd_);
  }
      
  // Use (continuous) TT for intervals instead of (jumpy) UT
  const UT&
  UT::operator+=(double dt) {
    setTT(getTT()+dt/DAY); 
    return *this;
  }
  const UT&
  UT::operator-=(double dt) {
    setTT(getTT()-dt/DAY); 
    return *this;
  }
  double
  UT::operator-(const UT& rhs) const {
    return (getTT() - rhs.getTT()) * DAY;
  }
  // Helpful routines for I/O of sexagesimal quantities
  // Convert string to decimal degrees
  double 
  dmsdeg(const string& s)
  {
    vector<double> args(3,0.);
    string sub;
    int i,iarg,j,sgn;

    /*strip off leading whitespace and look for - sign*/
    for (i=0; i<s.size() && isspace(s[i]); i++);
    if (i>=s.size()) 
      throw AstrometryError("Empty input string to dmsdeg");
    if (s[i]=='-') {
      sgn = -1;
      i++;
    }
    else {
      sgn = 1;
    }

    for (int iarg=0; iarg<3; iarg++) {
      string sub(s,i);
      std::istringstream is(sub);
      is >> args[iarg];
      if (!is)
	throw AstrometryError("DMS format error in string <" 
			      + s + ">");
      size_t j=s.find_first_not_of("0123456789.+",i);
      if (j!=string::npos) j=s.find_first_not_of(" :",j);
      if (j==string::npos) break;
      i=j;
    };
    args[1] += args[2]/60.;
    args[0] += args[1]/60.;
    return( args[0]*sgn );
  }

  double 
  hmsdeg(const string& s)
  {
    return(15. * dmsdeg(s));
  }

  string
  dmsWriter(double degr, int decimalPlaces, bool isDegrees)
  {
    std::ostringstream os;
    int sign;
    if (degr<0.) {
      sign=-1;
      degr *= -1.;
    } else {
      sign=+1;
    }

    int ideg = static_cast<int> (floor(degr));
    degr = (degr-ideg)*60.;
    int imin = static_cast<int> (floor(degr));
    degr = (degr-imin)*60.;

    // Check for rounding up seconds to 60s
    if (degr > 60. - 0.5*pow(10,-1.*decimalPlaces)) {
      degr = 0.;
      imin += 1;
      if (imin>=60) {
	imin -= 60;
	ideg += 1;
      }
    }
    if (sign>0) 
      if (isDegrees) {
	if (ideg<100)
	  os << " +" << std::setw(2) << std::setfill('0') << ideg;
	else
	  os << "+" << std::setw(3) << std::setfill('0') << ideg;
      } else {
	// Hours - leave off plus sign is usual convention
	os << std::setw(2) << std::setfill('0') << ideg;
      }
    else
      if (isDegrees) {
	if (ideg<100)
	  os << " -" << std::setw(2) << std::setfill('0') << ideg;
	else
	  os << "-" << std::setw(3) << std::setfill('0') << ideg;
      } else {
	// hours: negative is an aberration.
	os << "-" << std::setw(2) << std::setfill('0') << ideg;
      }
    os << ":" 
       << std::setw(2) << std::setfill('0') << imin
       << ":" 
       << std::setprecision(decimalPlaces) << std::setfill('0')
       << std::fixed << std::setw(decimalPlaces+3) << degr;
    return os.str();
  }

  string
  degdms(double degr, int decimalPlaces) {
    return dmsWriter(degr, decimalPlaces, true);
  }
  string
  deghms(double degr, int decimalPlaces) {
    return dmsWriter(degr/15., decimalPlaces, false);
  }

  bool hasDot(string& s) {
    return s.find('.')!=string::npos;
  }

  ///////////////////////////////////////////////////////////////
  // Manipulation of astronomical (UT) time/date stamps
  ///////////////////////////////////////////////////////////////

  void
  UT::set(double jd_) {
    jd=jd_; 
    ymdValid=false;
  }
  // ??? Add validation of date specs to below???
  void
  UT::set(int y_, int m_, double d_) {
    if (m_<1 || m_>12)
      throw AstrometryError("Bad month specification");
    if (d_<1. || d_>32.)
      throw AstrometryError("Bad day specification");
    y=y_; mo=m_; d=d_;
    ymdValid=true;
    buildJD();
  }
  void
  UT::set(int y_, int m_, int d_, int h_, int min_, double s_) {
    if (m_<1 || m_>12)
      throw AstrometryError("Bad month specification");
    if (d_<1 || d_>31)
      throw AstrometryError("Bad day specification");
    if (h_<0 || h_>23)
      throw AstrometryError("Bad hour specification");
    if (min_<0 || min_>59)
      throw AstrometryError("Bad minute specification");
    if (s_<0. || s_>61.) //note leap seconds allowed
      throw AstrometryError("Bad second specification");
    y=y_; mo=m_; 
    d=d_ + (h_+(min_+s_/60.)/60.)/24.;
    /**/cerr << "y,m,d: " << y << " " << mo << " " << d << endl;
    ymdValid=true; 
    buildJD();
  }

  void
  UT::buildJD() {
    if (!ymdValid)
      throw AstrometryError("Attempt to find JD for uninitialized UT");
    // Following is taken from Skycalc:
    /* From Meeus' Astronomical Formulae for Calculators. */
    int y0, m0;
    long A, B;
  
    if(mo <= 2) {
      y0 = y - 1;
      m0 = mo + 12;
    }
    else {
      y0 = y;
      m0 = mo;
    }

    A = (long) (y0 / 100.);
    B = 2 - A + (long) (A / 4.);

    jd = (long) (365.25 * y0) + (long) (30.6001 * (m0 + 1)) + d +
      1720994.5;

    if(y > 1583) jd += B;
    /* Not quite right, since Gregorian calendar first
       adopted around Oct 1582.  But fine for modern. */
  }

  void
  UT::buildYMD() const {
    if (ymdValid) return;
    if (jd==0.)
      throw AstrometryError("Attempt to find YMD for uninitialized UT");
    // Following is taken from Skycalc:
    /* from Jean Meeus, Astronomical Formulae for Calculators,
       published by Willman-Bell Inc.
       Avoids a copyrighted routine from Numerical Recipes.
       Tested and works properly from the beginning of the
       Gregorian calendar era (1583) to beyond 3000 AD. */

    double jdtmp;
    long alpha;
    long Z;
    long A, B, C, D, E;
    double F;
    int intd;

    jdtmp = jd + 0.5;
    Z = (long) jdtmp;

    F = jdtmp - Z;

    if(Z < 2299161) A = Z;
    else {
      alpha = (long) ((Z - 1867216.25) / 36524.25);
      A = Z + 1 + alpha - (long) (alpha / 4);
    }

    B = A + 1524;
    C = ((B - 122.1) / 365.25);
    D =  (365.25 * C);
    E =  ((B - D) / 30.6001);

    intd = B - D - (long)(30.6001 * E);
    if(E < 13.5) mo = E - 1;
    else mo = E - 13;
    if(mo  > 2.5)  y = C - 4716;
    else y = C - 4715;
  
    d = intd + F;
  }
  
  void 
  UT::getYMDHMS(int& y_, int& m_, int& d_, 
		int& h, int& mn, double& s) const {
    buildYMD(); 
    y_=y; m_=mo; 
    d_ = static_cast<int> (floor(d));
    double remainder=d-d_;
    remainder *= 24.;
    h = static_cast<int> (floor(remainder));
    remainder -= h;
    remainder *= 60.;
    mn = static_cast<int> (floor(remainder));
    remainder -= mn;
    s = remainder*60.;
  }

  void
  UT::writeYMD(ostream& os) const {
    buildYMD();
    char oldfill = os.fill('0');
    os << setw(4) << std::right << y
       << " " << setw(2) << std::right << mo
       << " " << setw(os.precision()+3) << d; // note 1s = 1e-5 day
    os << std::setfill(oldfill);
  }

  void
  UT::writeYMDHMS(ostream& os) const {
    int yy, mm, dd, hh, mn;
    double s;
    getYMDHMS(yy,mm,dd,hh,mn,s);
    double h=hh + (mn +s/60.)/60.;
    string hs=deghms(15.*h, 1);
    char oldfill = os.fill('0');
    os << setw(4) << std::right << y
       << " " << setw(2) << std::right << mo
       << " " << setw(2) << std::right << dd
       << std::setfill(oldfill)
       << " " << hs;
  }

  ostream& operator<<(ostream& os, const UT& rhs) {
    rhs.writeYMD(os);
    return os;
  }
  istream& operator>>(istream& is, UT& rhs) {
    rhs.read(is);
    return is;
  }

  // Read a date+time from stream.  Following forms are accepted, number
  // of digits is alway irrelevant
  // jd  (can tell by large number here)
  // year mo dd.dddd  (decimal place on date is required)
  // year mo dd hh:mm:ss.sss (time spec must be present but can be truncated).
  void
  UT::read(istream& is) {
    double jd1;
    string buff;
    ymdValid = false;
    jd=0.;
    if (!(is >> buff)) return;
    {
      std::istringstream iss(buff);
      if (!(iss >> jd1)) {
	is.setstate(is.rdstate() | std::ios_base::failbit);
	return;
      }
      if (hasDot(buff) || jd1>10000.) {
	// This must be a JD, not a year
	set(jd1);
	return;
      }
    }
    // Must be in YMD format then
    y = static_cast<int> (jd1);
    if (!(is >> mo))   return;
    if (!(is >> buff)) return;
    if ( hasDot(buff)) {
      // Decimal date indicates this is fractional-date format
      std::istringstream iss(buff);
      if (!(iss >> d)) {
	is.setstate(is.rdstate() | std::ios_base::failbit);
	return;
      }
      // ?? validate date??
    } else {
      // Integer date means a time follows
      int dd;
      {
	std::istringstream iss(buff);
	if (!(iss >> dd)) {
	  is.setstate(is.rdstate() | std::ios_base::failbit);
	  return;
	}
	// ?? validate date??
      }
      if (!(is >> buff)) return;
      double deg;
      try {
	deg=hmsdeg(buff);
      } catch (AstrometryError& e) {
	// Format error in the hms read
	is.setstate(is.rdstate() | std::ios_base::failbit);
	return;
      }
      d = dd + deg/360.;
    }
    ymdValid= true;
    buildJD();
    return;
  }

  ////////////////////////////////////////////////////////////////
  //  3d Cartesian coordinate systems
  ////////////////////////////////////////////////////////////////


  // First some common 3d rotations:
  // Produce new Cartesian coords about a new coord system that
  // has z axis at inclination to current z axis, and its x axis
  // along ascending node at given longitude in THIS system.
  Vector3
  rotateToPole(const Vector3& xin,
	       double inclination,
	       double ascendingNode,
	       Matrix33* partials) {
    double sn=sin(ascendingNode);
    double cn=cos(ascendingNode);
    double si=sin(inclination);
    double ci=cos(inclination);
    Vector3 dest;
    dest[0] =  xin[0]*cn    + xin[1]*sn;
    dest[1] = -xin[0]*ci*sn + xin[1]*ci*cn + xin[2]*si;
    dest[2] =  xin[0]*si*sn - xin[1]*si*cn + xin[2]*ci;

    if (partials) {
      (*partials)(0,0) = cn;
      (*partials)(0,1) = sn;
      (*partials)(0,2) = 0.;
      (*partials)(1,0) = -ci*sn;
      (*partials)(1,1) = ci*cn;
      (*partials)(1,2) = si;
      (*partials)(2,0) = si*sn;
      (*partials)(2,1) = -si*cn;
      (*partials)(2,2) = ci;
    }
    return dest;
  }

  Vector3
  rotateFromPole(const Vector3& xin,
		 double inclination,
		 double ascendingNode,
		 Matrix33* partials){
    double sn=sin(ascendingNode);
    double cn=cos(ascendingNode);
    double si=sin(inclination);
    double ci=cos(inclination);
    Vector3 dest;
    dest[0] =  xin[0]*cn - xin[1]*ci*sn + xin[2]*si*sn;
    dest[1] =  xin[0]*sn + xin[1]*ci*cn - xin[2]*si*cn;
    dest[2] =              xin[1]*si    + xin[2]*ci;

    if (partials) {
      (*partials)(0,0) = cn;
      (*partials)(1,0) = sn;
      (*partials)(2,0) = 0.;
      (*partials)(0,1) = -ci*sn;
      (*partials)(1,1) = ci*cn;
      (*partials)(2,1) = si;
      (*partials)(0,2) = si*sn;
      (*partials)(1,2) = -si*cn;
      (*partials)(2,2) = ci;
    }
    return dest;
  }

  // Caches for transformations to/from fixed frames
  Matrix33*
  SphericalEcliptic::icrs2ecliptic = nullptr;
  Matrix33*
  SphericalInvariable::icrs2invariable = nullptr;

  void
  SphericalEcliptic::setup() {
    if (icrs2ecliptic) return;
    // Use rotation routines to build ecliptic rotation
    icrs2ecliptic = new Matrix33;
    Vector3 dummy;
    rotateToPole(dummy,
		 EclipticInclination,
		 EclipticNode,
		 icrs2ecliptic);
  }

  void
  SphericalInvariable::setup() {
    if (icrs2invariable) return;
    // Use rotation routines to build ecliptic rotation
    icrs2invariable = new Matrix33;
    Vector3 dummy;
    rotateToPole(dummy,
		 InvariableInclination,
		 InvariableNode,
		 icrs2invariable);
  }
  
  ////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////
  void
  CartesianCoords::write(std::ostream& os) const {
    os << fixed << setprecision(9) << x[0]
       << " " << x[1]
       << " " << x[2];
  }
  void
  CartesianCoords::read(std::istream& is) {
    is >> x[0] >>  x[1] >> x[2];
  }
  std::ostream& operator<<(std::ostream& os, 
			   const CartesianCoords& rhs) {
    rhs.write(os);
    return os;
  }
  std::istream& operator>>(std::istream& is, 
			   CartesianCoords& rhs) {
    rhs.read(is);
    return is;
  }
  
  void
  CartesianCoords::convertFrom(const CartesianCoords& rhs) {
    Vector3 xeq = rhs.convertToICRS();
    this->convertFromICRS(xeq);
  }
  void
  CartesianCoords::convertFrom(const CartesianCoords& rhs,
			       Matrix33& partials) {
    Matrix33 partialRHS;
    Vector3 xeq = rhs.convertToICRS(&partialRHS);
    Matrix33 partial3d;
    this->convertFromICRS(xeq, &partial3d);
    partials = partial3d * partialRHS;
  }

  // Cartesian Coordinate transformations:
  void
  CartesianICRS::convertFromICRS(const Vector3& xeq, 
				     Matrix33* partials) {
    setVector(xeq);
    if (partials) partials->setToIdentity();
  }
  Vector3
  CartesianICRS::convertToICRS(Matrix33* partials) const {
    Vector3 xeq = getVector();
    if (partials) partials->setToIdentity();
    return xeq;
  }

  void
  CartesianEcliptic::convertFromICRS(const Vector3& xeq, 
				   Matrix33* partials) {
    SphericalEcliptic::setup();  // Check rotation matrix
    setVector( *SphericalEcliptic::icrs2ecliptic * xeq);
    if (partials) *partials = *SphericalEcliptic::icrs2ecliptic;
  }

  Vector3
  CartesianEcliptic::convertToICRS(Matrix33* partials) const {
    SphericalEcliptic::setup();  // Check rotation matrix
    if (partials) *partials = SphericalEcliptic::icrs2ecliptic->transpose();
    return SphericalEcliptic::icrs2ecliptic->transpose() * getVector();
  }

  void
  CartesianInvariable::convertFromICRS(const Vector3& xeq, 
				   Matrix33* partials) {
    SphericalInvariable::setup();  // Check rotation matrix
    setVector( *SphericalInvariable::icrs2invariable * xeq);
    if (partials) *partials = *SphericalInvariable::icrs2invariable;
  }

  Vector3
  CartesianInvariable::convertToICRS(Matrix33* partials) const {
    SphericalInvariable::setup();  // Check rotation matrix
    if (partials) *partials = SphericalInvariable::icrs2invariable->transpose();
    return SphericalInvariable::icrs2invariable->transpose() * getVector();
  }

  void
  CartesianCustom::convertFromICRS(const Vector3& xeq, 
				   Matrix33* partials) {
    x =  ref->orient.fromICRS(xeq - ref->origin.getVector());
    if (partials) *partials=ref->orient.m();
  }
  Vector3
  CartesianCustom::convertToICRS(Matrix33* partials) const {
    if (partials) *partials=ref->orient.m().transpose();
    return  ref->orient.toICRS(x) + ref->origin.getVector();
  }

  CartesianICRS::CartesianICRS(const CartesianCoords& rhs) {
    convertFrom(rhs);
  }
  CartesianICRS::CartesianICRS(const CartesianCoords& rhs, 
			       Matrix33& partials) {
    convertFrom(rhs, partials);
  }

  CartesianEcliptic::CartesianEcliptic(const CartesianCoords& rhs) {
    convertFrom(rhs);
  }
  CartesianEcliptic::CartesianEcliptic(const CartesianCoords& rhs, 
				       Matrix33& partials) {
    convertFrom(rhs, partials);
  }

  CartesianInvariable::CartesianInvariable(const CartesianCoords& rhs) {
    convertFrom(rhs);
  }
  CartesianInvariable::CartesianInvariable(const CartesianCoords& rhs, 
					   Matrix33& partials) {
    convertFrom(rhs, partials);
  }

  CartesianCustom::CartesianCustom(const ReferenceFrame& ref_,
				   bool shareFrame): ref(&ref_), ownFrame(false) {
    if (!shareFrame) {
      ref = new ReferenceFrame(ref_);
      ownFrame = true;
    } 
  }

  CartesianCustom::CartesianCustom(const Vector3& x_,
				   const ReferenceFrame& ref_,
				   bool shareFrame): ref(&ref_),
						     ownFrame(false),
						     CartesianCoords(x_) {
    if (!shareFrame) {
      ref = new ReferenceFrame(ref_);
      ownFrame = true;
    } 
  }
  
  CartesianCustom::CartesianCustom(const CartesianCoords& rhs,
				   const ReferenceFrame& ref_,
				   bool shareFrame): ref(&ref_),
						     ownFrame(false) {
    if (!shareFrame) {
      ref = new ReferenceFrame(ref_);
      ownFrame = true;
    } 
    convertFrom(rhs);
  }
  CartesianCustom::CartesianCustom(const CartesianCoords& rhs, 
				   const ReferenceFrame& ref_,
				   Matrix33& partials,
				   bool shareFrame): ref(&ref_),
						     ownFrame(false) {
    if (!shareFrame) {
      ref = new ReferenceFrame(ref_);
      ownFrame = true;
    } 
    convertFrom(rhs, partials);
  }

  ///////////////////////////////////////////////////////////////////////////
  // Spherical coordinate systems
  ///////////////////////////////////////////////////////////////////////////

  SphericalCoords::SphericalCoords(): x(0.) {}
  SphericalCoords::SphericalCoords(double lon, double lat): x() {
    setLonLat(lon,lat);
  }
  SphericalCoords::SphericalCoords(const Vector3& x_) {
    setUnitVector(x_);
  }
  void
  SphericalCoords::setUnitVector(const Vector3& x_) {
    x=x_;
  }
  void
  SphericalCoords::setLonLat(const Vector2& x_) {
    setLonLat(x_[0], x_[1]);
  }
  void
  SphericalCoords::setLonLat(double lon, double lat) {
    double clat=cos(lat);
    x[0] = clat*cos(lon);
    x[1] = clat*sin(lon);
    x[2] = sin(lat);
  }

  double
  SphericalCoords::distance(const SphericalCoords& rhs) const {
    Vector3 x1=convertToICRS();
    Vector3 x2=rhs.convertToICRS();
    x1 -= x2;
    double chord=0.5*sqrt(x1[0]*x1[0]+x1[1]*x1[1]+x1[2]*x1[2]);
    return 2*asin(chord);
  }
    
  Matrix32
  SphericalCoords::derivsFrom2d() const {
    Matrix32 dxdl;
    double rho=hypot(x[1],x[0]);
    if (rho<=0.)
      throw AstrometryError("Taking LatLon derivs at pole");
    dxdl(0,0) = -x[1];
    dxdl(0,1) = -x[2]*x[0]/rho;
    dxdl(1,0) =  x[0];
    dxdl(1,1) = -x[2]*x[1]/rho;
    dxdl(2,0) = 0.;
    dxdl(2,1) = rho;
    return dxdl;
  }

  Matrix23
  SphericalCoords::derivsTo2d() const {
    Matrix23 project;
    double rho=hypot(x[0],x[1]);
    if (rho<=0.)
      throw AstrometryError("getLatLon spherical coord derivs at pole");
    project(0,0) = -x[1]/(rho*rho);
    project(0,1) =  x[0]/(rho*rho);
    project(0,2) = 0.;
    project(1,0) = -x[0]*x[2]/rho;
    project(1,1) = -x[1]*x[2]/rho;
    project(1,2) =  rho;
    return project;
  }

  void
  SphericalCoords::getLonLat(double& lon, double& lat) const {
    // ?? normalize by radius to avoid rounding problems?
    lat = asin(x[2]);
    lon = atan2(x[1],x[0]);
  }

  Vector2
  SphericalCoords::getLonLat() const {
    Vector2 l;
    getLonLat(l[0],l[1]);
    return l;
  }

  // Default I/O is to put have lat/lon in hms/dms format
  // Decimal places gives fraction of arcsec
  void
  SphericalCoords::write(std::ostream& os, int decimalPlaces) const {
    double lon, lat;
    getLonLat(lon,lat);
    while (lon<0.) lon+= TPI;
    string buffer = deghms(lon/DEGREE, decimalPlaces+1);
    os << buffer;
    buffer = degdms(lat/DEGREE,decimalPlaces);
    os << " " << buffer;
  }
  void
  SphericalCoords::read(std::istream& is) {
    double lon, lat;
    string buffer;
    double deg;
    x.setZero();
    if (!(is >> buffer)) return;
    try {
      deg=hmsdeg(buffer);
    } catch (AstrometryError& e) {
      // Format error in the hms read
      is.setstate(is.rdstate() | std::ios_base::failbit);
      return;
    }
    lon = deg*DEGREE;
    if (!(is >> buffer)) return;
    try {
      deg=dmsdeg(buffer);
    } catch (AstrometryError& e) {
      // Format error in the hms read
      is.setstate(is.rdstate() | std::ios_base::failbit);
      return;
    }
    lat = deg*DEGREE;
    setLonLat(lon,lat);
  }

  void
  SphericalCoords::convertFrom(const SphericalCoords& rhs) {
    Vector3 xeq = rhs.convertToICRS();
    this->convertFromICRS(xeq);
  }
  void
  SphericalCoords::convertFrom(const SphericalCoords& rhs,
			       Matrix33& partials) {
    Matrix33 partialRHS;
    Vector3 xeq = rhs.convertToICRS(&partialRHS);
    Matrix33 partialLHS;
    this->convertFromICRS(xeq, &partialLHS);
    partials = partialLHS * partialRHS;
  }
  void
  SphericalCoords::convertFrom(const SphericalCoords& rhs,
			       Matrix32& partials) {
    Matrix33 partial3d;
    convertFrom(rhs,partial3d);
    partials = partial3d * rhs.derivsFrom2d();
  }
  void
  SphericalCoords::convertFrom(const SphericalCoords& rhs,
			       Matrix23& partials) {
    Matrix33 partial3d;
    convertFrom(rhs,partial3d);
    partials = this->derivsTo2d()*partial3d;
  }
  void
  SphericalCoords::convertFrom(const SphericalCoords& rhs,
			       Matrix22& partials) {
    Matrix23 partial23;
    convertFrom(rhs,partial23);
    partials = partial23 * rhs.derivsFrom2d();
  }

  Matrix22
  SphericalCoords::convertWithCovariance(const SphericalCoords& rhs,
					 const Matrix22& cov) {
    Matrix22 dThisdRHS;
    convertFrom(rhs, dThisdRHS);
    return dThisdRHS * cov * dThisdRHS.transpose();
  }
  void
  SphericalICRS::convertFromICRS(const Vector3& xeq, 
				     Matrix33* partials) {
    setUnitVector(xeq);
    if (partials) {
      partials->setToIdentity();
    }
  }
  Vector3
  SphericalICRS::convertToICRS(Matrix33* partials) const {
    Vector3 xeq = getUnitVector();
    if (partials) {
      partials->setToIdentity();
    }
    return xeq;
  }

  void
  SphericalEcliptic::convertFromICRS(const Vector3& xeq, 
				   Matrix33* partials) {
    setup(); // Check rotation matrix
    setUnitVector( (*icrs2ecliptic) * xeq);
    if (partials) *partials = *icrs2ecliptic;
  }

  Vector3
  SphericalEcliptic::convertToICRS(Matrix33* partials) const {
    setup();
    if (partials) *partials = icrs2ecliptic->transpose();
    return icrs2ecliptic->transpose() * getUnitVector();
  }

  void
  SphericalInvariable::convertFromICRS(const Vector3& xeq, 
				   Matrix33* partials) {
    setup(); // Check rotation matrix
    setUnitVector( (*icrs2invariable) * xeq);
    if (partials) *partials = *icrs2invariable;
  }
  
  Vector3
  SphericalInvariable::convertToICRS(Matrix33* partials) const {
    setup();
    if (partials) *partials = icrs2invariable->transpose();
    return icrs2invariable->transpose() * getUnitVector();
  }

  void
  Orientation::buildMatrix() {
    Vector3 xyz=pole.getUnitVector();
    double x=xyz[0];
    double y=xyz[1];
    double z=xyz[2];
    double rho=hypot(y,x);
    double xrho,yrho;
    if (rho<=0.) {
      // For pole at ICRS pole, define approach from -y axis
      // so phi=0 puts new x axis at old x axis.
      xrho=0.;
      yrho=-1.;
    } else {
      xrho = x/rho;
      yrho = y/rho;
    }
    double cphi=cos(zrot);
    double sphi=sin(zrot);
    rMatrix(0,0) = -xrho*sphi*z - yrho*cphi;
    rMatrix(0,1) = +xrho*cphi   - yrho*sphi*z;
    rMatrix(0,2) = rho*sphi;
    rMatrix(1,0) = -xrho*cphi*z + yrho*sphi;
    rMatrix(1,1) = -xrho*sphi   - yrho*cphi*z;
    rMatrix(1,2) = rho*cphi;
    rMatrix(2,0) = x;
    rMatrix(2,1) = y;
    rMatrix(2,2) = z;
  }

  Vector3 
  Orientation::fromICRS(const Vector3& x) const {
    return rMatrix*x;
  }

  Vector3
  Orientation::toICRS(const Vector3& x) const {
    // Use fact that rMatrix inverse is its transpose:
    return rMatrix.transpose()*x;
  }

  Orientation::Orientation(const SphericalCoords& pole_, 
			   double pa): zrot(-pa) {
    pole.convertFrom(pole_);
    buildMatrix();
  }
  void
  Orientation::setPole(const SphericalCoords& pole_) {
    pole.convertFrom(pole_);
    buildMatrix();
  }
  void
  Orientation::set(const SphericalCoords& pole_, 
		   double pa) {
    zrot=-pa;
    pole.convertFrom(pole_);
    buildMatrix();
  }
  void
  Orientation::setZRot(double zrot_) {
    zrot=zrot_;
    buildMatrix();
  }

  // Normal N->E position angle is backwards in the
  // E=x N=y system
  void 
  Orientation::alignToICRS(double phi) {
    setZRot(-phi);
  }

  // Give the angle between ICRS North vector at x0 and the
  // vector from x0 to some inclined pole.
  double
  NorthAngle(SphericalICRS xx,
	     double inclination,
	     double ascendingNode) {
    double xp=sin(inclination)*sin(ascendingNode);
    double yp=-sin(inclination)*cos(ascendingNode);
    double zp=cos(inclination);
    Vector3 xyz=xx.getUnitVector();
    double x0=xyz[0];
    double y0=xyz[1];
    double z0=xyz[2];

    // Cross product x0 x p:
    double xpx = y0*zp - z0*yp;
    double xpy = z0*xp - x0*zp;
    double xpz = x0*yp - y0*xp;
    double xpr = xpx*xpx + xpy*xpy + xpz*xpz;
    if (xpr<=0.) {
      // xx is at the new pole, punt:
      return 0.;
    }
    xpr = sqrt(xpr);

    // Cross product x0 x z
    double xzx = y0;
    double xzy = -x0;
    double xzr = hypot(xzx, xzy);
    if (xzr<=0.) {
      // xx is on z axis, punt:
      return 0.;
    }
   
    double theta = acos((xpx*xzx + xpy*xzy)/(xpr*xzr));
    // resolve sign ambiguity:
    theta = abs(theta);
    if (xp*y0-x0*yp < 0.) theta*=-1.;
    return theta;
  }

  void 
  Orientation::alignToEcliptic(double phi) {
    setZRot(NorthAngle(pole,
		       EclipticInclination,
		       EclipticNode)
	    - phi);
  }
  void 
  Orientation::alignToInvariable(double phi) {
    setZRot(NorthAngle(pole,
		       InvariableInclination,
		       InvariableNode)
	    -phi);
  }

  void
  Orientation::write(std::ostream& os) const {
    os << pole << " "
       << fixed << setprecision(8) << getPA()/DEGREE;
  }

  void
  Orientation::read(std::istream& is) {
    string paString;
    is >> pole >> paString;
    if (!is) return;
    if (nocaseEqual(paString, "Ecliptic")) {
      alignToEcliptic();
    } else {
      std::istringstream iss(paString);
      double PA;
      if (!(iss >> PA)) 
	throw AstrometryError("Bad PA for Orientation: " + paString);
      alignToICRS(PA*DEGREE);
    }
  }

  void
  ReferenceFrame::write(std::ostream& os) const {
    os << orient << endl;
    os << fixed << std::setprecision(8)
       << origin[0] << " " << origin[1] << " " << origin[2];
  }
  void
  ReferenceFrame::read(std::istream& is) {
    // The information should be on consecutive lines
    is >> orient >> origin[0] >> origin[1] >> origin[2];
  }
  std::ostream& operator<<(std::ostream& os, const ReferenceFrame& rhs) {
    rhs.write(os);
    return os;
  }
  std::istream& operator>>(std::istream& is, ReferenceFrame& rhs) {
    rhs.read(is);
    return is;
  }


  ///////////////////////////////////////////////////////////////////////
  // Base class of systems with arbitrary pole:
  ///////////////////////////////////////////////////////////////////////

  void
  SphericalCustomBase::convertFromICRS(const Vector3& xeq, 
				       Matrix33* partials) {
    setUnitVector( orient->fromICRS(xeq));
    if (partials) *partials=orient->m();
  }
  Vector3
  SphericalCustomBase::convertToICRS(Matrix33* partials) const {
    if (partials) *partials=orient->m().transpose();
    return( orient->toICRS(getUnitVector()));
  }
  

  SphericalCustomBase::SphericalCustomBase(const Orientation& o,
					   bool shareOrient): orient(&o), ownOrient(false) {
    if (!shareOrient) {
      orient = new Orientation(o);
      ownOrient = true;
    } 
  }

  SphericalCustomBase::SphericalCustomBase(double xi, double eta, 
					   const Orientation& o,
					   bool shareOrient): orient(&o), ownOrient(false) {
    if (!shareOrient) {
      orient = new Orientation(o);
      ownOrient = true;
    } 
    setLonLat(xi,eta);
  }

  SphericalCustomBase::SphericalCustomBase(const SphericalCustomBase& rhs): 
    SphericalCoords(rhs),
    orient(rhs.orient),
    ownOrient(false) {
    if (rhs.ownOrient) {
      orient = new Orientation(*rhs.orient);
      ownOrient = true;
    } 
  }

  SphericalCustomBase& 
  SphericalCustomBase::operator=(const SphericalCustomBase& rhs) {
    if (this == &rhs) return *this; // Self-assignment
    x = rhs.x;
    if (ownOrient) delete orient;
    if (rhs.ownOrient) {
      orient = new Orientation(*rhs.orient);
      ownOrient = true;
    } else {
      orient = rhs.orient;
      ownOrient = false;
    }
    return *this;
  }

  SphericalCustomBase::SphericalCustomBase(const SphericalCoords& rhs,
					   const Orientation& o,
					   bool shareOrient): orient(&o), ownOrient(false) {
    if (!shareOrient) {
      orient = new Orientation(o);
      ownOrient = true;
    }
    convertFrom(rhs);
  }
  SphericalCustomBase::SphericalCustomBase(const SphericalCoords& rhs, 
					   const Orientation& o,
					   Matrix33& partials,
					   bool shareOrient): orient(&o), ownOrient(false) {
    if (!shareOrient) {
      orient = new Orientation(o);
      ownOrient = true;
    }
    convertFrom(rhs, partials);
  }

  SphericalCustomBase::SphericalCustomBase(const CartesianCustom& rhs):
    orient(&(rhs.getFrame()->orient)),
    ownOrient(false) {
    if (rhs.ownFrame) {
      orient = new Orientation(rhs.getFrame()->orient);
      ownOrient = true;
    }
    projectIt(rhs);
  }
  SphericalCustomBase::SphericalCustomBase(const CartesianCustom& rhs,
					   Matrix33& partials):
    orient(&(rhs.getFrame()->orient)),
    ownOrient(false) {
    if (rhs.ownFrame) {
      orient = new Orientation(rhs.getFrame()->orient);
      ownOrient = true;
    }
    projectIt(rhs, &partials);
  }

  void
  Gnomonic::setLonLat(double xi, double eta) {
    x[2] = pow(1.+xi*xi+eta*eta , -0.5);
    //if (x[2] < 0.5){
    //   cerr << xi << "  " << eta << "  " << x[2] << endl;
    //}
    x[0] = xi*x[2];
    x[1] = eta*x[2];
  }
  void
  Gnomonic::getLonLat(double& xi, double& eta) const {
    if (x[2]<=0.){ 
      //cerr << x[0] << "  " << x[1] << "  " << x[2] << endl;
      throw AstrometryError("Gnomonic projection invalid for z<=0.");
    }
    xi = x[0]/x[2];
    eta= x[1]/x[2];
  }

  Matrix23
  Gnomonic::derivsTo2d() const {
    Matrix23 project(0.);
    if (x[2]<=0.)
      throw AstrometryError("Gnomonic projection derivs invalid for z<=0.");
    double invz=1./x[2];
    project(0,0) = invz;
    project(0,2) = -x[0]*invz*invz;
    project(1,1) = invz;
    project(1,2) = -x[1]*invz*invz;
    return project;
  }
  Matrix32
  Gnomonic::derivsFrom2d() const {
    Matrix32 project;
    project(0,0) = 1 - x[0]*x[0];
    project(0,1) =   - x[0]*x[1];
    project(1,0) =   - x[0]*x[1];
    project(1,1) = 1 - x[1]*x[1];
    project(2,0) =   - x[0]*x[2];
    project(2,1) =   - x[1]*x[2];
    return project;
  }

  // Tangent-plane I/O will use degrees, degrees
  void
  Gnomonic::write(std::ostream& os, int decimalPlaces) const {
    double lon, lat;
    getLonLat(lon,lat);
    string buffer = degdms(lon/DEGREE, decimalPlaces);
    os << buffer;
    buffer = degdms(lat/DEGREE,decimalPlaces);
    os << " " << buffer;
  }
  void
  Gnomonic::read(std::istream& is) {
    double lon, lat;
    string buffer;
    double deg;
    x.setZero();
    if (!(is >> buffer)) return;
    try {
      deg=dmsdeg(buffer);
    } catch (AstrometryError& e) {
      // Format error in the hms read
      is.setstate(is.rdstate() | std::ios_base::failbit);
      return;
    }
    lon = deg*DEGREE;
    if (!(is >> buffer)) return;
    try {
      deg=dmsdeg(buffer);
    } catch (AstrometryError& e) {
      // Format error in the hms read
      is.setstate(is.rdstate() | std::ios_base::failbit);
      return;
    }
    lat = deg*DEGREE;
    setLonLat(lon,lat);
  }

  ///////////////////////////////////////////////////////////////////////

  // Projection from 3d down to 2d.  Dimension of partials says whether
  // desired derivs are to the unit vectors or to the 2d representation.
  void
  SphericalCoords::projectIt(const CartesianCoords& rhs,
			     Matrix33* partials) {
    double r=rhs.radius();
    if (r<=0.) 
      throw AstrometryError("Projecting null 3d coordinates to spherical");
    x = rhs.getVector();
    x *= 1./r;
    if (!partials) return;
    (*partials)(0,0) =( 1. - x[0]*x[0])/r;
    (*partials)(0,1) =(    - x[0]*x[1])/r;
    (*partials)(0,2) =(    - x[0]*x[2])/r;
    (*partials)(1,0) =(    - x[1]*x[0])/r;
    (*partials)(1,1) =( 1. - x[1]*x[1])/r;
    (*partials)(1,2) =(    - x[1]*x[2])/r;
    (*partials)(2,0) =(    - x[2]*x[0])/r;
    (*partials)(2,1) =(    - x[2]*x[1])/r;
    (*partials)(2,2) =( 1. - x[2]*x[2])/r;
  }
  void
  SphericalCoords::projectIt(const CartesianCoords& rhs,
			     Matrix23* partials) {
    if (!partials) projectIt(rhs, (Matrix23*) 0);
    // We want derivatives of 2d coords w.r.t 3d
    Matrix33 partial3d;
    projectIt(rhs, &partial3d);
    *partials = derivsTo2d() * partial3d;
  }

  // Here is the list of spherical coordinate transformation
  // constructors and projections from their related Cartesian
  // systems.

  SphericalICRS::SphericalICRS(const SphericalCoords& rhs) {
    convertFrom(rhs);
  }
  SphericalICRS::SphericalICRS(const SphericalCoords& rhs,
			       Matrix33& partials) {
    convertFrom(rhs, partials);
  }
  SphericalICRS::SphericalICRS(const CartesianICRS& rhs) {
    projectIt(rhs);
  }
  SphericalICRS::SphericalICRS(const CartesianICRS& rhs,
			       Matrix33& partials) {
    projectIt(rhs, &partials);
  }
  SphericalICRS::SphericalICRS(const CartesianICRS& rhs,
			       Matrix23& partials) {
    projectIt(rhs, &partials);
  }

  SphericalEcliptic::SphericalEcliptic(const SphericalCoords& rhs) {
    convertFrom(rhs);
  }
  SphericalEcliptic::SphericalEcliptic(const SphericalCoords& rhs,
			       Matrix33& partials) {
    convertFrom(rhs, partials);
  }
  SphericalEcliptic::SphericalEcliptic(const CartesianEcliptic& rhs) {
    projectIt(rhs);
  }
  SphericalEcliptic::SphericalEcliptic(const CartesianEcliptic& rhs,
			       Matrix33& partials) {
    projectIt(rhs, &partials);
  }

  SphericalInvariable::SphericalInvariable(const SphericalCoords& rhs) {
    convertFrom(rhs);
  }
  SphericalInvariable::SphericalInvariable(const SphericalCoords& rhs,
			       Matrix33& partials) {
    convertFrom(rhs, partials);
  }
  SphericalInvariable::SphericalInvariable(const CartesianInvariable& rhs) {
    projectIt(rhs);
  }
  SphericalInvariable::SphericalInvariable(const CartesianInvariable& rhs,
			       Matrix33& partials) {
    projectIt(rhs, &partials);
  }


#ifdef USE_YAML
  // YAML (de-)serialization routines.
  // First these are the structure that holds the names and creators for
  // all possible types of Spherical coordinates, and a routine
  // that will set up a standard set and be called when any use of
  // deserialization is made.

  std::map<string, SphericalCoords::Creator>
  SphericalCoords::typemap;
  
  void
  SphericalCoords::initializeTypemap() {
    if (typemap.empty()) {
      // Load the known SphericalCoordinate classes
      typemap[SphericalICRS::type()]       = SphericalICRS::create;
      typemap[SphericalEcliptic::type()]   = SphericalEcliptic::create;
      typemap[SphericalInvariable::type()] = SphericalInvariable::create;
      typemap[SphericalCustom::type()]     = SphericalCustom::create;
      typemap[Gnomonic::type()]            = Gnomonic::create;
    }
  }

  // And the factory itself:
  SphericalCoords*
  SphericalCoords::deserialize(const YAML::Node& node) {
    initializeTypemap();
    if (!node.IsMap() || !node["Type"])
      throw AstrometryError("Missing YAML key <Type> in SphericalCoords::deserialize");
    string t = node["Type"].as<string>();
    auto it = typemap.find(t);
    if (it==typemap.end())
      throw AstrometryError("Unknown SphericalCoord type <" + t + "> in YAML");
    // Call the type's factory on the node data:
    return (*(it->second))(node);
  }

  // (De-)serialization routines for all declared SphericalCoord types
  void
  SphericalICRS::serialize(YAML::Emitter& os) const {
    double lon, lat;
    getLonLat(lon, lat);
    while (lon<0.) lon+= TPI;
    os << YAML::BeginMap
       << YAML::Key << "Type" << YAML::Value << type()
       << YAML::Key << "RA"   << YAML::Value << lon/DEGREE
       << YAML::Key << "Dec"  << YAML::Value << lat/DEGREE
       << YAML::EndMap;
  }
  SphericalCoords*
  SphericalICRS::create(const YAML::Node& node) {
    if (!node.IsMap() || node["Type"].as<string>()!=type())
      throw AstrometryError("Missing/bad YAML key <Type> for " + type());
    double lon=0., lat=0.;   // Defaults to vernal equinox
    if (node["RA"])
      lon = node["RA"].as<double>() * DEGREE;
    if (node["Dec"])
      lat = node["Dec"].as<double>() * DEGREE;
    return new SphericalICRS(lon,lat);
  }

  void
  SphericalEcliptic::serialize(YAML::Emitter& os) const {
    double lon,lat;
    getLonLat(lon, lat);
    while (lon<0.) lon+= TPI;
    os << YAML::BeginMap
       << YAML::Key << "Type" << YAML::Value << type()
       << YAML::Key << "Lon"   << YAML::Value << lon/DEGREE
       << YAML::Key << "Lat"  << YAML::Value << lat/DEGREE
       << YAML::EndMap;
  }
  SphericalCoords*
  SphericalEcliptic::create(const YAML::Node& node) {
    if (!node.IsMap() || node["Type"].as<string>()!=type())
      throw AstrometryError("Missing/bad YAML key <Type> for " + type());
    double lon=0., lat=0.;   // Defaults to vernal equinox
    if (node["Lon"])
      lon = node["Lon"].as<double>() * DEGREE;
    if (node["Lat"])
      lat = node["Lat"].as<double>() * DEGREE;
    return new SphericalEcliptic(lon,lat);
  }

  void
  SphericalInvariable::serialize(YAML::Emitter& os) const {
    double lon,lat;
    getLonLat(lon, lat);
    while (lon<0.) lon+= TPI;
    os << YAML::BeginMap
       << YAML::Key << "Type" << YAML::Value << type()
       << YAML::Key << "Lon"   << YAML::Value << lon/DEGREE
       << YAML::Key << "Lat"  << YAML::Value << lat/DEGREE
       << YAML::EndMap;
  }
  SphericalCoords*
  SphericalInvariable::create(const YAML::Node& node) {
    if (!node.IsMap() || node["Type"].as<string>()!=type())
      throw AstrometryError("Missing/bad YAML key <Type> for " + type());
    double lon=0., lat=0.;   // Defaults to vernal equinox
    if (node["Lon"])
      lon = node["Lon"].as<double>() * DEGREE;
    if (node["Lat"])
      lat = node["Lat"].as<double>() * DEGREE;
    return new SphericalInvariable(lon,lat);
  }
  
  void
  Orientation::serialize(YAML::Emitter& os) const {
    double ra, dec;
    pole.getRADec(ra,dec);
    while (ra<0.) ra+= TPI;
    os << YAML::BeginMap
       << YAML::Key << "RA"   << YAML::Value << ra/DEGREE
       << YAML::Key << "Dec"  << YAML::Value << dec/DEGREE
       << YAML::Key << "PA"  << YAML::Value << getPA()/DEGREE
       << YAML::EndMap;
  }
  Orientation*
  Orientation::create(const YAML::Node& node) {
    // Note there is no defaulting allowed here
    if (!node.IsMap() || !node["RA"] || !node["Dec"] || !node["PA"])
      throw AstrometryError("Missing/bad YAML keys for Orientation");
    double ra = stod(node["RA"].as<string>()) * DEGREE;
    double dec = stod(node["Dec"].as<string>()) * DEGREE;
    double pa = stod(node["PA"].as<string>()) * DEGREE;
    return new Orientation(SphericalICRS(ra,dec), pa);
  }
  
  void
  SphericalCustom::serialize(YAML::Emitter& os) const {
    double lon,lat;
    getLonLat(lon, lat);
    while (lon<0.) lon+= TPI;
    os << YAML::BeginMap
       << YAML::Key << "Type" << YAML::Value << type()
       << YAML::Key << "Orientation"   << YAML::Value << *getOrient()
       << YAML::Key << "Lon"   << YAML::Value << lon/DEGREE
       << YAML::Key << "Lat"  << YAML::Value << lat/DEGREE
       << YAML::EndMap;
  }
  SphericalCoords*
  SphericalCustom::create(const YAML::Node& node) {
    if (!node.IsMap() || node["Type"].as<string>()!=type())
      throw AstrometryError("Missing/bad YAML key <Type> for " + type());
    if (!node["Orientation"])
      throw AstrometryError("Missing/bad YAML key <Orientation> for " + type());
    Orientation* orient = Orientation::create(node["Orientation"]);
    double lon=0., lat=0.;   // Defaults to "equinox" position
    if (node["Lon"])
      lon = node["Lon"].as<double>() * DEGREE;
    if (node["Lat"])
      lat = node["Lat"].as<double>() * DEGREE;
    auto out = new SphericalCustom(lon, lat, *orient);
    delete orient;
    return out;
  }

  void
  Gnomonic::serialize(YAML::Emitter& os) const {
    double xi,eta;
    try {
      getLonLat(xi, eta);
    } catch (AstrometryError& m) {
      // Uninitialized coords might throw projection error here.  Just get 0's
      xi = 0.;
      eta = 0.;
    }
    os << YAML::BeginMap
       << YAML::Key << "Type" << YAML::Value << type()
       << YAML::Key << "Orientation"   << YAML::Value << *getOrient()
       << YAML::Key << "Xi"   << YAML::Value << xi/DEGREE
       << YAML::Key << "Eta"  << YAML::Value << eta/DEGREE
       << YAML::EndMap;
  }
  SphericalCoords*
  Gnomonic::create(const YAML::Node& node) {
    if (!node.IsMap() || node["Type"].as<string>()!=type())
      throw AstrometryError("Missing/bad YAML key <Type> for " + type());
    Orientation* orient=0;
    if (!node["Orientation"]) {
      // For backwards compatibility, look for RA/Dec/PA
      if (!node["RA"] || !node["Dec"] || !node["PA"])
	throw AstrometryError("Missing/bad YAML key <Orientation> for " + type());
      double ra = stod(node["RA"].as<string>());
      double dec = stod(node["Dec"].as<string>());
      double pa = stod(node["PA"].as<string>());
      orient = new Orientation(SphericalICRS(ra*DEGREE,dec*DEGREE),pa*DEGREE);
    } else {
      orient = Orientation::create(node["Orientation"]);
    }
    double xi=0., eta=0.;   // Defaults to pole position
    if (node["Xi"])
      xi = stod(node["Xi"].as<string>()) * DEGREE;
    if (node["Eta"])
      eta = stod(node["Eta"].as<string>()) * DEGREE;
    auto out = new Gnomonic(xi, eta, *orient);
    delete orient;
    return out;
  }

#endif

} // namespace astrometry
