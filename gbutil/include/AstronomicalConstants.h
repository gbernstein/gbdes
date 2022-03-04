// Include file with common "magic numbers" for Astronomy
#ifndef ASTROCONST_H
#define ASTROCONST_H

#ifndef PI
#define PI 3.14159265358979323
#endif

// conversion constants:
const double RadToArcsec=3600.*180./PI;
const double AU=1.495978707e11;	    //Astronomical Unit, m (DE430 defined)
const double Parsec=AU*RadToArcsec;
const double ABZeropoint=3.63e-23;         //AB=0 mag, in W/Hz/m^2
const double Jansky=1e-26;	           //convert Jy to W/Hz/s/m^2
const double Planck=6.626068e-34;          //Planck's const in Js.
const double SpeedOfLight=2.99792458e8;    //speed of light in m/s
const double Boltzmann=1.38065e-23;	   //Boltzmann's const, J/K
const double Micron=1e-6;		   //micron, in m
const double HubbleLengthMpc=SpeedOfLight*1e-5; //c/H0, in h^-1 Mpc
const double HubbleDistance=HubbleLengthMpc*1e6*Parsec; //c/H0, in h^-1 m
const double RecombinationRedshift=1088.;

// Orbital Mechanics constants, units of AU, solar mass, and Year
const double  TPI           = 2.*PI;
const double  DEGREE	    = PI/180.;               // One degree, in rad
const double  GM            = 4.*PI*PI/1.000037773533;  //solar gravitation
const double  ARCSEC	    = PI/180./3600.;
const double  MILLIARCSEC   = 0.001*ARCSEC;
const double  ARCMIN	    = PI/180./60.;
const double  YEAR          = 1.;                    //Julian year
const double  DAY	    = 1./365.25;	     //Julian day, 86400 s
const double  HOUR          = DAY/24.;
const double  MINUTE        = HOUR/60.;
const double  TIMESEC       = MINUTE/60.;   // Note SECOND is name conflict with MKL
const double  METER         = 1./AU;
const double  SpeedOfLightAU= SpeedOfLight/AU/TIMESEC; //in AU/YR
// Gauss's constant, DE430 defined value, = sqrt(GM_sun) to w/in errors
const double  GaussK        = 0.01720209895 / DAY;         //

//Obliquity of ecliptic at J2000
const double  EclipticInclination =23.43928*DEGREE;  
const double  EclipticNode        =0.;
//Inclination and ascending node of invariable plane in the J2000
//equatorial system: ??? check this
const double  InvariableInclination=(23+22.11/3600.)*DEGREE;
const double  InvariableNode       = (3+(52.+23.7/60.)/60.)*DEGREE;

const double MercuryGM = 6.55371264e-06;
const double VenusGM =   9.66331433e-05;
const double EarthMoonGM=1.20026937e-04;
const double MarsGM =    1.27397978e-05;
const double JupiterGM = 3.76844407e-02 + 7.80e-6; // Include satellites
const double SaturnGM =  1.12830982e-02 + 2.79e-6;
const double UranusGM =  1.72348553e-03 + 0.18e-6;
const double NeptuneGM = 2.03318556e-03 + 0.43e-6;
const double SolarSystemGM= GM + MercuryGM + VenusGM + EarthMoonGM +
  MarsGM + JupiterGM + SaturnGM + UranusGM + NeptuneGM; // Not quite everything here but ok.

const double EarthMass        = 3.00349e-6;	//Earth mass in Solar units
const double MJD0             =2400000.5;	//Offset for modified Julian dates
const double JD2000           =2451545.0;       //JD of ICRS reference epoch


#endif  // ASTROCONST_H

