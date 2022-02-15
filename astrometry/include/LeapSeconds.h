// This function is declared in Astrometry.h but I have separated
// its definition into this file because the list will have to be
// maintained as new leap seconds occur.

// Data for calculation of leap seconds and TT-UTC.
// These tabular data were obtained from 
//     ftp://maia.usno.navy.mil/ser7/tai-utc.dat
// Add leap seconds at TOP of table as they occur.
#ifndef LEAPSECOND_H
#define LEAPSECOND_H

#include <vector>
#include "AstronomicalConstants.h"
namespace astrometry {
  double UT::TAIminusUTC(double jd_) {
    static std::vector<double> leapJD;
    static std::vector<double> dT;
    static bool initialized=false;
    if (!initialized) {
      leapJD.push_back(2457754.5); dT.push_back(37.0);	//Jan 1  2017
      leapJD.push_back(2457204.5); dT.push_back(36.0);	//Jun 30 2015
      leapJD.push_back(2456109.5); dT.push_back(35.0);	//Jun 30 2012
      leapJD.push_back(2454832.5); dT.push_back(34.0);	//Jan 1 2009
      leapJD.push_back(2453736.5); dT.push_back(33.0);	//Jan 1 2006
      leapJD.push_back(2451179.5); dT.push_back(32.0);	//Jan 1 1999
      leapJD.push_back(2450630.5); dT.push_back(31.0);
      leapJD.push_back(2450083.5); dT.push_back(30.0);
      leapJD.push_back(2449534.5); dT.push_back(29.0);
      leapJD.push_back(2449169.5); dT.push_back(28.0);
      leapJD.push_back(2448804.5); dT.push_back(27.0);
      leapJD.push_back(2448257.5); dT.push_back(26.0);
      leapJD.push_back(2447892.5); dT.push_back(25.0);
      leapJD.push_back(2447161.5); dT.push_back(24.0);
      leapJD.push_back(2446247.5); dT.push_back(23.0);
      leapJD.push_back(2445516.5); dT.push_back(22.0);
      leapJD.push_back(2445151.5); dT.push_back(21.0);
      leapJD.push_back(2444786.5); dT.push_back(20.0);
      leapJD.push_back(2444239.5); dT.push_back(19.0);
      leapJD.push_back(2443874.5); dT.push_back(18.0);
      leapJD.push_back(2443509.5); dT.push_back(17.0);
      leapJD.push_back(2443144.5); dT.push_back(16.0);
      leapJD.push_back(2442778.5); dT.push_back(15.0);
      leapJD.push_back(2442413.5); dT.push_back(14.0);
      leapJD.push_back(2442048.5); dT.push_back(13.0);
      leapJD.push_back(2441683.5); dT.push_back(12.0);
      leapJD.push_back(2441499.5); dT.push_back(11.0);
      leapJD.push_back(2441317.5); dT.push_back(10.0);	//Jan 1 1972
      initialized = true;
    }
    for (int i=0; i<leapJD.size(); i++)
      if (jd_>=leapJD[i]) return dT[i]*TIMESEC;
    return dT.back()*TIMESEC;
    // Note conversion to AstronomicalConstants.h units with TIMESEC'S
  }
} //namespace astrometry
#endif
