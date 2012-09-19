/* 	$Id: JPL.h,v 1.1.1.1 2009/11/02 03:43:48 garyb Exp $	 */
// Definitions used by Hoffman's DE405 C code.
// An alteration of ephem_types.h to reflect embedding in C++.

#ifndef JPL_H
#define JPL_H

/******************************************************************************/
/**                                                                          **/
/**  SOURCE FILE: ephem_types.h                                              **/
/**                                                                          **/
/**    This file contains C macros and type definitions for this version of  **/
/**    the JPL ephemeris utilities. Note that the macro EPHEMERIS, defined   **/
/**    below, must be consistent with the ephemeris being read in order for  **/
/**    any of the C utilities to work properly. Note also that the addition  **/
/**    of new ephemerides to the FTP site might require new values for the   **/
/**    macro ARRAY_SIZE.                                                     **/
/**                                                                          **/
/**   Programmer: David Hoffman/EG5                                          **/
/**               NASA, Johnson Space Center                                 **/
/**               Houston, TX 77058                                          **/
/**               e-mail: david.a.hoffman1@jsc.nasa.gov                      **/
/**                                                                          **/
/******************************************************************************/

/**==========================================================================**/
/**  C Macro Definitions                                                     **/
/**                                                                          **/
/**    The name of the ephemeris must be defined by the user. At the time    **/
/**    this software was created, the only  ephemerides available to the     **/
/**    public were DE200 and DE405. The ephemeris name is used to size the   **/ 
/**    records in the binary data file. These records are sized to accom-    **/
/**    odate a complete set of  interpolating coefficients in the form of    **/
/**    an array of double precision numbers. The parameter ARRAY_SIZE gives  **/
/**    the number of such coefficients in each record; it is used to size    **/
/**    the header records defined below. Finally, some C pre-processor       **/
/*     macros are defined for convenience.                                   **/
/**==========================================================================**/

namespace ephem {
#define EPHEMERIS 405                /* Note the obvious: input XXX for DEXXX */

#if EPHEMERIS==200  
  const int ARRAY_SIZE=826;
#elif EPHEMERIS==405
  const int  ARRAY_SIZE=1018;
#endif

/**==========================================================================**/
/**  C Record Type Definitions                                               **/
/**==========================================================================**/

   /*-------------------------------------------------------------------------*/
   /* Define the content of binary header records                             */
   /*-------------------------------------------------------------------------*/

   struct RecOneType {
         char label[3][84];
         char constName[400][6];
       double timeData[3];
     long int numConst;
       double AU;
       double EMRAT;
     long int coeffPtr[12][3];
     long int DENUM;
     long int libratPtr[3];
   };

   struct RecTwoType {
     double constValue[400];
   }; 

   /*-------------------------------------------------------------------------*/
   /* Define binary header record formats                                     */
   /*-------------------------------------------------------------------------*/

   struct HeaderOneType {
       RecOneType data;
       char pad[ARRAY_SIZE*sizeof(double) - sizeof(RecOneType)];
   };
  
   struct HeaderTwoType {
     RecTwoType data;
     char pad[ARRAY_SIZE*sizeof(double) - sizeof(RecTwoType)];
   };

}  // namespace ephem
#endif  // JPL_H

