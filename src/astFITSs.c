#include <stdio.h>
#include <string.h>
#include "astsys.h"
#include "astLib.h"
int astFITSs ( struct WCS wcsij,
               FRAMETYPE frame, struct EPOCH eqx, double djm,
               char* ctype1, char* crpix1, char* crval1,
               char* ctype2, char* crpix2, char* crval2,
               char* cd1_1, char* cd1_2,
               char* cd2_1, char* cd2_2,
               char* radecsys, char* equinox, char* mjdobs )
/*
**  - - - - - - - - -
**   a s t F I T S s
**  - - - - - - - - -
**
**  Express a world coordinate system transformation as FITS header
**  strings.
**
**  GIVEN:
**    wcsij     struct     WCS transformation, pixels to sky
**    frame     FRAMETYPE  type of sky coordinate system
**    eqx       struct     equinox (mean RA/Decs only)
**    djm       double     epoch of observation (MJD: JD-2400000.5)
**
**  RETURNED (arguments):
**    ctype1    char[81]   FITS header record:  CTYPE1
**    crpix1    char[81]   FITS header record:  CRPIX1
**    crval1    char[81]   FITS header record:  CRVAL1
**    ctype2    char[81]   FITS header record:  CTYPE1
**    crpix2    char[81]   FITS header record:  CRPIX2
**    crval2    char[81]   FITS header record:  CRVAL2
**    cd1_1     char[81]   FITS header record:  CD1_1
**    cd1_2     char[81]   FITS header record:  CD1_2
**    cd2_1     char[81]   FITS header record:  CD2_1
**    cd2_2     char[81]   FITS header record:  CD2_2
**    radecsys  char[81]   FITS header record:  RADECSYS
**    equinox   char[81]   FITS header record:  EQUINOX
**    mjdobs    char[81]   FITS header record:  MJD-OBS
**
**  RETURNED (function value):
**              int        0 = OK
**                        -1 = illegal frame
**
**  Defined in astLib.h:
**    WCS       struct     WCS transformation parameters
**    FRAMETYPE enum       coordinate system IDs
**    EPOCH     struct     B or J epoch
**
**  To guide potential non-Gemini-TCS users of this routine, here are
**  the present definitions of the above.
**
**  struct WCS {
**        double ab0[2];      / * Celestial coordinates at x=y=0 * /
**        double coeffs[6];   / * Affine transformation coefficients * /
**  };
**
**  typedef enum { AZEL_MNT = 0,    / * Mount Az/El, pre-flexure * /
**                 AZEL_TOPO = 1,   / * Topocentric Az/El * /
**                 APPT = 2,        / * Geocentric apparent RA/Dec * /
**                 FK5 = 3,         / * IAU 1976 RA/Dec, any equinox * /
**                 FK4 = 4          / * Pre IAU 1976 RA/Dec, any equinox * /
**               } FRAMETYPE;
**
**  struct EPOCH {
**     double year;     / * Epoch:  Byear (B), TT Jyear (J) * /
**     char   type;     / * Type of epoch ('B', 'J' or ' ') * /
**  };
**
**  Called:  astFITSv
**
**  Notes:
**
**  1  The output strings are exactly 81 characters long.  They have
**     trailing nulls but no trailing spaces.  (If used to write a
**     FITS file, trailing spaces must be appended.)
**
**  2  The keywords EQUINOX and/or RADECSYS may not be required.  An
**     empty string is returned in such cases.
**
**  3  The timescale for the epoch of observation, djm, is not defined
**     by the current FITS standard.  TAI is a suitable choice;  UTC
**     is not, strictly speaking, because expressing UTC as MJD leads
**     to ambiguities during leap seconds.
**
**  4  The present function returns strings.  Another function, astFITSv,
**     returns values.
**
**  5  Here is an example of the set of 13 strings returned by this routine:
**
**     CTYPE1  = 'RA---TAN'           / gnomonic projection
**     CRPIX1  =           984.048024 / pixel i-coordinate at rotator axis
**     CRVAL1  =       36.00179953283 / RA or Az at rotator axis
**     CTYPE2  = 'DEC--TAN'           / gnomonic projection
**     CRPIX2  =           438.647662 / pixel j-coordinate at rotator axis
**     CRVAL2  =       44.99795763567 / Dec or El at rotator axis
**     CD1_1   = -0.00001007236440932 / xi rotation/skew/scale matrix element
**     CD1_2   =  0.00000000574196532 / xj rotation/skew/scale matrix element
**     CD2_1   = -0.00000000499131497 / yi rotation/skew/scale matrix element
**     CD2_2   =  0.00001012093669430 / yj rotation/skew/scale matrix element
**     RADECSYS= 'FK5     '           / type of RA/Dec
**     EQUINOX =            2000.0000 / epoch of mean equator & equinox
**     MJD-OBS =     49560.6433912037 / epoch of observation (TAI MJD)
**
**  P.T.Wallace   22 November 1998
**
**  Copyright RAL 1998.  All rights reserved.
*/

{
   char ctype1_v[9], ctype2_v[9], radecsys_v[9];
   double crpix1_v, crval1_v, crpix2_v, crval2_v,
          cd1_1_v, cd1_2_v, cd2_1_v, cd2_2_v,
          equinox_v, mjdobs_v;
   int j;


/* Get the values. */
   if ( (j = astFITSv ( wcsij, frame, eqx, djm,
                       ctype1_v, &crpix1_v, &crval1_v,
                       ctype2_v, &crpix2_v, &crval2_v,
                       &cd1_1_v, &cd1_2_v, &cd2_1_v, &cd2_2_v,
                       radecsys_v, &equinox_v, &mjdobs_v )) ) return j;

/* Build the strings. */
   sprintf ( ctype1,
             "CTYPE1  = '%8s'           / gnomonic projection",
             ctype1_v );
   sprintf ( crpix1,
             "CRPIX1  = %20.6f / pixel i-coordinate at rotator axis",
             crpix1_v );
   sprintf ( crval1,
             "CRVAL1  = %20.11f / RA or Az at rotator axis",
             crval1_v );
   sprintf ( ctype2,
             "CTYPE2  = '%8s'           / gnomonic projection",
             ctype2_v );
   sprintf ( crpix2,
             "CRPIX2  = %20.6f / pixel j-coordinate at rotator axis", 
             crpix2_v );
   sprintf ( crval2,
             "CRVAL2  = %20.11f / Dec or El at rotator axis",
             crval2_v );
   sprintf ( cd1_1,
             "CD1_1   = %20.17f / xi rotation/skew/scale matrix element",
             cd1_1_v );
   sprintf ( cd1_2,
             "CD1_2   = %20.17f / xj rotation/skew/scale matrix element",
             cd1_2_v );
   sprintf ( cd2_1,
             "CD2_1   = %20.17f / yi rotation/skew/scale matrix element",
             cd2_1_v );
   sprintf ( cd2_2,
             "CD2_2   = %20.17f / yj rotation/skew/scale matrix element",
             cd2_2_v );
   sprintf ( radecsys,
             "RADECSYS= '%8s'           / type of RA/Dec",
             radecsys_v );
   sprintf ( equinox,
             "EQUINOX = %20.4f / epoch of mean equator & equinox",
             equinox_v );
   sprintf ( mjdobs,
             "MJD-OBS = %20.10f / epoch of observation (TAI MJD)",
             mjdobs_v );

/* Empty any records which are superfluous for this case. */
   if ( !*radecsys_v ) radecsys[0] = (char) '\0';
   if ( equinox_v < 0.0 ) equinox[0] = (char) '\0';

   return 0;
}
