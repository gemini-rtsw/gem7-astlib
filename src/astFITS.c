#include <stdio.h>
#include <string.h>
#include <slalib.h>
#include "astsys.h"
#include "astLib.h"
int astFITS ( struct WCS wcsij,
              FRAMETYPE frame, struct EPOCH eqx, double tt,
              char* ctype1, char* crpix1, char* crval1, char* cdelt1,
              char* cunit1, char* ctype2, char* crpix2, char* crval2,
              char* cdelt2, char* cunit2,
              char* pc001001, char* pc001002,
              char* pc002001, char* pc002002,
              char* radecsys, char* equinox, char* mjdobs )
/*
**  - - - - - - - -
**   a s t F I T S
**  - - - - - - - -
**
**  Express a world coordinate system transformation as FITS headers.
**
**  GIVEN:
**    wcsij     struct     WCS transformation, pixels to sky
**    frame     FRAMETYPE  type of sky coordinate system
**    eqx       struct     equinox (mean RA/Decs only)
**    tt        double     epoch (TT MJD)
**
**  RETURNED (arguments):
**    ctype1    char[80]   FITS header record:  CTYPE1
**    crpix1    char[80]   FITS header record:  CRPIX1
**    crval1    char[80]   FITS header record:  CRVAL1
**    cdelt1    char[80]   FITS header record:  CDELT1
**    cunit1    char[80]   FITS header record:  CUNIT1
**    ctype2    char[80]   FITS header record:  CTYPE1
**    crpix2    char[80]   FITS header record:  CRPIX2
**    crval2    char[80]   FITS header record:  CRVAL2
**    cdelt2    char[80]   FITS header record:  CDELT2
**    cunit2    char[80]   FITS header record:  CUNIT2
**    pc001001  char[80]   FITS header record:  PC001001
**    pc001002  char[80]   FITS header record:  PC001002
**    pc002001  char[80]   FITS header record:  PC002001
**    pc002002  char[80]   FITS header record:  PC002002
**    radecsys  char[80]   FITS header record:  RADECSYS
**    equinox   char[80]   FITS header record:  EQUINOX
**    mjdobs    char[80]   FITS header record:  MJD-OBS
**
**  RETURNED (function value):
**              int        0 = OK
**                        -1 = illegal frame
**
**  Defined in astsys.h:
**    R2D       double     radians to degrees
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
**  Called:  slaDranrm
**
**  Notes:
**
**  1  The output strings are null-terminated without trailing spaces.
**     However, the recipient arrays must contain room for the maximum
**     81 characters to allow for any subsequent spacefilling.
**
**  2  Unused keywords produce empty strings.
**
**  3  Here is an example of the set of 17 strings output by this routine:
**
**     CTYPE1  = 'RA---TAN'           / TAN projection
**     CRPIX1  =     984.044158264385 / pixel i-coordinate at rotator axis
**     CRVAL1  =      36.000410440627 / RA at rotator axis
**     CDELT1  =      -0.000010086847 / degrees per 1st-axis unit
**     CUNIT1  = 'deg     '           / radians per degree
**     CTYPE2  = 'DEC--TAN'           / TAN projection
**     CRPIX2  =     438.644955473382 / pixel j-coordinate at rotator axis
**     CRVAL2  =     +44.997610683246 / Dec at rotator axis
**     CDELT2  =      +0.000010110504 / degrees per 2nd-axis unit
**     CUNIT2  = 'deg     '           / radians per degree
**     PC001001=       0.905380195975 / xi rotation/skew matrix element
**     PC001002=      -0.424601814334 / xj rotation/skew matrix element
**     PC002001=       0.420605831504 / yi rotation/skew matrix element
**     PC002002=       0.907243481379 / yj rotation/skew matrix element
**     RADECSYS= 'FK5     '           / type of RA/Dec
**     EQUINOX =    2000.000000000000 / epoch of mean equator & equinox
**     MJD-OBS =   49560.643763703702 / epoch of observation (TT MJD)
**
**  P.T.Wallace   2 May 1997
**
**  Copyright RAL 1997.  All rights reserved.
*/

{
   double a, b, c, d, e, f, w1, w2, det;

/* Copy the WCS coefficients. */
   a = wcsij.coeffs[0];
   b = wcsij.coeffs[1];
   c = wcsij.coeffs[2];
   d = wcsij.coeffs[3];
   e = wcsij.coeffs[4];
   f = wcsij.coeffs[5];
   w1 = sqrt ( b * b + c * c );
   w2 = sqrt ( e * e + f * f );

/* Initialize all the strings. */
   strncpy ( ctype1,
"CTYPE1  =                      / TAN projection", 81 );
   strncpy ( crpix1,
"CRPIX1  =                      / pixel i-coordinate at rotator axis", 81 );
   strncpy ( crval1,
"CRVAL1  =                      / ", 81 );
   strncpy ( cdelt1,
"CDELT1  =                      / degrees per 1st-axis unit", 81 );
   strncpy ( cunit1,
"CUNIT1  = 'deg      '          / radians per degree", 81 );
   strncpy ( ctype2,
"CTYPE2  =                      / TAN projection", 81 );
   strncpy ( crpix2,
"CRPIX2  =                      / pixel j-coordinate at rotator axis", 81 );
   strncpy ( crval2,
"CRVAL2  =                      / ", 81 );
   strncpy ( cdelt2,
"CDELT2  =                      / degrees per 2nd-axis unit", 81 );
   strncpy ( cunit2,
"CUNIT2  = 'deg      '          / radians per degree", 81 );
   strncpy ( pc001001,
"PC001001=                      / xi rotation/skew matrix element", 81 );
   strncpy ( pc001002,
"PC001002=                      / xj rotation/skew matrix element", 81 );
   strncpy ( pc002001,
"PC002001=                      / yi rotation/skew matrix element", 81 );
   strncpy ( pc002002,
"PC002002=                      / yj rotation/skew matrix element", 81 );
   strncpy ( radecsys,
"RADECSYS=                      / type of RA/Dec", 81 );
   strncpy ( equinox,
"EQUINOX =                      / epoch of mean equator & equinox", 81 );
   strncpy ( mjdobs,
"MJD-OBS =                      / epoch of observation (TT MJD)", 81 );

/* Axis types. */
   switch ( frame ) {

   case AZEL_TOPO:
      strcpy ( crval1 + 33, "Az at rotator axis" );
      strncpy ( ctype1 + 10, "'AZ---TAN'", 10 );
      strcpy ( crval2 + 33, "El at rotator axis" );
      strncpy ( ctype2 + 10, "'EL---TAN'", 10 );
      break;

   case APPT:
   case FK5:
   case FK4:
      strcpy ( crval1 + 33, "RA at rotator axis" );
      strncpy ( ctype1 + 10, "'RA---TAN'", 10 );
      strcpy ( crval2 + 33, "Dec at rotator axis" );
      strncpy ( ctype2 + 10, "'DEC--TAN'", 10 );

      w1 = -w1;        /* For RA/Dec cases, CDELT1 will be -ve. */
      break;

   default:
      return -1;
   }

/* Pixel coordinates of reference point. */
   det = c * e - b * f;
   sprintf ( crpix1 + 10, "%20.12f", ( a * f - c * d ) / det );
   crpix1[30] = (char) ' ';
   sprintf ( crpix2 + 10, "%20.12f", ( b * d - a * e ) / det );
   crpix2[30] = (char) ' ';

/* Sky coordinates of reference point. */
   sprintf ( crval1 + 10, "%20.12f", slaDranrm ( wcsij.ab0[0] ) * R2D );
   crval1[30] = (char) ' ';
   sprintf ( crval2 + 10, "%+20.12f", wcsij.ab0[1] * R2D );
   crval2[30] = (char) ' ';

/* Scales. */
   sprintf ( cdelt1 + 10, "%20.12f", R2D * w1 );
   cdelt1[30] = (char) ' ';
   sprintf ( cdelt2 + 10, "%+20.12f", R2D * w2 );
   cdelt2[30] = (char) ' ';

/* PC-matrix. */
   sprintf ( pc001001 + 10, "%20.12f", b / w1 );
   pc001001[30] = (char) ' ';
   sprintf ( pc001002 + 10, "%20.12f", c / w1 );
   pc001002[30] = (char) ' ';
   sprintf ( pc002001 + 10, "%20.12f", e / w2 );
   pc002001[30] = (char) ' ';
   sprintf ( pc002002 + 10, "%20.12f", f / w2 );
   pc002002[30] = (char) ' ';


/* RA/Dec system. */
   switch ( frame ) {

   case AZEL_TOPO:
      radecsys[0] = (char) '\0';
      break;

   case APPT:
      strncpy ( radecsys + 10, "'GAPPT   '", 10 );
      break;

   case FK5:
      strncpy ( radecsys + 10, "'FK5     '", 10 );
      break;

   case FK4:
      strncpy ( radecsys + 10, "'FK4     '", 10 );
      break;
   }

/* Equinox. */
   switch ( frame ) {

   case FK4:
   case FK5:
      sprintf ( equinox + 10, "%20.12f",
                slaEpco ( frame == FK4 ? 'B' : 'J' , eqx.type, eqx.year ) );
      equinox[30] = (char) ' ';
      break;

   default:
      equinox[0] = (char) '\0';
   }

/* Epoch. */
   sprintf ( mjdobs + 10, "%20.12f", tt );
   mjdobs[30] = (char) ' ';

   return 0;
}
