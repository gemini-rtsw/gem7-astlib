#include <string.h>
#include <slalib.h>
#include "astsys.h"
#include "astLib.h"
int astFITSv ( struct WCS wcsij,
               FRAMETYPE frame, struct EPOCH eqx, double djm,
               char* ctype1, double* crpix1, double* crval1,
               char* ctype2, double* crpix2, double* crval2,
               double* cd1_1, double* cd1_2,
               double* cd2_1, double* cd2_2,
               char* radecsys, double* equinox, double* mjdobs )
/*
**  - - - - - - - - -
**   a s t F I T S v
**  - - - - - - - - -
**
**  Express a world coordinate system transformation as FITS keyword
**  values.
**
**  GIVEN:
**    wcsij     struct     WCS transformation, pixels to sky
**    frame     FRAMETYPE  type of sky coordinate system
**    eqx       struct     equinox (mean RA/Decs only)
**    djm       double     epoch of observation (MJD: JD-2400000.5)
**
**  RETURNED (arguments):
**    ctype1    char[9]    value for FITS keyword CTYPE1
**    crpix1    double*    value for FITS keyword CRPIX1
**    crval1    double*    value for FITS keyword CRVAL1
**    ctype2    char[9]    value for FITS keyword CTYPE1
**    crpix2    double*    value for FITS keyword CRPIX2
**    crval2    double*    value for FITS keyword CRVAL2
**    cd1_1     double*    value for FITS keyword CD1_1
**    cd1_2     double*    value for FITS keyword CD1_2
**    cd2_1     double*    value for FITS keyword CD2_1
**    cd2_2     double*    value for FITS keyword CD2_2
**    radecsys  char[9]    value for FITS keyword RADECSYS
**    equinox   double*    value for FITS keyword EQUINOX
**    mjdobs    double*    value for FITS keyword MJD-OBS
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
**  Called:  slaDranrm, slaEpco
**
**  Notes:
**
**  1  The string arguments, ctype1, ctype2 and radecsys, are exactly
**     nine characters long (i.e. eight characters plus the terminating
**     null).  They are left-justified, with trailing spaces.
**
**  2  The keywords EQUINOX and/or RADECSYS may not be required.  In
**     such cases, a negative value is returned for equinox and an
**     empty string for radecsys.
**
**  3  The timescale for the epoch of observation, djm, is not defined
**     by the current FITS standard.  TAI is a suitable choice;  UTC
**     is not, strictly speaking, because expressing UTC as MJD leads
**     to ambiguities during leap seconds.
**
**  4  The present function returns values.  Another function, astFITSs,
**     returns strings.  (Note that the argument djm is simply copied to
**     the argument mjd-obs as it is.  It is included merely to make the
**     calls to the two functions match.)
**
**  5  Here is an example of the set of 13 values returned by this routine:
**
**     ctype1    "RA---TAN"            gnomonic projection
**     crpix1    984.048024            pixel i-coordinate at rotator axis
**     crval1    36.00179953283        RA or El at rotator axis
**     ctype2    "DEC--TAN"            gnomonic projection
**     crpix2    438.647662            pixel j-coordinate at rotator axis
**     crval2    44.99795763567        Dec or El at rotator axis
**     cd1_1     -0.00001007236440932  xi rotation/skew/scale matrix element
**     cd1_2     0.00000000574196532   xj rotation/skew/scale matrix element
**     cd2_1     -0.00000000499131497  yi rotation/skew/scale matrix element
**     cd2_2     0.00001012093669430   yj rotation/skew/scale matrix element
**     radecsys  "FK5     "            type of RA/Dec
**     equinox   2000.0000             epoch of mean equator & equinox
**     mjd-obs   49560.6433912037      epoch of observation (TAI MJD)
**
**  P.T.Wallace   22 November 1998
**
**  Copyright RAL 1998.  All rights reserved.
*/

{
   double a, b, c, d, e, f, det;


/* Axis types. */
   switch ( frame ) {

   case AZEL_TOPO:
      strcpy ( ctype1, "AZ---TAN" );
      strcpy ( ctype2, "EL---TAN" );
      break;

   case APPT:
   case FK5:
   case FK4:
      strcpy ( ctype1, "RA---TAN" );
      strcpy ( ctype2, "DEC--TAN" );
      break;

   default:
      return -1;
   }

/* Pixel coordinates of reference point. */
   a = wcsij.coeffs[0];
   b = wcsij.coeffs[1];
   c = wcsij.coeffs[2];
   d = wcsij.coeffs[3];
   e = wcsij.coeffs[4];
   f = wcsij.coeffs[5];
   det = c * e - b * f;
   *crpix1 = ( a * f - c * d ) / det;
   *crpix2 = ( b * d - a * e ) / det;

/* Sky coordinates of reference point. */
   *crval1 = slaDranrm ( wcsij.ab0[0] ) * R2D;
   *crval2 = wcsij.ab0[1] * R2D;

/* CD-matrix. */
   *cd1_1 = b * R2D;
   *cd1_2 = c * R2D;
   *cd2_1 = e * R2D;
   *cd2_2 = f * R2D;

/* RA/Dec system. */
   switch ( frame ) {

   case AZEL_MNT:
   case AZEL_TOPO:
      radecsys[0] = (char) '\0';
      break;

   case APPT:
      strcpy ( radecsys, "GAPPT   " );
      break;

   case FK5:
      strcpy ( radecsys, "FK5     " );
      break;

   case FK4:
      strcpy ( radecsys, "FK4     " );
      break;
   }

/* Equinox. */
   switch ( frame ) {

   case FK4:
   case FK5:
      *equinox = slaEpco ( frame == FK4 ? 'B' : 'J' , eqx.type, eqx.year );
      break;

   default:
      *equinox = -1.0;
   }

/* Epoch. */
   *mjdobs = djm;

   return 0;
}
