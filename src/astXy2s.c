#include "astLib.h"
int astXy2s ( double x, double y, FRAMETYPE frame,
              struct EPOCH equinox, double wavel, int ichop,
              double *a, double *b )
/*
**  - - - - - - - -
**   a s t X y 2 s
**  - - - - - - - -
**
**  Transform a focal-plane x/y position into a specified celestial
**  coordinate system.
**
**  GIVEN:
**    x,y       double     x/y in focal plane (mm)
**    frame     FRAMETYPE  type of coordinate system
**    equinox   struct     equinox (mean RA/Decs only)
**    wavel     double     wavelength (microns)
**    ichop     int        chop state (0=A, 1=B, 2=C)
**
**  RETURNED (argument):
**    a,b       double*    RA/Dec or Az/El (radians)
**
**  RETURNED (function value):
**              int        0 = OK
**                        -1 = error from astGettr
**                        -2 = error from astXy2sq
**
**  Defined in astLib.h:
**    WCS       struct     WCS transformation parameters
**    FRAMETYPE enum       coordinate system IDs
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
**  Called:  astGettr, astXy2sq
**
**  Notes:
**
**   1  The x/y coordinates are in the Cassegrain focal-plane and
**      rotate with the instrument mount.
**
**   2  Azimuths are north-through-east.
**
**  P.T.Wallace  30 January 1997
**
**  Copyright RAL 1997.  All rights reserved.
*/
{
   struct WCS wcs;
   double time;


/* Obtain WCS parameters. */
   if ( astGettr ( frame, equinox, wavel, ichop, &wcs, &time ) ) return -1;

/* Focal-plane coordinates to celestial coordinates. */
   if ( astXy2sq ( x, y, wcs, a, b ) ) return -2;

   return 0;
}
