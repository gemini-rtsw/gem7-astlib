#include "astLib.h"
int astS2xy ( double a, double b, FRAMETYPE frame,
              struct EPOCH equinox, double wavel, int ichop,
              double *x, double *y )
/*
**  - - - - - - - -
**   a s t S 2 x y
**  - - - - - - - -
**
**  Transform celestial coordinates in a specified celestial coordinate
**  system into a focal-plane x/y position.
**
**  GIVEN:
**    a,b       double*    RA/Dec or Az/El (radians)
**    frame     FRAMETYPE  type of coordinate system
**    equinox   struct     equinox (mean RA/Decs only)
**    wavel     double     wavelength (microns)
**    ichop     int        chop state (0=A, 1=B, 2=C)
**
**  RETURNED (argument):
**    x,y       double     x/y in focal plane (mm)
**
**  RETURNED (function value):
**              int        0 = OK
**                        -1 = error from astGettr
**                        -2 = error from astInvtr
**                        -3 = error from astS2xyq
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
**  Called:  astGettr, astInvtr, astS2xyq
**
**  Notes:
**
**   1  Azimuths are north-through-east.
**
**   2  The x/y coordinates are in the Cassegrain focal-plane and
**      rotate with the instrument mount.
**
**  P.T.Wallace  30 January 1997
**
**  Copyright RAL 1997.  All rights reserved.
*/
{
   struct WCS wcs, iwcs;
   double time;


/* Obtain WCS parameters. */
   if ( astGettr ( frame, equinox, wavel, ichop, &wcs, &time ) ) return -1;

/* Invert the transformation. */
   if ( astInvtr ( wcs, &iwcs ) ) return -2;

/* Celestial coordinates to focal-plane coordinates. */
   if ( astS2xyq ( a, b, iwcs, x, y ) ) return -3;

   return 0;
}
