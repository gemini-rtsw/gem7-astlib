#include <slalib.h>
#include "astLib.h"
int astS2xyq ( double a, double b, struct WCS iwcs, double *x, double *y )
/*
**  - - - - - - - - -
**   a s t S 2 x y q
**  - - - - - - - - -
**
**  Transform celestial coordinates into a focal-plane x/y position
**  using pre-computed WCS parameters.
**
**  GIVEN:
**    a,b     double     RA/Dec or Az/El (radians)
**    iwcs    struct     WCS parameters for the inverse transformation
**
**  RETURNED (argument):
**    x,y     double*    x/y in focal plane (mm)
**
**  RETURNED (function value):
**            int        0 = OK
**                      -1 = error, star too far from axis
**                      -2 = error, antistar on tangent plane
**                      -3 = error, antistar too far from axis
**
**  Defined in astLib.h:
**    WCS       struct     WCS transformation parameters
**
**  To guide potential non-Gemini-TCS users of this routine, here is
**  the present definition of the above.
**
**  struct WCS {
**        double ab0[2];      / * Celestial coordinates at x=y=0 * /
**        double coeffs[6];   / * Affine transformation coefficients * /
**  };
**
**  Called:  slaDs2tp, slaXy2xy
**
**  Notes:
**
**   1  Azimuths are north-through-east.
**
**   2  The x/y coordinates are in the Cassegrain focal-plane and
**      rotate with the instrument mount.
**
**  P.T.Wallace   20 March 1997
**
**  Copyright RAL 1997.  All rights reserved.
*/
{
   int j;
   double xi, eta;


/* Celestial coordinates to standard coordinates. */
   slaDs2tp ( a, b, iwcs.ab0[0], iwcs.ab0[1], &xi, &eta, &j );
   if ( j ) return -j;

/* Standard coordinates to focal-plane coordinates. */
   slaXy2xy ( xi, eta, iwcs.coeffs, x, y );

   return 0;
}
