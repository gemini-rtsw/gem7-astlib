#include <slalib.h>
#include "astLib.h"
int astXy2sq ( double x, double y, struct WCS wcs, double *a, double *b )
/*
**  - - - - - - - - -
**   a s t X y 2 s q
**  - - - - - - - - -
**
**  Transform a focal-plane x/y position into celestial coordinates
**  using pre-computed WCS parameters.
**
**  GIVEN:
**    x,y     double     x/y in focal plane (mm)
**    wcs     struct     WCS parameters
**
**  RETURNED (argument):
**    a,b     double*    RA/Dec or Az/El (radians)
**
**  RETURNED (function value):
**            int        0 = OK (always)
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
**  Called:  slaXy2xy, slaDtp2s
**
**  Notes:
**
**   1  The x/y coordinates are in the Cassegrain focal-plane and
**      rotate with the instrument mount.
**
**   2  Azimuths are north-through-east.
**
**  P.T.Wallace   20 March 1997
**
**  Copyright RAL 1997.  All rights reserved.
*/
{
   double xi, eta;


/* Focal-plane coordinates to standard coordinates. */
   slaXy2xy ( x, y, wcs.coeffs, &xi, &eta );

/* Standard coordinates to celestial coordinates. */
   slaDtp2s ( xi, eta, wcs.ab0[0], wcs.ab0[1], a, b );

   return 0;
}
