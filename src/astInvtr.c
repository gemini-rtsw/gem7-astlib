#include <slalib.h>
#include "astLib.h"
int astInvtr ( struct WCS wcs, struct WCS *iwcsp )
/*
**  - - - - - - - - -
**   a s t I n v t r
**  - - - - - - - - -
**
**  Invert a world coordinate system transformation, from x/y-to-sky
**  to sky-to-x/y or vice versa.
**
**  GIVEN:
**    wcs       struct     WCS transformation
**
**  RETURNED (argument):
**    iwcsp     struct*    pointer to inverse WCS transformation
**
**  RETURNED (function value):
**              int        0 = OK
**                        -1 = (bad status from slaInvf)
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
**  Called:  slaInvf
**
**  P.T.Wallace   20 March 1997
**
**  Copyright RAL 1997.  All rights reserved.
*/

{
   int j;

/* Copy the field centre. */
   iwcsp->ab0[0] = wcs.ab0[0];
   iwcsp->ab0[1] = wcs.ab0[1];

/* Invert the affine transformation. */
   slaInvf ( wcs.coeffs, iwcsp->coeffs, &j );
   if ( j ) return -1;

   return 0;
}
