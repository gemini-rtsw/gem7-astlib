#include <slalib.h>
#include "astLib.h"
int astFitij ( int n, double fpxy[][2], double pixij[][2],
               double cij[6], double *pixis, double *pixjs,
               double *perp, double *orient )
/*
**  - - - - - - - - -
**   a s t F i t i j
**  - - - - - - - - -
**
**  Fit instrumental i/j coordinates to focal-plane x/y coordinates.
**
**  GIVEN:
**    n         int             number of sample points (at least 3)
**    fpxy      double [n][2]   sets of x/y coordinates
**    pixij     double [n][2]   corresponding i/j coordinates
**
**  RETURNED (arguments):
**    cij       double[]        affine transformation, i/j to x/y
**    pixis     double*         i scale, x units per i unit
**    pixjs     double*         j scale, y units per j unit
**    perp      double*         i/j non-perpendicularity (radians)
**    orient    double*         orientation of i/j wrt x/y
**
**  RETURNED (function value):
**              int        0 = OK
**                        -1 = bug
**                        -2 = not enough points
**                        -3 = singular
**
** Notes:
**
**  1  The coefficients cij describe the affine transformation between
**     instrumental coordinates i/j (for example pixel coordinates in the
**     case of a CCD) and focal-plane coordinates x/y.  Writing cij[0-5]
**     as a,b,c,d,e,f:
**
**       x = a + bi + cj
**       y = d + ei + fj
**
**  2  The coefficients can be used by the astXtndtr function to extend
**     the WCS transformation so that it starts from instrumental i/j
**     instead of focal-plane x/y.  This transformation then relates
**     instrumental coordinates (e.g. pixels) directly to sky
**     coordinates (e.g. J2000 RA/Dec).
**
**  3  For i/j in pixels and x/y in mm, the scales pixis and pixjs are
**     in units of mm per pixel.  For a detector with square pixels,
**     these two numbers should be equal.
**
**  4  The nonperpendicularity, perp, should be zero for the usual case
**     where the instrumental i/j coordinates are orthogonal.
**
**  5  The orientation, orient, is the angle of the +j axis
**     anticlockwise of +y under the following circumstances:
**
**        The i and j axes appear right-handed.
**        The x and y axes appear left-handed.
**
**  6  The first two elements of the cij array are the pixel i/j
**     of the origin of the x/y coordinate system.  The x/y coordinates
**     corresponding to the centre of the pixel array can be found by
**     directly applying the i/j to x/y transformation (see note 1).
**
**  P.T.Wallace   15 November 1997
**
**  Copyright RAL 1997.  All rights reserved.
*/
{
   int j;
   double xz, yz;
 

   slaFitxy ( 6, n, fpxy, pixij, cij, &j );
   if ( j ) return j;
   slaDcmpf ( cij, &xz, &yz, pixis, pixjs, perp, orient );
   *orient *= -1.0;

   return 0;
}
