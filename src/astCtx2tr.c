#include <slalib.h>
#include <timeLib.h>
#include "astLib.h"
int astCtx2tr ( struct WCS_CTX ctx, FRAMETYPE frame,
                struct EPOCH equinox, double wavel, int ichop,
                struct WCS *wcsp, double *time )
/*
**  - - - - - - - - - -
**   a s t C t x 2 t r
**  - - - - - - - - - -
**
**  Generate a world coordinate system transformation.
**
**  GIVEN:
**    ctx       struct     WCS context
**    frame     FRAMETYPE  type of coordinate system
**    equinox   struct     equinox (mean RA/Decs only)
**    wavel     double     wavelength (microns)
**    ichop     int        chop state (0=A, 1=B, 2=C)
**
**  RETURNED (arguments):
**    wcsp      struct*    pointer to WCS structure
**    time      double*    Gemini raw time at which the transformation
**                         was correct
**
**  RETURNED (function value):
**              int        0 = OK
**                        -1 = illegal frame
**                        -2 = illegal ichop value
**                        -3 = internal error (bad status from timeThenD)
**                        -4 = internal error (bad status from astCoco)
**                        -5 = internal error (bad status from slaDs2tp)
**                        -6 = internal error (bad status from slaFitxy)
**
**  Defined in astLib.h:
**    WCS_CTX   struct     WCS context
**    WCS       struct     WCS transformation parameters
**    FRAMETYPE enum       coordinate system IDs
**    TELP      struct     telescope parameters
**    EPOCH     struct     B or J epoch
**    PMPXRV    struct     proper motion, parallax, radial velocity
**
**  Defined in timelib.h:
**    TT        timescale  ID for Terrestrial Time
**
**  Note:  The frame cannot be mount az/el because it refers to the
**         telescope rather than the sky.
**
**  To guide potential non-Gemini-TCS users of this routine, here are
**  the present definitions of the above.
**
**  struct WCS_CTX {
**        double ab0[2];      / * Pre-flexure Az/El * /
**        struct TELP tel;    / * Telescope parameters * /
**        double aoprms[15];  / * Apparent-to-observed parameters * /
**        double m2xy[3][2];  / * M2 tip/tilt (3 chop states) * /
**        double time;        / * Gemini raw time * /
**  };
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
**  struct TELP {
**        double fl;       / * Focal length (metres) * /
**        double rma;      / * Rotator orientation (radians) * /
**        double an;       / * Azimuth axis tilt NS (radians) * /
**        double aw;       / * Azimuth axis tilt EW (radians) * /
**        double pnpae;    / * Az/El nonperpendicularity (radians) * /
**        double ca;       / * LR collimation (radians) * /
**        double ce;       / * UD collimation (radians) * /
**        double pox;      / * not required by this routine * /
**        double poy;      / * not required by this routine * /
**  };
**
**  struct EPOCH {
**     double year;     / * Epoch:  Byear (B), TT Jyear (J) * /
**     char   type;     / * Type of epoch ('B', 'J' or ' ') * /
**  };
**
**  struct PMPXRV {
**     int    pm;       / * False = proper motion inertially-zero * /
**     double pmRA;     / * RA proper motion (radians/year) * /
**     double pmDec;    / * Dec proper motion (radians/year) * /
**     double px;       / * Parallax (arcsec) * /
**     double rv;       / * Radial velocity (km/s) * /
**  } ;
**
**   typedef enum {
**                  TAI,      / * International Atomic Time    * /
**                  UTC,      / * Coordinated Universal Time   * /
**                  UT1,      / * Universal Time               * /
**                  TT,       / * Terrestrial Time             * /
**                  TDB,      / * Barycentric Dynamical Time   * /
**                  GMST,     / * Greenwich Mean Sidereal Time * /
**                  LAST      / * Local Apparent Sidereal Time * /
**                } timescale;
**
**  Called:  timeThenD, astCoco, slaAtmdsp, slaDs2tp, slaFitxy
**
**  P.T.Wallace   22 March 1999
**
**  Copyright RAL 1999.  All rights reserved.
**
**  21 Oct 02 - extend context to include optical distortion of telescope (CJM)
*/

#define GRID 100.0    /* Spacing for sample points (mm) */

{
/* Null proper-motion. */
   struct PMPXRV pm = { 0 };

/* Sample points in the focal plane (x/y in mm) */
   double xy[5][2] = {
      {  0.0,  0.0 },          /* First point must be the origin */
      {-GRID,-GRID },
      {-GRID, GRID },
      { GRID,-GRID },
      { GRID, GRID }
   };

/* Corresponding WCS standard coordinates */
   double xieta[5][2];

   double date, a, b;
   int i, j;
   double distort[6] ;
   struct WCS wcstp ;

/* Disallow mount az/el. */
   if ( frame == AZEL_MNT ) return -1;

/* Adjust the collimation coefficients for M2 tip/tilt. */
   if ( ichop >= 0 && ichop <= 2 ) {
      ctx.tel.ca -= ctx.m2xy[ichop][0];
      ctx.tel.ce -= ctx.m2xy[ichop][1];
   } else {
      return -2;
   }

/* Convert the raw time to TT. */
   if ( timeThenD ( ctx.time, TT, &date ) ) return -3;

/* Adjust the refraction constants for colour. */
   a = ctx.aoprms[10];
   b = ctx.aoprms[11];
   slaAtmdsp ( ctx.aoprms[5],
               ctx.aoprms[6],
               ctx.aoprms[7],
               ctx.aoprms[8],
               a, b, wavel,
               & ctx.aoprms[10], &ctx.aoprms[11] );

/* Look at five points in the focal plane. */
   for ( i = 0; i < 5; i++ ) {

   /* The point. */
      ctx.tel.pox = xy[i][0];
      ctx.tel.poy = xy[i][1];

   /* World coordinates of the point projected onto the sky. */
      if ( astCoco ( ctx.ab0[0], ctx.ab0[1], pm,
                AZEL_MNT, equinox, equinox,
                frame, equinox,
                date, ctx.aoprms, ctx.tel,
                &a, &b ) ) return -4;

   /* The field centre. */
      if ( i == 0 ) {
         wcstp.ab0[0] = a;
         wcstp.ab0[1] = b;
      }

   /* To standard coordinates. */
      slaDs2tp ( a, b, wcstp.ab0[0], wcstp.ab0[1],
                 &xieta[i][0], &xieta[i][1], &j );
      if ( j ) return -5;
   }

/* Fit 6-coefficient linear model to the x/y and xi/eta samples. */
   slaFitxy ( 6, 5, xieta, xy, wcstp.coeffs, &j );
   if ( j ) return -6;

/* Extend to include optical distortion */
   astGetdistortion(distort);
   astXtndtr (distort, wcstp, wcsp) ;

/* Export the timestamp. */
   *time = ctx.time;

   return 0;
}
