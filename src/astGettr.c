#include "astLib.h"
int astGettr ( FRAMETYPE frame, struct EPOCH equinox, double wavel,
               int ichop, struct WCS *wcsp, double *time )
/*
**  - - - - - - - - -
**   a s t G e t t r
**  - - - - - - - - -
**
**  Obtain the current world coordinate system transformation (x/y to
**  sky).
**
**  GIVEN:
**    frame     FRAMETYPE  type of coordinate system
**    equinox   struct     equinox (mean RA/Decs only)
**    wavel     double     wavelength (microns)
**    ichop     int        chop state (0=A, 1=B, 2=C)
**
**  RETURNED (argument):
**    wcsp      struct*    pointer to WCS structure (x/y to sky)
**    time      double*    Gemini raw time at which the transformation
**                         was correct
**
**  RETURNED (function value):
**              int        0 = OK
**                        -1 = bad status from astGetctx
**                        -2 = bad status from astCtxtr
**
**  Defined in astLib.h:
**    WCS_CTX   struct     WCS context
**    WCS       struct     WCS transformation parameters
**    FRAMETYPE enum       coordinate system IDs
**    TELP      struct     telescope parameters
**    EPOCH     struct     B or J epoch
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
**        double pox;      / * not used by this routine * /
**        double poy;      / * not used by this routine * /
**  };
**
**  struct EPOCH {
**        double year;     / * Epoch:  Byear (B), TT Jyear (J) * /
**        char   type;     / * Type of epoch ('B', 'J' or ' ') * /
**  };
**
**  Called:  astGetctx, astCtx2tr
**
**  P.T.Wallace   13 November 1997
**
**  Copyright RAL 1997.  All rights reserved.
*/
{
/* WCS context */
   struct WCS_CTX ctx;


/* Grab the context. */
   if ( astGetctx ( &ctx ) ) return -1;

/* Generate the WCS transformation. */
   if ( astCtx2tr ( ctx, frame, equinox, wavel, ichop, wcsp, time ) )
      return -2;

   return 0;
}
