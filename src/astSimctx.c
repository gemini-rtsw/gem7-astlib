#include <slalib.h>
#include <timeLib.h>
#include "astLib.h"
int astSimctx ( double tai,
                double elongm, double phim, double hm,
                double xp, double yp,
                double tdk, double pmb, double rh, double tlr,
                double wl, struct TELP tel, double m2xy[3][2],
                double a1, double b1, FRAMETYPE cosys, struct EPOCH equinox,
                struct WCS_CTX *ctxp )
/*
**  - - - - - - - - - -
**   a s t S i m c t x
**  - - - - - - - - - -
**
**  Generate a simulated context from which world coordinate system
**  transformations can be obtained.  The mechanical orientation of
**  the rotator is known.
**
**  This function can be used online, in an IOC, or offline, in a
**  Unix host.
**
**  Normally, the context is obtained from the telescope control
**  system and applies to the circumstances prevailing at the time.
**  The present routine allows a context to be generated for a
**  different set of circumstances.  This technique has both online
**  and offline uses.  Online, it can be used to simulate a future
**  observation, for example so that guide probe coordinates can
**  be predicted.  Offline, it has applications in reconstructing
**  past observations and in pre-observing planning.
**
**  GIVEN (arguments):
**
**   Time and place
**    tai       double       TAI (MJD)
**    elongm    double       east longitude (mean, radians)
**    phim      double       latitude (mean geodetic, radians)
**    hm        double       height above reference spheroid (metres)
**    xp        double       polar motion x-coordinate (radians)
**    yp        double       polar motion y-coordinate (radians)
**    tdk       double       ambient temperature (deg K)
**    pmb       double       pressure (millibar)
**    rh        double       relative humidity (range 0-1)
**    tlr       double       tropospheric lapse rate (degK/metre)
**
**   Telescope
**    wl        double       effective wavelength (microns)
**    tel       TELP         telescope-dependent parameters
**    m2xy      double[3][2] M2 tip/tilts for chop A/B/C
**
**   Target
**    a1        double       right ascension or azimuth (radians)
**    b1        double       declination or elevation (radians)
**    cosys     FRAMETYPE    coordinate system
**    equinox   EPOCH        equinox (mean RA/Dec only)
**
**  RETURNED (argument):
**    ctxp      WCS_CTX*     WCS context
**
**  RETURNED (function value):
**              int          0 = OK 
**                          -1 = error calling timeThenD
**                          -2 = error calling astCoco
**
**  Defined in astLib.h:
**    WCS_CTX   struct       WCS context
**    FRAMETYPE enum         coordinate system IDs
**    EPOCH     struct       B or J epoch
**    PMPXRV    struct       proper motion, parallax, radial velocity
**    TELP      struct       telescope parameters
**
**  To guide potential non-Gemini-TCS users of this routine, here are
**  the present definitions of the above.
**
**  typedef enum { AZEL_MNT = 0,      / * Mount Az/El, pre-flexure * /
**                 AZEL_TOPO = 1,     / * Topocentric Az/El * /
**                 APPT = 2,          / * Geocentric apparent RA/Dec * /
**                 FK5 = 3,           / * IAU 1976 RA/Dec, any equinox * /
**                 FK4 = 4            / * Pre IAU 1976 RA/Dec, any equinox * /
**               } FRAMETYPE;
**
**  struct EPOCH {
**        double year;     / * Epoch:  Byear (B), TT Jyear (J) * /
**        char   type;     / * Type of epoch ('B', 'J' or ' ') * /
**  };
**
**  struct PMPXRV {
**        int    pm;       / * False = proper motion inertially-zero * /
**        double pmRA;     / * RA proper motion (radians/year) * /
**        double pmDec;    / * Dec proper motion (radians/year) * /
**        double px;       / * Parallax (arcsec) * /
**        double rv;       / * Radial velocity (km/s) * /
**  };
**
**  struct TELP {
**        double fl;       / * Focal length (metres) * /
**        double rma;      / * Rotator orientation (radians) * /
**        double an;       / * Azimuth axis tilt NS (radians) * /
**        double aw;       / * Azimuth axis tilt EW (radians) * /
**        double pnpae;    / * Az/El nonperpendicularity (radians) * /
**        double ca;       / * LR collimation (radians) * /
**        double ce;       / * UD collimation (radians) * /
**        double pox;      / * not used by present routine * /
**        double poy;      / * not used by present routine * /
**  };
**
**  struct WCS_CTX {
**        double ab0[2];      / * Pre-flexure Az/El * /
**        struct TELP tel;    / * Telescope parameters * /
**        double aoprms[15];  / * Apparent-to-observed parameters * /
**        double m2xy[3][2];  / * M2 tip/tilt (3 chop states) * /
**        double time;        / * Gemini raw time * /
**  };
**
**  Called:  timeTai2raw, timeThenD, slaEpj, slaAoppa, slaPolmo, astCoco
**
**  Notes:
**
**    1  For the present function to work, the time system must be
**       operating.  For online use, this means that timeInit must
**       have been called.  For offline use, timeOffline must have
**       been called.  An attempt to use the present function with
**       previously having started the time system produces an
**       error status.
**
**    2  The structure element aoprms[12] is not used by the other
**       WCS functions and is set to zero.
**
**    3  When setting up the tel telescope-parameter-structure, note
**       that the rotator angle is the mechanical one, i.e. with respect
**       to the vertical rather than north.  When starting from the
**       orientation of the field in the focal plane, either use the
**       function astRot to determine the rotator angle that will be
**       needed, or call the companion routine, astSimctx_r, which
**       directly addresses this case.
**
**    4  The target coordinates [a1,b1] are for the epoch of observation
**       (i.e. tai).
**
**    5  If the pressure is not known, it can be estimated using the
**       expression pmb=1013.25*exp(-hm/(29.3*tsl)), where tsl is the
**       approximate air temperature at sea-level, in degrees K.
**
**  P.T.Wallace   23 April 1998
**
**  Copyright RAL 1998.  All rights reserved.
*/
{
/* Parameters for no proper motion */
   struct PMPXRV pm = { 0, 0.0, 0.0, 0.0, 0.0 };

/* Current epoch */
   struct EPOCH epoch = { 0.0, 'J' };

   double raw, tt, ut1, elong, phi;


/* Express the time in various forms. */
   raw = timeTai2raw ( tai );
   if ( timeThenD ( raw, TT, &tt ) ) return -1;
   epoch.year = slaEpj ( tt );
   if ( timeThenD ( raw, UT1, &ut1 ) ) return -1;

/* Copy the telescope parameters. */
   ctxp->tel = tel;

/* Generate the apparent-to-observed parameters. */
   slaAoppa ( ut1, 0.0, elongm, phim, hm, xp, yp,
              tdk, pmb, rh, wl, tlr, ctxp->aoprms );
   slaPolmo ( elongm, phim, xp, yp, &elong, &phi, &ctxp->aoprms[14] );

/* Generate the pre-flexure mount coordinates. */
   if ( astCoco ( a1, b1, pm,
                  cosys, equinox, epoch,
                  AZEL_MNT, equinox, tt,
                  ctxp->aoprms, tel,
                  &ctxp->ab0[0], &ctxp->ab0[1] ) ) return -2;

/* Copy the m2 tip tilts. */
   ctxp->m2xy[0][0] = m2xy[0][0];
   ctxp->m2xy[0][1] = m2xy[0][1];
   ctxp->m2xy[1][0] = m2xy[1][0];
   ctxp->m2xy[1][1] = m2xy[1][1];
   ctxp->m2xy[2][0] = m2xy[2][0];
   ctxp->m2xy[2][1] = m2xy[2][1];

/* Copy the timestamp. */
   ctxp->time = raw;

   return 0;
}
