#include <slalib.h>
#include "astLib.h"
int astRot ( FRAMETYPE cosys, struct EPOCH eqx,
             double ra, double de,
             double date, double aoprms[], struct TELP tel,
             double aia, double pai, double *rma )
/*
**  - - - - - - -
**   a s t R o t
**  - - - - - - -
**
**  Predict the rotator orientation required to achieve a given
**  field orientation.
**
**  GIVEN (arguments):
**    cosys     FRAMETYPE  reference frame...
**    eqx       EPOCH      ...for specifying field orientation
**    ra        double     RA or azimuth in that frame (radians)
**    de        double     Dec or elevation in that frame (radians)
**    date      double     epoch of observation (TT MJD)
**    aoprms    double[15] star-independent appt-to-observed parameters
**    tel       TELP       telescope-dependent parameters
**    aia       double     IAA: instrument alignment angle (Notes 1,2)
**    pai       double     IPA: desired sky PA of IPD (Notes 1,3)
**
**  RETURNED (argument):
**    rma       double     required rotator mechanical angle (range +/- pi)
**
**  RETURNED (function value):
**              int        0 = OK
**                        -1 = unreachable position
**                        -2 = invalid cosys or other astCoco error
**                        -3 = failed to converge
**
**  Defined in astLib.h:
**    FRAMETYPE enum       coordinate system IDs
**    EPOCH     struct     B or J epoch
**    PMPXRV    struct     proper motion, parallax, radial velocity
**    TELP      struct     telescope parameters
**
**  To guide potential non-Gemini-TCS users of this routine, here are the
**  present definitions of the above.
**
**  typedef enum { AZEL_MNT = 0,   / * Mount Az/El, pre-flexure (illegal) * /
**                 AZEL_TOPO = 1,  / * Topocentric Az/El * /
**                 APPT = 2,       / * Geocentric apparent RA/Dec * /
**                 FK5 = 3,        / * IAU 1976 RA/Dec, any equinox * /
**                 FK4 = 4         / * Pre IAU 1976 RA/Dec, any equinox * /
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
**        double pox;      / * Pointing-origin x-component (mm) * /
**        double poy;      / * Pointing-origin y-component (mm) * /
**  };
**
**  Called: astCoco, slaEpj, slaDbear, slaDrange
**
**  Notes:
**
**    1  The orientation of the instrument is specified through some
**       chosen direction called the "instrument principal direction"
**       (IPD).
**
**    2  The instrument alignment angle, aia, is the orientation of
**       the IPD with respect to the rotator +y axis, in radians;  it
**       increases anticlockwise on the sky.  The usual name of the
**       angle is IAA.
**
**    3  The desired position-angle, pai, is that of the IPD at the
**       pointing-origin, as seen projected on the sky, in radians.
**       It is zero when the IPD points north (or up) in the specified
**       reference frame, and it increases anti-clockwise on the sky.
**       The usual name of the angle is IPA.
**
**    4  When using this routine, a potential source of confusion is
**       that there are two sorts of "now":
**
**          date        epoch of observation
**          aoprms[13]  sidereal time
**
**       In principle, the date and the sidereal time ought to be
**       consistent.  However, in practice this is unimportant:  the
**       date is used for calculating precession etc, while the sidereal
**       time is used for converting between RA/Dec and Az/El.  Note,
**       however, that the sidereal time will need to be accurate if the
**       field is rapidly rotating as the telescope tracks.
**
**    5  In most circumstances, one coordinate system will be used for
**       (a) supplying the target coordinates, (b) controlling the
**       telescope tracking and (c) specifying the field orientation.
**       In particular, it will be usual for star positions to be
**       supplied as J2000 RA/Dec, for the telescope to be tracking
**       in J2000 RA/Dec and for the field to be oriented with
**       respect to north in J2000 RA/Dec.  However, there may also
**       be circumstances where two or three different coordinate
**       systems are involved:  the star position may be supplied as
**       B1950 RA/Dec, and the position angle may be fixed in Az/El
**       for example.  Where this is so, preliminary calls to the
**       astCoco function will be required before the present routine
**       can be used.
**
**    6  The aoprms array contains the star-independent apparent-to-
**       observed parameters:
**
**       [0]      geodetic latitude (radians)
**       [1,2]    sine and cosine of geodetic latitude
**       [3]      magnitude of diurnal aberration vector
**       [4]      height (metres)
**       [5]      ambient temperature (degK)
**       [6]      pressure (mB)
**       [7]      relative humidity (0-1)
**       [8]      effective wavelength (microns)
**       [9]      tropospheric lapse rate (degK/metre)
**       [10,11]  refraction constants A and B (radians)
**       [12]     longitude + eqn of equinoxes + sidereal DUT (radians)
**       [13]     local apparent sidereal time (radians)
**       [14]     polar-motion adjustment to azimuth (radians)
**
**       It can be constructed by calling the SLALIB routine slaAoppa
**       to generate elements 0-13 and then slaPolmo to generate the
**       final element.  The sidereal time (element 13) can be
**       selectively updated by calling slaAoppat.  Adjustment
**       of the refraction constants can be accomplished most
**       efficiently by using slaAtmdsp to provide new values for
**       elements 8, 10 and 11.
**
**       Element [12] is not used by the present routine.  It is
**       included for compatibility with slaAoppat.
**
**    7  The function always uses the telescope parameters and so a
**       valid structure must be supplied.  However, for many purposes
**       nominal values for the various telescope parameters will be
**       perfectly acceptable.
**
**    8  The transformations between topocentric and mount coordinates
**       assume tangent-plane geometry in the focal plane.
**
**    9  Azimuths are with respect to terrestrial rather than celestial
**       north.
**
**  P.T.Wallace   23 April 1998
**
**  Copyright RAL 1998.  All rights reserved.
*/

#define NITS 100                  /* Maximum number of iterations */
#define DMAX 1e-6                 /* Convergence criterion (radians) */
#define DELTA 1.0                 /* Trial shift in pointing-origin (mm) */

{
   struct PMPXRV pm = { 0, 0.0, 0.0, 0.0, 0.0 };
   struct EPOCH ep;
   double x, y, dx, dy, x1, y1, x2, y2, az, el, ra_py, de_py,
          ra_my, de_my, eps, dtheta, factor;
   int it, j;


/* Express the TT date as a Julian epoch. */
   ep.type = 'J';
   ep.year = slaEpj ( date );

/* The actual pointing-origin. */
   x = tel.pox;
   y = tel.poy;

/* Two trial pointing-origins either side. */
   dx = DELTA * sin ( aia );
   dy = DELTA * cos ( aia );
   x1 = x + dx;
   y1 = y + dy;
   x2 = x - dx;
   y2 = y - dy;

/* Initialize iteration count. */
   it = 0;

/* Loop until converged. */
   do {

   /* Error exit if iteration limit reached. */
      if ( it >= NITS ) return -3;

   /* Transform the source position to pre-flexure mount az/el using  */
   /* the current rotator mechanical angle and given pointing-origin. */
      tel.pox = x;
      tel.poy = y;
      j = astCoco ( ra, de, pm, cosys, eqx, ep,
                    AZEL_MNT, eqx, date, aoprms, tel,
                    &az, &el );
      if ( j == 1 ) return -1;
      if ( j ) return -2;

   /* Transform the mount position back to sky coordinates, but this */
   /* time using a different pointing origin, displaced in +IPD.     */
      tel.pox = x1;
      tel.poy = y1;
      if ( astCoco ( az, el, pm, AZEL_MNT, eqx, ep,
                     cosys, eqx, date, aoprms, tel,
                     &ra_py, &de_py ) ) return -2;

   /* The same thing for a -IPD displacement. */
      tel.pox = x2;
      tel.poy = y2;
      if ( astCoco ( az, el, pm, AZEL_MNT, eqx, ep,
                     cosys, eqx, date, aoprms, tel,
                     &ra_my, &de_my ) ) return -2;

   /* Calculate the position-angle between the two points and hence */
   /* how far away from the desired field orientation we are.       */
      eps = slaDrange ( pai - slaDbear ( ra_my, de_my, ra_py, de_py )
                                 * ( cosys == AZEL_TOPO ? -1.0 : 1.0 ) );
      

   /* Adjust the correction to speed up convergence. */
      factor = it ? factor * ( dtheta + eps ) / dtheta : 1.0;
      dtheta = eps * factor;

   /* Update the rotator prediction. */
      tel.rma += dtheta;

   /* Increment the iteration count. */
      it++;

   /* Round again if the last error estimate was still significant. */
   } while ( fabs ( eps ) > DMAX );

/* Return the final rotator position-angle. */
   *rma = slaDrange ( tel.rma );
   return 0;
}
