#include "astLib.h"
int astCoco ( double a1, double b1, struct PMPXRV pmotion,
              FRAMETYPE cosys1, struct EPOCH equinox1, struct EPOCH epoch1,
              FRAMETYPE cosys2, struct EPOCH equinox2,
              double date, double aoprms[], struct TELP tel,
              double *a2, double *b2 )
/*
**  - - - - - - - -
**   a s t C o c o
**  - - - - - - - -
**
**  Transform a source position from one coordinate system to another.
**
**  GIVEN (arguments):
**
**   Source
**    a1        double     right ascension or azimuth (radians)
**    b1        double     declination or elevation (radians)
**    pmotion   PMPXRV     proper motion etc (see Notes 2 and 13)
**
**   Starting coordinate system
**    cosys1    FRAMETYPE  original coordinate system
**    equinox1  EPOCH      catalogue equinox (see Note 13)
**    epoch1    EPOCH      catalogue epoch (see Note 13)
**
**   Final coordinate system
**    cosys2    FRAMETYPE  final coordinate system
**    equinox2  EPOCH      equinox of final frame (see Note 13)
**
**   Circumstances
**    date      double     epoch of observation (TT MJD; see Note 1)
**    aoprms    double[15] star-independent appt-to-observed parameters
**                                                    (see Notes 5 and 13)
**    tel       TELP       telescope-dependent parameters (see Note 13)
**
**  RETURNED (arguments):
**    a2        double     final right ascension or azimuth (radians)
**    b2        double     final declination or elevation (radians)
**
**  RETURNED (function value):
**              int       +1 = unreachable position
**                         0 = OK
**                        -1 = invalid original frame
**                        -2 = invalid final frame
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
**        double pox;      / * Pointing-origin x-component (mm) * /
**        double poy;      / * Pointing-origin y-component (mm) * /
**  };
**
**  Called: astCocoR
**
**  Notes:
**
**    1  When using this routine, a potential source of confusion is
**       that there are two sorts of "now":
**
**          date        epoch of observation
**          aoprms[13]  sidereal time
**
**       In principle, the date and the sidereal time ought to be
**       consistent.  However, in practice this is unimportant:  the
**       date is used for calculating proper motion, precession etc,
**       while the sidereal time is used for converting between an
**       RA/Dec and an Az/El.  Note, however, that the sidereal time
**       must be accurate:  any vagueness will appear as corresponding
**       uncertainties in the rapidly-changing RA/Dec to/from Az/El
**       transformation.
**
**    2  The pmotion.pm flag distinguishes between the cases where a
**       proper motion in the (rotating) FK4 frame is being supplied and
**       the case where there is assumed to be no proper motion in a
**       (non-rotating) inertial frame.  False (0) means zero proper
**       motion in an inertial frame and is the correct value to use
**       for extragalactic sources.  True (not 0) is the correct choice
**       where there is a known proper motion, e.g. from a catalogue.
**       Note that in FK4 coordinates a distant object exhibits a
**       fictitious proper motion of up to about 0.5 arcseconds per
**       century.  Conversely, a star that happened to have a zero
**       proper motion in the FK4 system would in fact have a proper
**       motion of up to about 0.5 arcseconds per century with respect
**       to an inertial frame.  Setting pmotion.pm to false has no
**       effect for non-FK4 frames except that the parallax and radial
**       velocity are assumed zero.
**
**    3  The treatment of proper motion, parallax and radial velocity
**       is, strictly speaking, approximate in the case of FK4/FK5
**       frames of non-standard equinox (i.e. not B1950/J2000).  The
**       treatment here assumes that the proper motions are constant
**       in spherical coordinates and the radial velocity and parallax
**       are correct at B1950/J2000.
**
**    4  Input mean RA/Dec frames are barycentric (i.e. the coordinates
**       are unaffected by annual parallax).  Output mean RA/Dec frames
**       are geocentric (i.e. are subject to variations due to
**       annual parallax).
**
**    5  The aoprms array contains the star-independent apparent-to-observed
**       parameters:
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
**    6  Where the given position is one of the Az/Els, it is assumed
**       to refer to the sidereal time supplied in the aoprms array.
**       This means that the function cannot convert an Az/El from one
**       epoch to another.
**
**    7  There is no provision for different telescope parameters
**       applying to the original and final coordinate systems.
**       However, transformations involving different telescope
**       parameters (for example specifying two different pointing-
**       origins) can be carried out by calling the function twice,
**       first to transform into apparent place (say) and then to
**       transform back again with the new telescope parameters.
**
**    8  The transformations between topocentric and mount coordinates
**       assume tangent-plane geometry in the focal plane.
**
**    9  Azimuths are with respect to terrestrial rather than celestial
**       north.
**
**   10  The "unreachable position" status occurs when predicting
**       mount Az/El for a point near the zenith that cannot be reached
**       because of collimation error.  A reachable position nearby is
**       returned under these circumstances.
**
**   11  For transformations into mount Az/El, there may be cases near
**       the zenith where there are two solutions.  In such cases,
**       the function returns only the solution where the star and the
**       telescope are on the same side of the zenith.
**
**   12  Refraction corrections are disabled below the Gemini elevation
**       limit to reduce computation.  See the astCocoR routine for
**       further details.
**
**   13  Not all the arguments are referred to in all cases.  The
**       following rules apply:
**
**         argument    referred to if and only if
**
**         equinox1    cosys1 is FK4 or FK5
**         epoch1      cosys1 is FK4 or FK5
**         pmotion     cosys1 is FK4 or FK5
**         equinox2    cosys2 is FK4 or FK5
**         epoch2      cosys2 is FK4 or FK5
**         tel         cosys1 and/or cosys2 is AZEL_MNT
**         aoprms      both cosys1 and cosys2 are APPT, AZEL_TOPO or AZEL_MNT
**
**  P.T.Wallace   13 November 1997
**
**  Copyright RAL 1997.  All rights reserved.
*/

{
   double w1, w2;

/* Call the "with drift rates" version. */
   return astCocoR ( a1, b1, pmotion, 0.0, 0.0,
                       cosys1, equinox1, epoch1,
                       cosys2, equinox2, date, aoprms, tel,
                       a2, b2, &w1, &w2 );
}
