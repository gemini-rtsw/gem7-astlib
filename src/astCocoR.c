#include <string.h>
#include <slalib.h>
#include "astLib.h"
#include "astsys.h"
int astCocoR ( double a1, double b1, struct PMPXRV pmotion,
               double adot1, double bdot1,
               FRAMETYPE cosys1, struct EPOCH equinox1, struct EPOCH epoch1,
               FRAMETYPE cosys2, struct EPOCH equinox2,
               double date, double aoprms[], struct TELP tel,
               double *a2, double *b2,
               double *adot2, double *bdot2 )
/*
**  - - - - - - - - -
**   a s t C o c o R
**  - - - - - - - - -
**
**  Transform a source position and a drift-rate from one coordinate
**  system to another.
**
**  GIVEN (arguments):
**
**   Source
**    a1        double     right ascension or azimuth (radians)
**    b1        double     declination or elevation (radians)
**    pmotion   PMPXRV     proper motion etc (see Notes 2 and 13)
**
**   Drift-rate
**    adot1     double     drift-rate in RA or Az (radians per SI second)
**    bdot1     double     drift-rate in Dec or El (radians per SI second)
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
**    adot2     double     final drift-rate in RA or Az (radians per SI sec)
**    bdot2     double     final drift-rate in Dec or El (radians per SI sec)
**
**  RETURNED (function value):
**              int       +1 = unreachable position
**                         0 = OK
**                        -1 = invalid original frame
**                        -2 = invalid final frame
**
**  Defined in astsys.h:
**    AS2R      double     1 arcsecond in radians
**    S2R       double     1 second in radians
**    D90       double     90 degrees in radians
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
**  Called: slaEpj, slaEpco, slaPm, slaPrebn, slaEtrms, slaDcs2c,
**          slaDmxv, slaDcc2s, slaFk45z, slaPrec, slaMappa,
**          slaFk54z, slaMapqk, slaAmpqk, slaAopqk, slaOapqk,
**          slaDranrm, slaDrange
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
**       the function only the solution where the star and the
**       telescope are on the same side of the zenith.
**
**   12  The routine includes special precautions to avoid expensive
**       and pointless refraction calculations for cases below the Gemini
**       15 degree elevation limit.  These precautions consist of disabling
**       the refraction corrections when transforming to or from AZEL_MNT
**       positions at elevations less than about 14.5 degrees.
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
**   14  All combinations of target and tracking frames produce
**       drift-rates.  These are generated by differencing the
**       prediction for the given date and a few seconds in the future.
**
**  P.T.Wallace   20 March 1999
**
**  Copyright RAL 1999.  All rights reserved.
*/

/* Length of time used for drift-rate calculation (seconds) */
#define TDELTA 10.0

/* Same thing in radians of sidereal time */
#define STDELTA (TDELTA*1.002737909350795*S2R)

/* Elevation (about 14.5 deg) below which refraction disabled, and sine */
#define ELIM 0.253072742
#define SELIM 0.250380004

/* Minimum r-squared for collimation calculations */
#define TINY 1e-20

{
/* Internal frame identifier */
   typedef enum { MNT    = 0,   /* Mount Az/El */
                  TOPO   = 1,   /* Topocentric Az/El */
                  GAPP   = 2,   /* Geocentric apparent RA/Dec */
                  Jxxxx  = 3,   /* FK5 RA/Dec, any equinox */
                  Bxxxx  = 4,   /* FK4 RA/Dec, any equinox */
/* Frames not available as final go below this line */
                  J2000b = 5,   /* FK5 RA/Dec, J2000, barycentric */
                  J2000  = 6,   /* FK5 RA/Dec, J2000, geocentric */
                  OBS    = 7    /* Observed Az/El */
                } COCO_FRAME;
/* n.b. The Jxxxx and Bxxxx frames are barycentric for the      */
/* input coordinates and geocentric for the output coordinates. */

/* Number of different start and end frames */
#define NSTARTF 8
#define NENDF 5

/*
** Direct Transformations
** ----------------------
**
** These are the atomic operations which are carried out in the
** required sequence to accomplish the specified transformation.
**
**    J2000         <-+  --+  --+       --+  <-+
**                    |    |    |         |    |
**    J2000b   <-+  --+    |    |  <-+    |    |
**               |         |    |    |    |    |
**    Bxxxx    --+       <-+    |    |    |    |
**                              |    |    |    |
**    Jxxxx                   <-+  --+    |    |
**                                        |    |
**    GAPP     --+  <-+  --+  <-+       <-+  --+
**               |    |    |    |
**    TOPO     <-+  --+    *    *
**                         |    |
**    OBS      --+  <-+  <-+  --+
**               |    |
**    MNT      <-+  --+
**
** The two transformations marked * have refraction corrections disabled
** at elevations below the Gemini limit, to avoid unnecessary computation.
**
*/
   typedef enum { Bxxxx_to_J2000b,
                  Jxxxx_to_J2000b,
                  J2000b_to_J2000,
                  J2000_to_Bxxxx,
                  J2000_to_Jxxxx,
                  J2000_to_GAPP, GAPP_to_J2000,
                  GAPP_to_TOPO, TOPO_to_GAPP,
                  GAPP_to_OBS, OBS_to_GAPP,
                  OBS_to_MNT, MNT_to_OBS,
                  IDENTITY
                } COCO_TRANS;

/*
** State table defining transformation paths
** -----------------------------------------
**
** For any combination of current frame and final frame, the table
** shows which of the above atomic operations to carry out next and
** what frame we will be in afterwards.
**
*/
   static struct { COCO_TRANS next_trans;
                   COCO_FRAME new_frame;
                 } transform [ NSTARTF ]     /* from */
                             [ NENDF ] = {   /* to   */
     {
       { IDENTITY, MNT  },             /* MNT to MNT      */   /*terminal*/
       { MNT_to_OBS, OBS },            /* MNT to TOPO     */
       { MNT_to_OBS, OBS },            /* MNT to GAPP     */
       { MNT_to_OBS, OBS },            /* MNT to Jxxxx    */
       { MNT_to_OBS, OBS }             /* MNT to Bxxxx    */
     },
     {
       { TOPO_to_GAPP, GAPP },         /* TOPO to MNT     */
       { IDENTITY, TOPO },             /* TOPO to TOPO    */   /*terminal*/
       { TOPO_to_GAPP, GAPP },         /* TOPO to GAPP    */   /*terminal*/
       { TOPO_to_GAPP, GAPP },         /* TOPO to Jxxxx   */
       { TOPO_to_GAPP, GAPP }          /* TOPO to Bxxxx   */
     },
     {
       { GAPP_to_OBS, OBS },           /* GAPP to MNT     */
       { GAPP_to_TOPO, TOPO },         /* GAPP to TOPO    */   /*terminal*/
       { IDENTITY, GAPP },             /* GAPP to GAPP    */   /*terminal*/
       { GAPP_to_J2000, J2000 },       /* GAPP to Jxxxx   */
       { GAPP_to_J2000, J2000 }        /* GAPP to Bxxxx   */
     },
     {
       { Jxxxx_to_J2000b, J2000b },    /* Jxxxx to MNT    */
       { Jxxxx_to_J2000b, J2000b },    /* Jxxxx to TOPO   */
       { Jxxxx_to_J2000b, J2000b },    /* Jxxxx to GAPP   */
       { Jxxxx_to_J2000b, J2000b },    /* Jxxxx to Jxxxx  */
       { Jxxxx_to_J2000b, J2000b }     /* Jxxxx to Bxxxx  */
     },
     {
       { Bxxxx_to_J2000b, J2000b },    /* Bxxxx to MNT    */
       { Bxxxx_to_J2000b, J2000b },    /* Bxxxx to TOPO   */
       { Bxxxx_to_J2000b, J2000b },    /* Bxxxx to GAPP   */
       { Bxxxx_to_J2000b, J2000b },    /* Bxxxx to Jxxxx  */
       { Bxxxx_to_J2000b, J2000b }     /* Bxxxx to Bxxxx  */
     },
     {
       { J2000b_to_J2000, J2000 },     /* J2000b to MNT   */
       { J2000b_to_J2000, J2000 },     /* J2000b to TOPO  */
       { J2000b_to_J2000, J2000 },     /* J2000b to GAPP  */
       { J2000b_to_J2000, J2000 },     /* J2000b to Jxxxx */
       { J2000b_to_J2000, J2000 }      /* J2000b to Bxxxx */
     },
     {
       { J2000_to_GAPP, GAPP },        /* J2000 to MNT    */
       { J2000_to_GAPP, GAPP },        /* J2000 to TOPO   */
       { J2000_to_GAPP, GAPP },        /* J2000 to GAPP   */   /*terminal*/
       { J2000_to_Jxxxx, Jxxxx },      /* J2000 to Jxxxx  */   /*terminal*/
       { J2000_to_Bxxxx, Bxxxx }       /* J2000 to Bxxxx  */   /*terminal*/
     },
     {
       { OBS_to_MNT, MNT },            /* OBS to MNT      */   /*terminal*/
       { OBS_to_GAPP, GAPP },          /* OBS to TOPO     */
       { OBS_to_GAPP, GAPP },          /* OBS to GAPP     */   /*terminal*/
       { OBS_to_GAPP, GAPP },          /* OBS to Jxxxx    */
       { OBS_to_GAPP, GAPP }           /* OBS to Bxxxx    */
     }
   };

/* Target coordinates at epoch and TDELTA later on */
   double acoord [ 2 ], bcoord [ 2 ];

/* Current frame and ending frame of target coordinates */
   COCO_FRAME cframe, eframe;

/* Epoch of observation (J or B year) */
   struct EPOCH epoch2;

/* Vectors and matrices */
   double v1 [ 3 ], v2 [ 3 ], v3 [ 3 ], v4 [ 3 ], rmat [ 3 ] [ 3 ];

/* Azimuth, zenith distance, hour angle, declination, right ascension */
   double az, zd, ha, dc, ra;

/*
**  Star-independent mean-to-apparent parameters:
**
**    [0]       time interval for proper motion (Julian years)
**    [1-3]     barycentric position of the Earth (AU)
**    [4-6]     heliocentric direction of the Earth (unit vector)
**    [7]       (grav rad Sun)*2/(Sun-Earth distance)
**    [8-10]    abv: barycentric Earth velocity in units of c
**    [11]      sqrt(1-v**2) where v=modulus(abv)
**    [12-20]   precession/nutation (3,3) matrix
*/
   static double amprms [ 21 ];

/* Second-epoch appt-to-observed parameters (15th element not required) */
   double aoprms2 [ 14 ];

   int i, j, istat;
   double w1, w2, w3, w4, san, can, saw, caw, sansaw, cansaw,
          sancaw, cancaw, daz, pa, spa, cpa, xr, yr,
          xt, yt, zt, xm, ym, zm, rm2, rm, eta, xi, xi2, eta2p1, sdf,
          r2, r, f, xc, yc, zc;



/*---------------------------------------------------------------------------*/

/*
** -------------
** PRELIMINARIES
** -------------
*/

/* Present status to OK. */
   istat = 0;

/* Express date of observation as a Julian epoch. */
   epoch2.type = 'J';
   epoch2.year = slaEpj ( date);

/* Identify starting frame and adjust equinoxes/epochs to match. */
   switch ( cosys1 ) {

   case AZEL_MNT:
      cframe = MNT;
      break;

   case AZEL_TOPO:
      cframe = TOPO;
      break;

   case APPT:
      cframe = GAPP;
      break;

   case FK5:
      equinox1.year = slaEpco ( 'J', equinox1.type, equinox1.year );
      epoch1.year = slaEpco ( 'J', epoch1.type, epoch1.year );
      cframe = Jxxxx;
      break;

   case FK4:
      equinox1.year = slaEpco ( 'B', equinox1.type, equinox1.year );
      epoch1.year = slaEpco ( 'B', epoch1.type, epoch1.year );
      cframe = Bxxxx;
      break;

   default:
      return -1;
   }

/* Identify ending frame. */
   switch ( cosys2 ) {

   case AZEL_MNT:
      eframe = MNT;
      break;

   case AZEL_TOPO:
      eframe = TOPO;
      break;

   case APPT:
      eframe = GAPP;
      break;

   case FK5:
      equinox2.year = slaEpco ( 'J', equinox2.type, equinox2.year );
      eframe = Jxxxx;
      break; 

   case FK4:
      equinox2.year = slaEpco ( 'B', equinox2.type, equinox2.year );
      eframe = Bxxxx;
      break; 

   default:
      return -2;
   }

/*
** ------------
** SPACE MOTION
** ------------
*/

/* Proceed only if the given position is a mean RA/Dec. */
   if ( cosys1 == FK4 || cosys1 == FK5 ) {

   /* PM specified? */
      if ( pmotion.pm ) {

      /* Yes: apply it and update the epoch. */
         slaPm ( a1, b1, pmotion.pmRA, pmotion.pmDec,
                       pmotion.px, pmotion.rv, epoch1.year, epoch2.year,
                 &a1, &b1 );
         epoch1.year = slaEpco ( epoch1.type, epoch2.type, epoch2.year );

      } else {

      /* No: safe the parallax and radial velocity. */
         pmotion.px = pmotion.rv = 0.0;
      }
   }

/* Copy the starting coordinates into arrays. */
   acoord [ 0 ] = a1;
   bcoord [ 0 ] = b1;

/* Determine the coordinates for epoch TDELTA seconds after date. */
   acoord [ 1 ] = a1 + adot1 * TDELTA;
   bcoord [ 1 ] = b1 + bdot1 * TDELTA;

/*
** --------------
** TRANSFORMATION
** --------------
**
** The following code repeatedly transforms the two positions
** (acoord[0],bcoord[0] for date and acoord[1],bcoord[1] for
** date + TDELTA sec) into the next frame in the chain, until the
** tracking frame is reached.
*/

   do {
      switch ( transform[cframe][eframe].next_trans ) {

      case Bxxxx_to_J2000b:

      /* Transform to 1950 if necessary. */
         if ( equinox1.year != 1950.0 ) {

         /* Generate precession matrix. */
            slaPrebn ( equinox1.year, 1950.0, rmat );

         /* Generate E-terms vector for starting equinox and for 1950. */
            slaEtrms ( equinox1.year, v1 );
            slaEtrms ( 1950.0, v2 );

         /* Transform the two positions to B1950. */
            for ( i = 0; i < 2; i++ ) {

            /* To Cartesian. */
               slaDcs2c ( acoord [ i ], bcoord [ i ], v3 );

            /* Remove E-terms. */
               for ( j = 0; j < 3; j++ ) {
                  v3 [ j ] -= v1 [ j ];
               }

            /* Precess. */
               slaDmxv ( rmat, v3, v4 );

            /* Add E-terms. */
               for ( j = 0; j < 3; j++ ) {
                  v4 [ j ] += v2 [ j ];
               }

            /* To spherical. */
               slaDcc2s ( v4, acoord + i, bcoord + i );
            }
         }

      /* Transform the two positions to J2000. */
         for ( i = 0; i < 2; i++ ) {
            slaFk45z ( acoord [ i ], bcoord [ i ], epoch1.year,
                       &w1, &w2 );
            acoord [ i ] = w1;
            bcoord [ i ] = w2;
         }

         break;

      case Jxxxx_to_J2000b:

      /* Precess to 2000 if necessary. */
         if ( equinox1.year != 2000.0 ) {

         /* Generate the precession matrix. */
            slaPrec ( equinox1.year, 2000.0, rmat );

         /* Precess the two positions to J2000. */
            for ( i = 0; i < 2; i++ ) {
               slaDcs2c ( acoord [ i ], bcoord [ i ], v1 );
               slaDmxv ( rmat, v1, v2 );
               slaDcc2s ( v2, acoord + i, bcoord + i );
            }
         }

         break;

      case J2000b_to_J2000:

      /* Allow for any parallax. */
         if ( pmotion.px > 0.0 ) {

         /* Compute J2000-to-apparent data. */
            slaMappa ( 2000.0, date, amprms );

         /* Parallax in radians. */
            w1 = pmotion.px * AS2R;

         /* Target x,y,z in the J2000 barycentric frame. */
            slaDcs2c ( acoord [ 0 ], bcoord [ 0 ], v1 );

         /* Geocentric direction. */
            for ( i = 0; i < 3; i++ ) {
               v1 [ i ] -= w1 * amprms [ i + 1 ];
            }

         /* Back to spherical. */
            slaDcc2s ( v1, &w1, &w2 );

         /* Apply the adjustment to both sets of coordinates. */
            acoord [ 1 ] += w1 - acoord [ 0 ];
            bcoord [ 1 ] += w2 - bcoord [ 0 ];
            acoord [ 0 ] = w1;
            bcoord [ 0 ] = w2;
         }

         break;

      case J2000_to_Bxxxx:

      /* Transform both positions to B1950. */
         for ( i = 0; i < 2; i++ ) {
            slaFk54z ( acoord [ i ] , bcoord [ i ], epoch2.year,
                       &w1, &w2, &w3, &w4 );
            acoord [ i ] = w1;
            bcoord [ i ] = w2;
         }

      /* If necessary, precess to final equinox. */
         if ( equinox2.year != 1950.0 ) {

         /* Generate the precession matrix. */
            slaPrebn ( 1950.0, equinox2.year, rmat );

         /* Transform the two positions. */
            for ( i = 0; i < 2; i++ ) {
               slaDcs2c ( acoord [ i ], bcoord [ i ], v1 );
               slaDmxv ( rmat, v1, v2 );
               slaDcc2s ( v2, acoord + i, bcoord + i );
            }
         }

         break;

      case J2000_to_Jxxxx:

      /* If necessary, precess to final equinox. */
         if ( equinox2.year != 2000.0 ) {

         /* Generate the precession matrix. */
            slaPrec ( 2000.0, equinox2.year, rmat );

         /* Precess the two positions. */
            for ( i = 0; i < 2; i++ ) {
               slaDcs2c ( acoord [ i ], bcoord [ i ], v1 );
               slaDmxv ( rmat, v1, v2 );
               slaDcc2s ( v2, acoord + i, bcoord + i );
            }
         }

         break;

      case J2000_to_GAPP:

      /* Compute J2000-to-apparent data. */
         slaMappa ( 2000.0, date, amprms );

      /* Transform the two positions. */
         for ( i = 0; i < 2; i++ ) {
            slaMapqkz ( acoord [ i ], bcoord [ i ], amprms, &w1, &w2 );
            acoord [ i ] = w1;
            bcoord [ i ] = w2;
         }

         break;

      case GAPP_to_J2000:

      /* Compute J2000-to-apparent data. */
         slaMappa ( 2000.0, date, amprms );

      /* Transform the two positions. */
         for ( i = 0; i < 2; i++ ) {
            slaAmpqk ( acoord [ i ], bcoord [ i ], amprms, &w1, &w2 );
            acoord [ i ] = w1;
            bcoord [ i ] = w2;
         }

         break;

      case GAPP_to_TOPO:

      /* Save the refraction parameters. */
         w1 = aoprms [ 6 ];
         w2 = aoprms [ 10 ];
         w3 = aoprms [ 11 ];

      /* Eliminate refraction to make topocentric=observed. */
         aoprms [ 6 ] = 0.0;
         aoprms [ 10 ] = 0.0;
         aoprms [ 11 ] = 0.0;

      /* Generate the second-epoch AOP parameters. */
         for ( i = 0; i < 13; i++ ) {
            aoprms2 [ i ] = aoprms [ i ];
         }
         aoprms2 [ 13 ] = aoprms [ 13 ] + STDELTA;

      /* Azimuth correction due to polar motion. */
         daz = aoprms [ 14 ];

      /* Transform the two positions. */
         slaAopqk ( acoord [ 0 ], bcoord [ 0 ], aoprms,
                    &az, &zd, &ha, &dc, &ra );
         acoord [ 0 ] = az + daz;
         bcoord [ 0 ] = D90 - zd;
         slaAopqk ( acoord [ 1 ], bcoord [ 1 ], aoprms2,
                    &az, &zd, &ha, &dc, &ra );
         acoord [ 1 ] = az + daz;
         bcoord [ 1 ] = D90 - zd;

      /* Restore the refraction parameters. */
         aoprms [ 6 ] = w1;
         aoprms [ 10 ] = w2;
         aoprms [ 11 ] = w3;

         break;

      case TOPO_to_GAPP:

      /* Save the refraction parameters. */
         w1 = aoprms [ 6 ];
         w2 = aoprms [ 10 ];
         w3 = aoprms [ 11 ];

      /* Eliminate refraction to make topocentric=observed. */
         aoprms [ 6 ] = 0.0;
         aoprms [ 10 ] = 0.0;
         aoprms [ 11 ] = 0.0;

      /* Compute the second-epoch AOP parameters. */
         for ( i = 0; i < 13; i++ ) {
            aoprms2 [ i ] = aoprms [ i ];
         }
         aoprms2 [ 13 ] = aoprms [ 13 ] + STDELTA;

      /* Azimuth correction due to polar motion. */
         daz = aoprms [ 14 ];

      /* Transform the two positions. */
         slaOapqk ( "A", acoord [ 0 ] - daz, D90 - bcoord [ 0 ], aoprms,
                    acoord, bcoord );
         slaOapqk ( "A", acoord [ 1 ] - daz, D90 - bcoord [ 1 ], aoprms2,
                    acoord + 1, bcoord + 1 );

      /* Restore the refraction parameters. */
         aoprms [ 6 ] = w1;
         aoprms [ 10 ] = w2;
         aoprms [ 11 ] = w3;

         break;

      case GAPP_to_OBS:

      /* Save the refraction parameters. */
         w1 = aoprms [ 6 ];
         w2 = aoprms [ 10 ];
         w3 = aoprms [ 11 ];

      /* If below the Gemini elevation limit, disable refraction. */
         if ( ( sin ( bcoord [ 0 ] ) * aoprms [ 1 ] +
                cos ( bcoord [ 0 ] ) * aoprms [ 2 ] *
                       cos ( aoprms [ 13 ] - acoord [ 0 ] ) ) < SELIM ) {
            aoprms [ 6 ] = 0.0;
            aoprms [ 10 ] = 0.0;
            aoprms [ 11 ] = 0.0;
         }

      /* Generate the second-epoch AOP parameters. */
         for ( i = 0; i < 13; i++ ) {
            aoprms2 [ i ] = aoprms [ i ];
         }
         aoprms2 [ 13 ] = aoprms [ 13 ] + STDELTA;

      /* Azimuth correction due to polar motion. */
         daz = aoprms [ 14 ];

      /* Transform the two positions. */
         slaAopqk ( acoord [ 0 ], bcoord [ 0 ], aoprms,
                    &az, &zd, &ha, &dc, &ra );
         acoord [ 0 ] = az + daz;
         bcoord [ 0 ] = D90 - zd;
         slaAopqk ( acoord [ 1 ], bcoord [ 1 ], aoprms2,
                    &az, &zd, &ha, &dc, &ra );
         acoord [ 1 ] = az + daz;
         bcoord [ 1 ] = D90 - zd;

      /* Ensure the refraction parameters are as they were on entry. */
         aoprms [ 6 ] = w1;
         aoprms [ 10 ] = w2;
         aoprms [ 11 ] = w3;

         break;

      case OBS_to_GAPP:

      /* Save the refraction parameters. */
         w1 = aoprms [ 6 ];
         w2 = aoprms [ 10 ];
         w3 = aoprms [ 11 ];

      /* If below the Gemini elevation limit, disable refraction. */
         if ( bcoord [ 0 ] < ELIM ) {
            aoprms [ 6 ] = 0.0;
            aoprms [ 10 ] = 0.0;
            aoprms [ 11 ] = 0.0;
         }

      /* Compute the second-epoch AOP parameters. */
         for ( i = 0; i < 13; i++ ) {
            aoprms2 [ i ] = aoprms [ i ];
         }
         aoprms2 [ 13 ] = aoprms [ 13 ] + STDELTA;

      /* Azimuth correction due to polar motion. */
         daz = aoprms [ 14 ];

      /* Transform the two positions. */
         slaOapqk ( "A", acoord [ 0 ] - daz, D90 - bcoord [ 0 ], aoprms,
                    acoord, bcoord );
         slaOapqk ( "A", acoord [ 1 ] - daz, D90 - bcoord [ 1 ], aoprms2,
                    acoord + 1, bcoord + 1 );

      /* Ensure the refraction parameters are as they were on entry. */
         aoprms [ 6 ] = w1;
         aoprms [ 10 ] = w2;
         aoprms [ 11 ] = w3;

         break;

      case OBS_to_MNT:

      /* Functions of azimuth axis tilt north and west. */
         san = sin ( tel.an );
         can = cos ( tel.an );
         saw = sin ( tel.aw );
         caw = cos ( tel.aw );
         sansaw = san * saw;
         cansaw = can * saw;
         sancaw = san * caw;
         cancaw = can * caw;

      /* Pointing origin position wrt rotator axis, scaled to radians. */
         xr = tel.pox / tel.fl;
         yr = tel.poy / tel.fl;

      /* In the presence of distortion, with respect to */
      /* gnomonic-projection geometry, xr and yr would  */
      /* at this point require adjustment.              */

      /* Transform the two positions. */
         for ( i = 0; i < 2; i++ ) {

         /* Az (N=0,E=90) and El to RH Cartesian. */
            w2 = cos ( bcoord [ i ] );
            xt = - cos ( acoord [ i ] ) * w2;
            yt = sin ( acoord [ i ] ) * w2;
            zt = sin ( bcoord [ i ] );

         /* Tilt of the azimuth axis. */
            xm = can * xt + san * zt;
            ym = - sansaw * xt + caw * yt + cansaw * zt;
            zm = - sancaw * xt - saw * yt + cancaw * zt;

         /* Cos(El). */
            rm2 = xm * xm + ym * ym;
            if ( rm2 < TINY ) rm2 = TINY;
            rm = sqrt ( rm2 );

         /* Rotator PA, corrected for Az/El nonperpendicularity. */
            if ( !i ) {
               pa = tel.rma + tel.pnpae * rm;
               spa = sin ( pa );
               cpa = cos ( pa );
            }

         /* Pointing origin position wrt ideal axis. */
            eta = xr * spa - yr * cpa + tel.ce;
            xi = - xr * cpa - yr * spa
                 + tel.ca
                 + tel.pnpae * ( zm - eta * rm );

         /* Apply complete collimation correction. */
            xi2 = xi * xi;
            eta2p1 = eta * eta + 1.0;
            sdf = zm * sqrt ( xi2 + eta2p1 );
            r2 = rm2 * eta2p1 - zm * zm * xi2;
            if ( r2 < 0.0 ) {
               istat = 1;
               r2 = 0.0;
            }
            if ( r2 < TINY ) r2 = TINY;
            r = sqrt ( r2 );
            f = ( sdf * eta + r ) / ( eta2p1 * rm * sqrt ( r2 + xi2 ) );
            xc = f * ( xm * r + ym * xi );
            yc = f * ( ym * r - xm * xi );
            zc = ( sdf - eta * r ) / eta2p1;

         /* To spherical (N=0, E=90). */
            acoord [ i ] = slaDranrm ( ( xc == 0.0 && yc == 0.0 ) ?
                                       0.0 : atan2 ( yc, - xc ) );
            bcoord [ i ] = atan2 ( zc, sqrt ( xc * xc + yc * yc ) );
         }

         break;

      case MNT_to_OBS:

      /* Functions of azimuth axis tilt north and west. */
         san = sin ( tel.an );
         can = cos ( tel.an );
         saw = sin ( tel.aw );
         caw = cos ( tel.aw );
         sansaw = san * saw;
         cansaw = can * saw;
         sancaw = san * caw;
         cancaw = can * caw;

      /* Pointing origin position wrt rotator axis. */
         xr = tel.pox / tel.fl;
         yr = tel.poy / tel.fl;

      /* In the presence of distortion, with respect to */
      /* gnomonic-projection geometry, xr and yr would  */
      /* at this point require adjustment.              */

      /* Transform the two positions. */
         for ( i = 0; i < 2; i++ ) {

         /* Az (N=0,E=90) and El to RH Cartesian. */
            w1 = acoord [ i ];
            w2 = cos ( bcoord [ i ] );
            xc = - cos ( w1 ) * w2;
            yc = sin ( w1 ) * w2;
            zc = sin ( bcoord [ i ] );

         /* Rotator PA, corrected for Az/El nonperpendicularity. */
            r2 = xc * xc + yc * yc;
            if ( r2 < TINY ) r2 = TINY;
            r = sqrt ( r2 );
            if ( !i ) {
               pa = tel.rma + tel.pnpae * r;
               spa = sin ( pa );
               cpa = cos ( pa );
            }

         /* Pointing axis position wrt ideal axis. */
            xi = - xr * cpa - yr * spa + tel.ca + tel.pnpae * zc;
            eta = xr * spa - yr * cpa + tel.ce;

         /* Remove combined collimation. */
            f = sqrt ( 1.0 + xi * xi + eta * eta );
            xm = ( xc - ( xi * yc + eta * xc * zc ) / r ) / f;
            ym = ( yc + ( xi * xc - eta * yc * zc ) / r ) / f;
            zm = ( zc + eta * r ) / f;

         /* Remove tilt of the azimuth axis. */
            xt = can * xm - sansaw * ym - sancaw * zm;
            yt = caw * ym - saw * zm;
            zt = san * xm + cansaw * ym + cancaw * zm;

         /* To spherical (N=0, E=90). */
            acoord [ i ] = slaDranrm (
              ( ( xt == 0.0 && yt == 0.0 ) ?  0.0 : atan2 ( yt, - xt ) ) );
            bcoord [ i ] = atan2 ( zt, sqrt ( xt * xt + yt * yt ) );
         }

         break;

      default:

         break;
      }

   /* Update the frame. */
      cframe = transform[cframe][eframe].new_frame;

   } while ( cframe != eframe );

/* Final position. */
   *a2 = slaDranrm ( acoord [ 0 ] );
   *b2 = bcoord [ 0 ];

/* Final drift-rate. */
   *adot2 = slaDrange ( acoord [ 1 ] - acoord [ 0 ] ) / TDELTA;
   *bdot2 = slaDrange ( bcoord [ 1 ] - bcoord [ 0 ] ) / TDELTA;

   return istat;
}
