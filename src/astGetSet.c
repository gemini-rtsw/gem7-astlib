#if defined (vxWorks)
#include <semLib.h>
#endif
#include <slalib.h>
#include "astLib.h"

/*
** Static variables to store the context and the semaphore that protects
** the structure from simultanious access by both astSetctx and astGetctx.
**
** The context is initialised to a set of harmless numbers so that the ast
** routines can be called without a connection to the TCS.
**/

static double Ctxa[AST_CTXA_SIZE] = {
    0.0, 0.0, 0.70710678, 0.70710678,
    128000.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 1.0, 0.0, 0.0, 275.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 1.0, 0.0, 0.0, 0.0, 1.0
    };
#if defined (vxWorks)
static SEM_ID AstSemId;
#endif

int astGetctx ( struct WCS_CTX *ctxp )
/*
**  - - - - - - - - - -
**   a s t G e t c t x
**  - - - - - - - - - -
**
**  Obtain the context information needed for generating world
**  coordinate system transformations.
**
**  RETURNED (argument):
**    ctxp      struct*    WCS context
**
**  RETURNED (function value):
**              int        0 = OK 
**
**  Defined in astLib.h:
**    WCS_CTX   struct     WCS context
**    TELP      struct     telescope parameters
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
**  Called:  semMCreate, semTake, slaDcc2s, semGive
**
**  Notes:
**
**    1  The structure elements WCS_CTX.tel.pox and WCS_CTX.tel.poy
**       receive zero values, rather than, for example, the current
**       mount pointing origin coordinates.  This is because they
**       have no effect on the WCS, though other routines which use
**       the WCS_CTX structure do require them.
**
**    2  The structure element aoprms[12] is not used by the other
**       WCS functions and is set to zero.
**
**  D.L.Terrett   7 May 1998
**
**  Copyright RAL 1998.  All rights reserved.
*/

{
   double v0[3];

/* Take the semaphore (creating it if necessary). */
#if defined (vxWorks)
   if ( !AstSemId ) 
      AstSemId = semMCreate( SEM_Q_PRIORITY |
                             SEM_DELETE_SAFE |
                             SEM_INVERSION_SAFE);
   (void) semTake( AstSemId, WAIT_FOREVER );
#else
   if ( Ctxa[0] == 0.0 ) return -1;
#endif


/* Copy the timestamp. */
   ctxp->time = Ctxa[0];

/* Copy the telescope parameters. */
   ctxp->tel.fl = Ctxa[4];
   ctxp->tel.rma = Ctxa[5];
   ctxp->tel.an = Ctxa[6];
   ctxp->tel.aw = Ctxa[7];
   ctxp->tel.pnpae = Ctxa[8];
   ctxp->tel.ca = Ctxa[9];
   ctxp->tel.ce = Ctxa[10];
   ctxp->tel.pox = 0.0;         /* For safety: no bearing on WCS */
   ctxp->tel.poy = 0.0;         /* For safety: no bearing on WCS */

/* Copy the apparent-to-observed parameters. */
   ctxp->aoprms[0] =  Ctxa[11];
   ctxp->aoprms[1] = Ctxa[12];
   ctxp->aoprms[2] = Ctxa[13];
   ctxp->aoprms[3] = Ctxa[14];
   ctxp->aoprms[4] = Ctxa[15];
   ctxp->aoprms[5] = Ctxa[16];
   ctxp->aoprms[6] = Ctxa[17];
   ctxp->aoprms[7] = Ctxa[18];
   ctxp->aoprms[8] = Ctxa[19];
   ctxp->aoprms[9] = Ctxa[20];
   ctxp->aoprms[10] = Ctxa[21];
   ctxp->aoprms[11] = Ctxa[22];
   ctxp->aoprms[12] = 0.0;      /* For safety: no bearing on WCS */
   ctxp->aoprms[13] = Ctxa[23];
   ctxp->aoprms[14] = Ctxa[24];

/* Copy the m2 tip tilts. */
   ctxp->m2xy[0][0] = Ctxa[25];
   ctxp->m2xy[0][1] = Ctxa[26];
   ctxp->m2xy[1][0] = Ctxa[27];
   ctxp->m2xy[1][1] = Ctxa[28];
   ctxp->m2xy[2][0] = Ctxa[29];
   ctxp->m2xy[2][1] = Ctxa[30];

/* Grab the pre-flexure mount coordinates. */
   v0[0] = -Ctxa[1];
   v0[1] = Ctxa[2];
   v0[2] = Ctxa[3];

/* Release the semaphore. */
#if defined (vxWorks)
   (void) semGive( AstSemId );
#endif

/* Convert the mount pre-TF az/el to spherical coordinates. */
   slaDcc2s ( v0, & ctxp->ab0[0], & ctxp->ab0[1] );

   return 0;

}

int astGetpo ( struct PO *pop )
/*
**  - - - - - - - - -
**   a s t G e t p o
**  - - - - - - - - -
**
**  Obtain the current pointing origins.
**
**  RETURNED (argument):
**    pop       struct*    pointing origins
**
**  RETURNED (function value):
**              int        0 = OK 
**
**  Called:  semMCreate, semTake, semGive
**
**  Defined in astLib.h:
**    PO        struct     pointing origin structure
**
**  D.L.Terrett   7 May 1998
**
**  Copyright RAL 1998.  All rights reserved.
*/

{

/* Take the semaphore (creating it if necessary). */
#if defined (vxWorks)
   if ( !AstSemId ) 
      AstSemId = semMCreate( SEM_Q_PRIORITY |
                             SEM_DELETE_SAFE |
                             SEM_INVERSION_SAFE);
   (void) semTake( AstSemId, WAIT_FOREVER );
#else
   if ( Ctxa[0] == 0.0 ) return -1;
#endif

/* Copy the timestamp. */
   pop->time = Ctxa[0];

/* Copy the pointing origins. */
   pop->mx = Ctxa[31];
   pop->my = Ctxa[32];
   pop->ax = Ctxa[33];
   pop->ay = Ctxa[34];
   pop->bx = Ctxa[35];
   pop->by = Ctxa[36];
   pop->cx = Ctxa[37];
   pop->cy = Ctxa[38];

/* Release the semaphore. */
#if defined (vxWorks)
   (void) semGive( AstSemId );
#endif

   return 0;

}

int astSetctx ( double ctxa[AST_CTXA_SIZE] )
/*
**  - - - - - - - - - -
**   a s t S e t c t x
**  - - - - - - - - - -
**
**  Set the context information needed for generating world
**  coordinate system transformations.
**
**  This routine is intended only to be called by the EPICS code that
**  copied data from the TCS.
**
**  GIVEN (argument):
**    ctxa      double*    array containing the context
**
**  Called:  semMCreate, semTake, semGive
**
**  D.L.Terrett   7 May 1998
**
**  Copyright RAL 1998.  All rights reserved.
*/

{
   int i;

/* Take the semaphore (creating it if necessary). */
#if defined (vxWorks)
   if ( !AstSemId ) 
      AstSemId = semMCreate( SEM_Q_PRIORITY |
                             SEM_DELETE_SAFE |
                             SEM_INVERSION_SAFE);
   (void) semTake( AstSemId, WAIT_FOREVER );
#endif

/* Copy the context structure array. */
   for ( i = 0; i < AST_CTXA_SIZE; i++ ) Ctxa[i] = ctxa[i];

/* Release the semaphore. */
#if defined (vxWorks)
   (void) semGive( AstSemId );
#endif

   return 0;
}

int astGetdistortion ( double distort[6] )
/*
**  - - - - - - - - - - - - - - - -
**   a s t G e t d i s t o r t i o n 
**  - - - - - - - - - - - - - - - -
**
**  Fetch the optical distortion coefficients that describe the
**  departures of the telescope field from a simple tangent
**  plane projection.
**
**  GIVEN (argument):
**    distort      double*   affine coefficients describing distortion 
**
**  Called:  semMCreate, semTake, semGive
**
**  C. J. Mayer  21 October 2002 
**
**  Copyright Observatory Sciences Ltd. 2002  All rights reserved.
*/

{
   int i;

/* Take the semaphore (creating it if necessary). */
#if defined (vxWorks)
   if ( !AstSemId ) 
      AstSemId = semMCreate( SEM_Q_PRIORITY |
                             SEM_DELETE_SAFE |
                             SEM_INVERSION_SAFE);
   (void) semTake( AstSemId, WAIT_FOREVER );
#endif

/* Copy the context structure array. */
   for ( i = 0; i < 6; i++ ) distort[i] = Ctxa[i+39] ;

/* Release the semaphore. */
#if defined (vxWorks)
   (void) semGive( AstSemId );
#endif

   return 0;
}

int astSetdistortion ( double distort[6] )
/*
**  - - - - - - - - - - - - - - - -
**   a s t S e t d i s t o r t i o n 
**  - - - - - - - - - - - - - - - -
**
**  Set the optical distortion coefficients that describe the
**  departures of the telescope field from a simple tangent
**  plane projection.
**
**  This routine would normally only be called when simulating the
**  WCS context. Most systems will simply call astSetctx after fetching 
**  the information from the TCS. If a system is simulating the WCS
**  context and dosn't call this routine then a simple tangent plane
**  projection is used. 
**
**  GIVEN (argument):
**    distort      double*   affine coefficients describing distortion 
**
**  Called:  semMCreate, semTake, semGive
**
**  C. J. Mayer  21 October 2002 
**
**  Copyright Observatory Sciences Ltd. 2002  All rights reserved.
*/

{
   int i;

/* Take the semaphore (creating it if necessary). */
#if defined (vxWorks)
   if ( !AstSemId ) 
      AstSemId = semMCreate( SEM_Q_PRIORITY |
                             SEM_DELETE_SAFE |
                             SEM_INVERSION_SAFE);
   (void) semTake( AstSemId, WAIT_FOREVER );
#endif

/* Copy the context structure array. */
   for ( i = 0; i < 6; i++ ) Ctxa[i+39] = distort[i];

/* Release the semaphore. */
#if defined (vxWorks)
   (void) semGive( AstSemId );
#endif

   return 0;
}
