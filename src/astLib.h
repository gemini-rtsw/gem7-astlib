#ifndef ASTLIBHDEF
#define ASTLIBHDEF

/*
**  - - - - - - - - -
**   a s t L i b . h
**  - - - - - - - - -
**
**  TCS coordinate transformation library macros and
**  prototype function declarations.
**
**  P.T.Wallace   12 November 1998
**
**  Copyright RAL 1998.  All rights reserved.
*/

/* Size of ast context array */
#define AST_CTXA_SIZE 45 

/* Tracking coordinate system (i.e. cosys values) */
typedef enum { AZEL_MNT = 0,        /* 0  Mount Az/El */
               AZEL_TOPO = 1,       /* 1  Topocentric Az/El */
               APPT = 2,            /* 2  Geocentric apparent RA/Dec */
               FK5 = 3,             /* 3  IAU 1976 RA/Dec, any equinox */
               FK4 = 4              /* 4  Pre IAU 1976 RA/Dec, any equinox */
             } FRAMETYPE ;

/* Epoch */
struct EPOCH {
      double year ;        /* Epoch in years */
      char   type ;        /* type of Epoch ('B' or 'J') */
} ;

/* Proper motion */
struct PMPXRV {
      int    pm ;      /* False = proper motion inertially-zero */
      double pmRA ;    /* RA proper motion (radians/year) */
      double pmDec ;   /* Dec proper motion (radians/year) */
      double px ;      /* Parallax (arcsec) */
      double rv ;      /* Radial velocity (km/s) */
};

/* Telescope parameters */
struct TELP {
      double fl;           /* Focal length (metres) */
      double rma;          /* Rotator orientation (radians) */
      double an;           /* Azimuth axis tilt NS (radians) */
      double aw;           /* Azimuth axis tilt EW (radians) */
      double pnpae;        /* Az/El nonperpendicularity (radians) */
      double ca;           /* LR collimation (radians) */
      double ce;           /* UD collimation (radians) */
      double pox;          /* Pointing-origin x-component (metres) */
      double poy;          /* Pointing-origin y-component (metres) */
} ;

/* WCS context */
   struct WCS_CTX {
       double ab0[2];      /* Pre-flexure Az/El */
       struct TELP tel;    /* Telescope parameters */
       double aoprms[15];  /* Apparent-to-observed parameters */
       double m2xy[3][2];  /* M2 tip/tilt (3 chop states) */
       double time;        /* Gemini raw time */
   };

/* WCS transformation */
   struct WCS {
      double ab0[2];       /* Celestial coordinates at x=y=0 */
      double coeffs[6];    /* Affine transformation coefficients */
   };

/* Pointing origin structure */
   struct PO {
      double mx;           /* Mount x */
      double my;           /* Mount y */
      double ax;           /* Source chop A x */
      double ay;           /* Source chop A y */
      double bx;           /* Source chop B x */
      double by;           /* Source chop B y */
      double cx;           /* Source chop C x */
      double cy;           /* Source chop C y */
      double time;         /* Gemini raw time */
   };

/* Function prototypes */
int astCoco ( double, double, struct PMPXRV,
              FRAMETYPE, struct EPOCH, struct EPOCH,
              FRAMETYPE, struct EPOCH,
              double, double*, struct TELP, double*, double* );
int astCocoR ( double, double, struct PMPXRV, double, double,
               FRAMETYPE, struct EPOCH, struct EPOCH,
               FRAMETYPE, struct EPOCH,
               double, double*, struct TELP, double*, double*,
               double*, double* );
int astCtx2tr ( struct WCS_CTX, FRAMETYPE, struct EPOCH,
                double, int, struct WCS*, double* );
int astFitij ( int, double[][2], double[][2], double[],
               double*, double*, double*, double* );
int astFITSs ( struct WCS, FRAMETYPE, struct EPOCH, double,
               char*, char*, char*, char*, char*, char*, char*,
               char*, char*, char*, char*, char*, char* );
int astFITSv ( struct WCS, FRAMETYPE, struct EPOCH, double,
               char*, double*, double*, char*, double*, double*, double*,
               double*, double*, double*, char*, double*, double* );
int astGetctx ( struct WCS_CTX* );
int astGetdistortion (double[]);
int astGetpo ( struct PO* );
int astGettr ( FRAMETYPE, struct EPOCH, double, int, struct WCS*, double* );
int astInvtr ( struct WCS, struct WCS* );
int astRot ( FRAMETYPE, struct EPOCH, double, double,
             double, double*, struct TELP, double, double, double* );
int astS2xy ( double, double, FRAMETYPE, struct EPOCH,
              double, int, double*, double* );
int astS2xyq ( double, double, struct WCS, double*, double* );
int astSetctx ( double[AST_CTXA_SIZE] );
int astSetdistortion (double[6]);
int astSimctx ( double, double, double, double, double, double, double,
                double, double, double, double, struct TELP, double[3][2],
                double, double, FRAMETYPE, struct EPOCH, struct WCS_CTX* );
int astSimctx_r ( double, double, double, double, double, double, double,
                  double, double, double, double, struct TELP, double[3][2],
                  double, double, FRAMETYPE, struct EPOCH,
                  double, double, FRAMETYPE, struct EPOCH, struct WCS_CTX* );
int astXtndtr ( double*, struct WCS, struct WCS* );
int astXy2s ( double, double, FRAMETYPE, struct EPOCH,
              double, int, double*, double* );
int astXy2sq ( double, double, struct WCS, double*, double* );

#endif
