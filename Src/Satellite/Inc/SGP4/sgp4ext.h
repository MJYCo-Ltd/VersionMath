#ifndef _sgp4ext_
#define _sgp4ext_
/*     ----------------------------------------------------------------
*
*                                 sgp4ext.h
*
*    this file contains extra routines needed for the main test program for sgp4.
*    these routines are derived from the astro libraries.
*
*                            companion code for
*               fundamentals of astrodynamics and applications
*                                    2007
*                              by david vallado
*
*       (w) 719-573-2600, email dvallado@agi.com
*
*    current :
*              20 apr 07  david vallado
*                           misc documentation updates
*    changes :
*              14 aug 06  david vallado
*                           original baseline
*       ----------------------------------------------------------------      */

#include <string.h>
#include <math.h>

#include "sgp4unit.h"


// ------------------------- function declarations -------------------------

double  sgn(double x);

double  angle(double vec1[3],double vec2[3]);

void newtonnu(double ecc, double nu,double& e0, double& m);

double asinh(double xval);

void rv2coe(double r[3], double v[3], double mu,
double& p, double& a, double& ecc, double& incl, double& omega, double& argp,
double& nu, double& m, double& arglat, double& truelon, double& lonper);

#endif

