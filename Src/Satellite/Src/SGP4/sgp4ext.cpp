/*     ----------------------------------------------------------------
*
*                               sgp4ext.cpp
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
*               7 may 08  david vallado
*                           fix sgn
*    changes :
*               2 apr 07  david vallado
*                           fix jday floor and str lengths
*                           updates for constants
*              14 aug 06  david vallado
*                           original baseline
*       ----------------------------------------------------------------      */
#include <SGP4/sgp4ext.h>
#include "sofa.h"

double  sgn
        (
          double x
        )
   {
     if (x < 0.0)
       {
          return -1.0;
       }
       else
       {
          return 1.0;
       }

   }  // end sgn

/* -----------------------------------------------------------------------------
*
*                           procedure angle
*
*  this procedure calculates the angle between two vectors.  the output is
*    set to 999999.1 to indicate an undefined value.  be sure to check for
*    this at the output phase.
*
*  author        : david vallado                  719-573-2600    1 mar 2001
*
*  inputs          description                    range / units
*    vec1        - vector number 1
*    vec2        - vector number 2
*
*  outputs       :
*    theta       - angle between the two vectors  -pi to pi
*
*  locals        :
*    temp        - temporary real variable
*
*  coupling      :
*    dot           dot product of two vectors
* --------------------------------------------------------------------------- */

double  angle
        (
          double vec1[3],
          double vec2[3]
        )
   {
     double dSmall, undefined, magv1, magv2, temp;
     dSmall     = 0.00000001;
     undefined = 999999.1;

     magv1 = iauPm(vec1);
     magv2 = iauPm(vec2);

     if (magv1*magv2 > dSmall*dSmall)
       {
         temp= iauPdp(vec1,vec2) / (magv1*magv2);
         if (fabs( temp ) > 1.0)
             temp= sgn(temp) * 1.0;
         return acos( temp );
       }
       else
         return undefined;
   }  // end angle


/* -----------------------------------------------------------------------------
*
*                           function asinh
*
*  this function evaluates the inverse hyperbolic sine function.
*
*  author        : david vallado                  719-573-2600    1 mar 2001
*
*  inputs          description                    range / units
*    xval        - angle value                                  any real
*
*  outputs       :
*    arcsinh     - result                                       any real
*
*  locals        :
*    none.
*
*  coupling      :
*    none.
*
* --------------------------------------------------------------------------- */

double  asinh
        (
          double xval
        )
   {
     return log( xval + sqrt( xval*xval + 1.0 ) );
   }  // end asinh


/* -----------------------------------------------------------------------------
*
*                           function newtonnu
*
*  this function solves keplers equation when the true anomaly is known.
*    the mean and eccentric, parabolic, or hyperbolic anomaly is also found.
*    the parabolic limit at 168?is arbitrary. the hyperbolic anomaly is also
*    limited. the hyperbolic sine is used because it's not double valued.
*
*  author        : david vallado                  719-573-2600   27 may 2002
*
*  revisions
*    vallado     - fix small                                     24 sep 2002
*
*  inputs          description                    range / units
*    ecc         - eccentricity                   0.0  to
*    nu          - true anomaly                   -2pi to 2pi rad
*
*  outputs       :
*    e0          - eccentric anomaly              0.0  to 2pi rad       153.02 ?
*    m           - mean anomaly                   0.0  to 2pi rad       151.7425 ?
*
*  locals        :
*    e1          - eccentric anomaly, next value  rad
*    sine        - sine of e
*    cose        - cosine of e
*    ktr         - index
*
*  coupling      :
*    asinh       - arc hyperbolic sine
*
*  references    :
*    vallado       2007, 85, alg 5
* --------------------------------------------------------------------------- */

void newtonnu
     (
       double ecc, double nu,
       double& e0, double& m
     )
     {
       double dsmall, sine, cose;

     // ---------------------  implementation   ---------------------
     e0= 999999.9;
     m = 999999.9;
     dsmall = 0.00000001;

     // --------------------------- circular ------------------------
     if ( fabs( ecc ) < dsmall  )
       {
         m = nu;
         e0= nu;
       }
       else
         // ---------------------- elliptical -----------------------
         if ( ecc < 1.0-dsmall  )
           {
             sine= ( sqrt( 1.0 -ecc*ecc ) * sin(nu) ) / ( 1.0 +ecc*cos(nu) );
             cose= ( ecc + cos(nu) ) / ( 1.0  + ecc*cos(nu) );
             e0  = atan2( sine,cose );
             m   = e0 - ecc*sin(e0);
           }
           else
             // -------------------- hyperbolic  --------------------
             if ( ecc > 1.0 + dsmall  )
               {
                 if ((ecc > 1.0 ) && (fabs(nu)+0.00001 < DPI-acos(1.0 /ecc)))
                   {
                     sine= ( sqrt( ecc*ecc-1.0  ) * sin(nu) ) / ( 1.0  + ecc*cos(nu) );
                     e0  = asinh( sine );
                     m   = ecc*sinh(e0) - e0;
                   }
                }
               else
                 // ----------------- parabolic ---------------------
                 if ( fabs(nu) < 168.0*DD2R  )
                   {
                     e0= tan( nu*0.5  );
                     m = e0 + (e0*e0*e0)/3.0;
                   }

     if ( ecc < 1.0  )
       {
         m = fmod( m,D2PI );
         if ( m < 0.0  )
             m = m + D2PI;
         e0 = fmod( e0,D2PI );
       }
   }  // end newtonnu


/* -----------------------------------------------------------------------------
*
*                           function rv2coe
*
*  this function finds the classical orbital elements given the geocentric
*    equatorial position and velocity vectors.
*
*  author        : david vallado                  719-573-2600   21 jun 2002
*
*  revisions
*    vallado     - fix special cases                              5 sep 2002
*    vallado     - delete extra check in inclination code        16 oct 2002
*    vallado     - add constant file use                         29 jun 2003
*    vallado     - add mu                                         2 apr 2007
*
*  inputs          description                    range / units
*    r           - ijk position vector            km
*    v           - ijk velocity vector            km / s
*    mu          - gravitational parameter        km3 / s2
*
*  outputs       :
*    p           - semilatus rectum               km
*    a           - semimajor axis                 km
*    ecc         - eccentricity
*    incl        - inclination                    0.0  to pi rad
*    omega       - longitude of ascending node    0.0  to 2pi rad
*    argp        - argument of perigee            0.0  to 2pi rad
*    nu          - true anomaly                   0.0  to 2pi rad
*    m           - mean anomaly                   0.0  to 2pi rad
*    arglat      - argument of latitude      (ci) 0.0  to 2pi rad
*    truelon     - true longitude            (ce) 0.0  to 2pi rad
*    lonper      - longitude of periapsis    (ee) 0.0  to 2pi rad
*
*  locals        :
*    hbar        - angular momentum h vector      km2 / s
*    ebar        - eccentricity     e vector
*    nbar        - line of nodes    n vector
*    c1          - v**2 - u/r
*    rdotv       - r dot v
*    hk          - hk unit vector
*    sme         - specfic mechanical energy      km2 / s2
*    i           - index
*    e           - eccentric, parabolic,
*                  hyperbolic anomaly             rad
*    temp        - temporary variable
*    typeorbit   - type of orbit                  ee, ei, ce, ci
*
*  coupling      :
*    mag         - magnitude of a vector
*    cross       - cross product of two vectors
*    angle       - find the angle between two vectors
*    newtonnu    - find the mean anomaly
*
*  references    :
*    vallado       2007, 126, alg 9, ex 2-5
* --------------------------------------------------------------------------- */

void rv2coe
     (
       double r[3], double v[3], double mu,
       double& p, double& a, double& ecc, double& incl, double& omega, double& argp,
       double& nu, double& m, double& arglat, double& truelon, double& lonper
     )
     {
       double undefined, dsmall, hbar[3], nbar[3], magr, magv, magn, ebar[3], sme,
              rdotv, infinite, temp, c1, hk, twopi, magh, halfpi, e;

       int i;
       char typeorbit[3];

     twopi  = D2PI;
     halfpi = 0.5 * DPI;
     dsmall  = 0.00000001;
     undefined = 999999.1;
     infinite  = 999999.9;

     // -------------------------  implementation   -----------------
     magr = iauPm( r );
     magv = iauPm( v );

     // ------------------  find h n and e vectors   ----------------
     iauPxp(r,v,hbar);
     magh = iauPm( hbar );
     if ( magh > dsmall )
       {
         nbar[0]= -hbar[1];
         nbar[1]=  hbar[0];
         nbar[2]=   0.0;
         magn = iauPm( nbar );
         c1 = magv*magv - mu /magr;
         rdotv = iauPdp( r,v );
         for (i= 0; i <= 2; ++i)
             ebar[i]= (c1*r[i] - rdotv*v[i])/mu;
         ecc = iauPm( ebar );

         // ------------  find a e and semi-latus rectum   ----------
         sme= ( magv*magv*0.5  ) - ( mu /magr );
         if ( fabs( sme ) > dsmall )
             a= -mu  / (2.0 *sme);
           else
             a= infinite;
         p = magh*magh/mu;

         // -----------------  find inclination   -------------------
         hk= hbar[2]/magh;
         incl= acos( hk );

         // --------  determine type of orbit for later use  --------
         // ------ elliptical, parabolic, hyperbolic inclined -------
         strcpy(typeorbit,"ei");
         if ( ecc < dsmall )
         {
             // ----------------  circular equatorial ---------------
             if  ((incl<dsmall) | (fabs(incl-DPI)<dsmall))
             {
                 strcpy(typeorbit,"ce");
             }
             else
             {
                 // --------------  circular inclined ---------------
                 strcpy(typeorbit,"ci");
             }
         }
         else
         {
             // - elliptical, parabolic, hyperbolic equatorial --
             if  ((incl<dsmall) | (fabs(incl-DPI)<dsmall))
             {
                 strcpy(typeorbit,"ee");
             }
         }

         // ----------  find longitude of ascending node ------------
         if ( magn > dsmall )
           {
             temp= nbar[0] / magn;
             if ( fabs(temp) > 1.0  )
                 temp= sgn(temp);
             omega= acos( temp );
             if ( nbar[1] < 0.0  )
                 omega= twopi - omega;
           }
           else
             omega= undefined;

         // ---------------- find argument of perigee ---------------
         if ( strcmp(typeorbit,"ei") == 0 )
           {
             argp = angle( nbar,ebar);
             if ( ebar[2] < 0.0  )
                 argp= twopi - argp;
           }
           else
             argp= undefined;

         // ------------  find true anomaly at epoch    -------------
         if ( typeorbit[0] == 'e' )
           {
             nu =  angle( ebar,r);
             if ( rdotv < 0.0  )
                 nu= twopi - nu;
           }
           else
             nu= undefined;

         // ----  find argument of latitude - circular inclined -----
         if ( strcmp(typeorbit,"ci") == 0 )
           {
             arglat = angle( nbar,r );
             if ( r[2] < 0.0  )
                 arglat= twopi - arglat;
             m = arglat;
           }
           else
             arglat= undefined;

         // -- find longitude of perigee - elliptical equatorial ----
         if(strcmp(typeorbit,"ee") == 0)
           {
             temp= ebar[0]/ecc;
             if ( fabs(temp) > 1.0  )
                 temp= sgn(temp);
             lonper= acos( temp );
             if ( ebar[1] < 0.0  )
                 lonper= twopi - lonper;
             if ( incl > halfpi )
                 lonper= twopi - lonper;
           }
           else
             lonper= undefined;

         // -------- find true longitude - circular equatorial ------
         if  (( magr>dsmall ) && ( strcmp(typeorbit,"ce") == 0 ))
           {
             temp= r[0]/magr;
             if ( fabs(temp) > 1.0  )
                 temp= sgn(temp);
             truelon= acos( temp );
             if ( r[1] < 0.0  )
                 truelon= twopi - truelon;
             if ( incl > halfpi )
                 truelon= twopi - truelon;
             m = truelon;
           }
           else
             truelon= undefined;

         // ------------ find mean anomaly for all orbits -----------
         if ( typeorbit[0] == 'e' )
         {
             newtonnu(ecc,nu,  e, m);
         }
         else
         {
             ecc = 0.;
         }
     }
      else
     {
        p    = undefined;
        a    = undefined;
        ecc  = undefined;
        incl = undefined;
        omega= undefined;
        argp = undefined;
        nu   = undefined;
        m    = undefined;
        arglat = undefined;
        truelon= undefined;
        lonper = undefined;
     }
   }  // end rv2coe
