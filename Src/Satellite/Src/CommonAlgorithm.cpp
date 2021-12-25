#include <cmath>
#include <cstdio>
#include <limits>
#include <cstring>
#include <Math/VecMat.h>
#include <Math/MathCommonAlgorithm.h>
#include <GisMath/GisMath.h>
#include <Satellite/CoorSys.h>

#include "CommonAlgorithm.h"
using namespace std;
#include "sofa.h"
#include "sofam.h"

using namespace Math;
using namespace Aerospace;
static const double eps_mach(numeric_limits<double>::epsilon());
#define MINUTES_PER_DAY 1440.
#define MINUTES_PER_DAY_SQUARED (MINUTES_PER_DAY * MINUTES_PER_DAY)
#define MINUTES_PER_DAY_CUBED (MINUTES_PER_DAY_SQUARED * MINUTES_PER_DAY)

#define AE 1.0
#define J1900 (2451545.5 - 36525. - 1.)

/// 岁差矩阵
void PrecMatrix(double Mjd_2, double dTemp[3][3])
{
    double zeta, z, theta;
    iauPrec76(DJ00, 0.0, DJM0, Mjd_2,&zeta,&z,&theta);

    iauRz(-zeta,dTemp);
    iauRy(theta,dTemp);
    iauRz(-z,dTemp);

}

/// 章动矩阵
void NutMatrix (double Mjd_TT,double & dpsi, double& eps, double& deps, double& dom,double dTemp[3][3])
{
    // Mean obliquity of the ecliptic

    eps = iauObl80(DJM0,Mjd_TT);

    // Nutation in longitude and obliquity

    iauNut80(DJM0,Mjd_TT,&dpsi,&deps,&dom);

    // Transformation from mean to true equator and equinox

    iauRx(eps,dTemp);
    iauRz(-dpsi,dTemp);
    iauRx(-eps-deps,dTemp);
}

/// 格林尼治真恒星时
double GAST (double Mjd_UT1, double dpsi, double deps, double eps, double om)
{
    return(iauAnp((iauGmst82(DJM0,Mjd_UT1) + EquationOfEquinoxes(dpsi,deps,eps,om))));
}

/// 地球旋转矩阵
void GHAMatrix(double Mjd_UT1, double dpsi, double deps, double eps, double om, double dTemp[3][3])
{
    iauRz(GAST(Mjd_UT1,dpsi,deps,eps,om),dTemp);
}

/// 极移矩阵
void PoleMatrix (const double& dX, const double& dY, double dTemp[3][3])
{
    iauRy(-dX,dTemp);
    iauRx(-dY,dTemp);
}

/// 计算时间
double CalP(double dT)
{
    return (0.0016568*sin((35999.37*dT+357.5)*DD2R)
            +0.0000224*sin((32964.5*dT+246)*DD2R)
            +0.0000138*sin((71998.7*dT+355)*DD2R)
            +0.0000048*sin((3034.9*dT+25)*DD2R)
            +0.0000047*sin((34777.3*dT+230)*DD2R));
}

/// 解开普勒方程
double KeplerFunc2(double ksi,double eta,double lamda)
{
    int JJ = 0;
    double Ew = lamda + ksi*sin(lamda) - eta*cos(lamda);
    double sEw = sin(Ew);
    double cEw = cos(Ew);
    double W0 = Ew - (lamda + ksi*sEw - eta*cEw);
    while(fabs(W0)>1e-8)
    {
        Ew = Ew - W0 / (1.0 - (ksi*cEw + eta*sEw));
        ++JJ;
        if(JJ>20)
        {
            break;
        }
        sEw = sin(Ew);
        cEw = cos(Ew);
        W0 = Ew - (lamda + ksi*sEw - eta*cEw);
    }
    return Ew;
}

/// 查找轨道的短周期项
void Short2(const CVector & MElem,double ZS[6],double & AS2)
{
    double A2 = 0.0016239400253297270;
    double M_a = MElem(0) / 6378245.;
    double M_i = MElem(2);
    //double M_omiga = MElem.o;
    double M_ksi = MElem(1)*cos(MElem(4));
    double M_eta = MElem(1)*sin(MElem(4));
    double M_lamda = fmod(MElem(4) + MElem(5) + (double)D2PI,(double)D2PI);

    double E2 = M_ksi*M_ksi + M_eta*M_eta;  // e^2
    double WE2 = 1.0 - E2;                  // 1 - e^2
    double GE2 = sqrt(WE2);                 // sqrt(1-e^2)
    double F1e = 1.0 / (1.0 + GE2);         // F1e = 1/( 1+sqrt(1-e^2) )
    // 这部分计算u，见《航天器轨道理论》P174下半页
    // 这里定义 ksi = e*cos(w) , eta = e*sin(w) 与刘林书上差个负号，相应的公式均加个负号
    double UE = KeplerFunc2(M_ksi,M_eta,M_lamda); // E + w
    double SU = sin(UE);                         // sin(E+w)
    double CU = cos(UE);                         // cos(E+w)
    double Ef = M_ksi*SU - M_eta*CU;             // UE-lamda = E - M
    double RA = 1.0/(1.0 - M_ksi*CU - M_eta*SU); // a/r
    double RA3 = RA * RA * RA;                   // (a/r)^3
    // 下面U的计算有两种，两种结果都正确,但二种方法在f=90和270时计算存在截断误差
    // 1.刘林《航天器轨道理论》书上的公式：
    double SfE = RA*Ef*( 1.0 + F1e*(-M_eta*SU - M_ksi*CU) );  // sin(u-UE)
    double W = SfE / sqrt( 1.0 - SfE*SfE );      // tan(u-UE)
    double fE = atan(W);                         // u-UE
    double UJMw = fE + Ef;                       // u-lmada = (u-UE) + (UE-lamda)
    double U = M_lamda + UJMw;                   // u
    double S1u = sin(U);           // sin(u)
    double C1u = cos(U);           // cos(u)
    // 2.《航天器轨道确定》书上的公式：
    //double S1u = RA*(M_ksi*(M_lamda-UE)*F1e + SU - M_eta);
    //double C1u = RA*(-M_eta*(M_lamda-UE)*F1e + CU - M_ksi);
    //double U = atan2(S1u,C1u);
    //double UJMw = U - M_lamda; // u-lamda
    /////////////////////////////////////////////////////
    double S2u = 2.0 * S1u * C1u;  // sin(2u)
    double C2u = 2.0*C1u*C1u-1.0;  // cos(2u)
    double SI = sin(M_i);          // sin(i)
    double S2 = SI * SI;           // sin(i)^2
    ZS[0] = A2/M_a * ( (2.0/3.0-S2) * (RA3-pow(WE2,-1.5)) + S2*RA3*C2u );   // a的短周期项as
    //--ZS[0] END----------------------------------------------------
    double P = M_a * WE2;   // p = a(1-e^2)
    double P2 = P * P;      // p^2
    double A2P2 = A2/P2;    // A2/p^2
    double F2e = (1.0+2.0*GE2)*F1e*F1e;     // F2(e)
    double F3e = 0.25*(1.0+3.0*GE2)*F1e*F1e*F1e; // F3(e)
    double S3u = S1u*C2u + C1u*S2u;  // sin(3u)
    double C3u = C1u*C2u - S1u*S2u;  // cos(3u)
    double S4u = S1u*C3u + C1u*S3u;  // sin(4u)
    double C4u = C1u*C3u - S1u*S3u;  // cos(4u)
    double S5u = S1u*C4u + C1u*S4u;  // sin(5u)
    double C5u = C1u*C4u - S1u*S4u;  // cos(5u)
    double CI = cos(M_i);           // cos(i)
    double W23i = 1.0 - 1.5*S2;     // 1-3/2*sin(i)^2
    double W25i = 2.0 - 2.5*S2;     // 2-5/2*sin(i)^2
    double Wkc1 = M_ksi * M_ksi;    // ksi^2
    double Wkc2 = M_eta * M_eta;    // eta^2
    double Wkc3 = Wkc1 - Wkc2;      // ksi^2 - eta^2
    double Wkc4 = M_ksi * M_eta;    // ksi*eta

    ZS[1] = 0.5*A2P2*SI*CI*(M_ksi*C1u-M_eta*S1u+C2u+(M_ksi*C3u+M_eta*S3u)/3.0+Wkc3*F2e/3.0);

    ZS[2] = -A2/P2*CI*(UJMw+0.5*(M_ksi*S1u-3.0*M_eta*C1u)-0.5*S2u-(M_ksi*S3u-M_eta*C3u)/6.0-Wkc4*F2e/3.0);

    ZS[3] = A2P2*(-W25i*UJMw*M_eta+((1.0+0.25*Wkc1+2.25*Wkc2)*C1u-Wkc4*S1u
                                    +0.5*M_ksi*C2u+M_eta*S2u+(Wkc1-3.0*Wkc2)*C3u/12.0+Wkc4*S3u/3.0+(1.0-E2*F2e/12.0)*M_ksi
                                    -(Wkc1-3.0*Wkc2)*M_ksi*F2e/12.0)+S2*((-1.25+0.375*Wkc1-3.125*Wkc2)*C1u+0.25*Wkc4*S1u
                                                                         +0.5*M_ksi*C2u-2.0*M_eta*S2u+(7.0/12.0+11.0/48.0*Wkc1+25.0/48.0*Wkc2)*C3u
                                                                         -7.0/24.0*Wkc4*S3u+0.375*(M_ksi*C4u+M_eta*S4u)+Wkc3*C5u/16.0+Wkc4*S5u/8.0-(1.25-E2*F2e/6.0)*M_ksi
                                                                         +(0.25*F2e-F3e/6.0)*(Wkc1-3.0*Wkc2)*M_ksi));

    ZS[4] = A2P2*(W25i*UJMw*M_ksi+((1.0+1.25*Wkc1+0.25*Wkc2)*S1u-2.0*Wkc4*C1u
                                   -0.5*M_eta*C2u-E2*S3u/12.0+(1.0-E2*F2e/12.0)*M_eta
                                   +(Wkc2-3.0*Wkc1)*M_eta*F2e/12.0)+S2*((-1.75-1.125*E2)*S1u+3.25*Wkc4*C1u
                                                                        +0.5*M_ksi*S2u+2.0*M_eta*C2u+(7.0/12.0+13.0/48.0*Wkc1+23.0/48.0*Wkc2)*S3u
                                                                        +5.0/24.0*Wkc4*C3u+0.375*(M_ksi*S4u-M_eta*C4u)
                                                                        +Wkc3*S5u/16.0-Wkc4*C5u/8.0-(1.25-E2*F2e/6.0)*M_eta+(0.25*F2e-F3e/6.0)*(3.0*Wkc1-Wkc2)*M_eta));

    ZS[5] = -CI*ZS[2]+A2P2*(W23i*UJMw+(M_ksi*S1u-M_eta*C1u)+F1e*((1.0-0.25*E2)*(M_ksi*S1u-M_eta*C1u)+0.5*Wkc3*S2u
                                                                 -Wkc4*C2u+(Wkc1-3.0*Wkc2)*M_ksi*S3u/12.0+(Wkc2-3.0*Wkc1)*M_eta*C3u/12.0))
            +A2P2*S2*(-(2.0*M_ksi*S1u-M_eta*C1u)+0.75*S2u-(M_ksi*S3u-M_eta*C3u)/6.0+0.5*Wkc4*F2e
                      +F1e*((-0.5+1.25*GE2+0.125*Wkc3)*M_ksi*S1u+(2.5+1.25*GE2-0.125*(7.0*Wkc1+5.0*Wkc2))*M_ksi*C1u
                            -0.75*Wkc3*S2u+1.5*Wkc4*C2u+(1.0+5.0*GE2/12.0-7.0*Wkc1/48.0+17.0*Wkc2/48.0)*M_ksi*S3u
                            -(1.0+5.0*GE2/12.0-19.0*Wkc1/48.0+5.0*Wkc2/48.0)*M_eta*C3u
                            +0.375*Wkc3*S4u-0.75*Wkc4*C4u+(Wkc1-3.0*Wkc2)*M_ksi*S5u/16.0
                            -(3.0*Wkc1-Wkc2)*M_eta*C5u/16.0-(0.25+(1.0+0.5*E2)*F2e/3.0)*Wkc4));
    //--ZS[1]--ZS[5] END----------------------------------------------------
    double A2a = 2.0*A2/M_a;   // 2*A2/a
    double RA4 = RA * RA3;     // (a/r)^4
    double RA5 = RA * RA4;     // (a/r)^5
    double E4 = E2 * E2;       // (1-e^2)^2
    double S4 = S2 * S2;       // sin(i)^4
    double F4e = F1e/GE2;      // F4(e)
    double F5e = 0.5*(5.0+3.0*GE2-2.0*E2)*F1e; // F5(e)
    double Wkc5 = Wkc1 * Wkc1;  // ksi^4
    double Wkc6 = Wkc2 * Wkc2;  // eta^4
    double Wkc7 = Wkc1 * Wkc2;  // ksi^2 * eta^2
    double AS21 = -2.0*(ZS[0]/M_a+A2P2*GE2*W23i)*ZS[0]-A2a*(SI*CI*RA3*(1.0-C2u))*ZS[1]
            -A2a*(RA4/GE2*(M_ksi*S1u-M_eta*C1u)*(1.0-1.5*S2*(1.0-C2u))+GE2*RA5*S2*S2u)*ZS[5]
            +A2a*(RA4*(1.0-1.5*S2*(1.0-C2u))*(C1u+F4e*(Wkc2*C1u-Wkc4*S1u))
                  +RA3*S2*S2u*pow(WE2,-1.5)*(-F5e*M_eta-2.0*S1u-0.5*(M_ksi*S2u-M_eta*C2u)+0.5*F1e*M_ksi
                                             *(4.0*(M_ksi*S1u-M_eta*C1u)+Wkc3*S2u-2.0*Wkc4*C2u)))*ZS[3]
            +A2a*(RA4*(1.0-1.5*S2*(1.0-C2u))*(S1u+F4e*(Wkc1*S1u-Wkc4*C1u))
                  +RA3*S2*S2u*pow(WE2,-1.5)*(F5e*M_ksi+2.0*C1u+0.5*(M_ksi*C2u+M_eta*S2u)-0.5*F1e*M_eta
                                             *(4.0*(M_ksi*S1u-M_eta*C1u)+Wkc3*S2u-2.0*Wkc4*C2u)))*ZS[4];
    double AS22 = A2P2*A2P2*M_a*GE2*(W23i*W23i*((16.0+19.0*E2)/9.0+35.0*E4/(18.0*WE2)+2.0*GE2/9.0)
                                     +S2*(1.0+2.0*E2/3.0)+S4*(25.0*E2/24.0-5.0/6.0+35.0*E4/(16.0*WE2)
                                                              +(Wkc5-6.0*Wkc7+Wkc6)/(32.0*WE2))+S2*(5.0/6.0-1.75*S2-2.0*W25i*F2e/3.0
                                                                                                    +E2*(7.0/3.0-3.5*S2)/WE2)*Wkc3);
    AS2 = AS21 - AS22;
    //--AS2-END-------------------------------------------------------------
    ZS[0] *= 6378245.;
    AS2 *= 6378245.;
}

/// 计算函数F = 1 - eta +(m/eta**2)*W(m/eta**2-l)
/// 被函数FindEta()调用

double F (double eta, double m, double l)
{
    // Constants
    const double eps = 100.0 * eps_mach;

    // Variables
    double  w,W,a,n,g;

    w = m/(eta*eta)-l;

    if (fabs(w)<0.1)
    { // Series expansion
        W = a = 4.0/3.0; n = 0.0;
        do
        {
            n += 1.0;
            a *= w*(n+2.0)/(n+1.5);
            W += a;
        }
        while (fabs(a) >= eps);
    }
    else
    {
        if (w > 0.0)
        {
            g = 2.0*asin(sqrt(w));
            W = (2.0*g - sin(2.0*g)) / pow(sin(g), 3);
        }
        else
        {
            g = 2.0*log(sqrt(-w)+sqrt(1.0-w));  // =2.0*arsinh(sqrt(-w))
            W = (sinh(2.0*g) - 2.0*g) / pow(sinh(g), 3);
        }
    }

    return ( 1.0 - eta + (w+l)*W );
}   // End of function F

//------------------------------------------------------------------------------
//
// FindEta
//
//   根据卫星的两个位置矢量和时间间隔，求扇形和三角形面积比eta。目的是由两组位置计算轨道根数。
// Input/Output:
//
//   r_a        Position at time t_a
//   r_a        Position at time t_b
//   tau        Normalized time (sqrt(GM)*(t_a-t_b))
//   <return>   Sector-triangle ratio
//
//------------------------------------------------------------------------------
double FindEta (const CVector& r_a, const CVector& r_b, double tau)
{
    // Constants

    const int maxit = 30;
    const double delta = 100.0*eps_mach;

    // Variables

    int    i;
    double kappa, m, l, s_a, s_b, eta_min, eta1, eta2, F1, F2, d_eta;


    // Auxiliary quantities(辅助量)

    s_a = CVecMat::Norm(r_a);
    s_b = CVecMat::Norm(r_b);

    kappa = sqrt ( 2.0*(s_a*s_b+CVecMat::Dot(r_a,r_b)) );

    m = tau*tau / pow(kappa,3);
    l = (s_a+s_b) / (2.0*kappa) - 0.5;

    eta_min = sqrt(m/(l+1.0));

    // Start with Hansen's approximation

    eta2 = ( 12.0 + 10.0*sqrt(1.0+(44.0/9.0)*m /(l+5.0/6.0)) ) / 22.0;
    eta1 = eta2 + 0.1;

    // Secant method

    F1 = F(eta1, m, l);
    F2 = F(eta2, m, l);

    i = 0;

    while (fabs(F2-F1) > delta)
    {
        d_eta = -F2*(eta2-eta1)/(F2-F1);
        eta1 = eta2; F1 = F2;

        while (eta2+d_eta<=eta_min)
        {
            d_eta *= 0.5;
        }

        eta2 += d_eta;
        F2 = F(eta2,m,l);
        ++i;

        if ( i == maxit )
        {
            cerr << "WARNING: Convergence problems in FindEta" << endl;
            break;
        }
    }

    return eta2;
}

double remaining_terms( const double ival)
{
    double rval = 0., z = 1;
    const double tolerance = 1e-30;
    int i = 2;

    do
    {
        z *= ival / (double)( i * (i + 1));
        rval += z;
        i += 2;
    } while( fabs( z) > tolerance);
    return( rval);
}

/// 计算岁差的几个角度值
void compute_ecliptic_precession_angles( const double epoch_from,
                                         const double epoch_to, double *eta, double *pi, double *p)
{
    const double big_t = (epoch_from - 2000.) / 100.;
    const double big_t2 = big_t * big_t;
    const double t = (epoch_to - epoch_from) / 100.;
    const double t2 = t * t, t3 = t2 * t;

    /* See Meeus,  _Astro Algorithms_, p 128: */
    *eta = (47.0029 - .06603 * big_t + .000598 * big_t2) * t
            + (-.03302 + .000598 * big_t) * t2
            + .000060 * t3;
    *pi =  (174.876384 * 3600.) + 3289.4789 * big_t + .60622 * big_t2
            - (869.8089 + .50491 * big_t) * t + .03536 * t2;
    *p = (5029.0966 + 2.22226 * big_t - .000042 * big_t2) * t
            + (1.11113 - .000042 * big_t) * t2
            - .000006 * t3;
    /* cvt arcseconds to radians: */
    *eta *= DAS2R;
    *pi *= DAS2R;
    *p *= DAS2R;
}

/// 转换轨道根数
void convert_elements( const double epoch_from, const double epoch_to,
                       double *incl, double *asc_node, double *arg_per)
{
    double pi, eta, p, phi, x, y, z;
    const double sin_incl = sin( *incl), cos_incl = cos( *incl);
    double sin_asc_node_minus_pi, cos_asc_node_minus_pi;

    compute_ecliptic_precession_angles( epoch_from, epoch_to, &eta, &pi, &p);
    /* See Meeus,  _Astro Algorithms_, p 147: */
    phi = pi + p;
    cos_asc_node_minus_pi = cos( *asc_node - pi);
    sin_asc_node_minus_pi = sin( *asc_node - pi);
    z = cos_incl * cos( eta) + sin_incl * sin( eta) *  cos_asc_node_minus_pi;
    y = sin_incl * sin_asc_node_minus_pi;
    x = -sin( eta) * cos_incl + cos( eta) * sin_incl * cos_asc_node_minus_pi;
    *asc_node = phi + atan2( y, x);
    /* avoid ill-conditioned cases: */
    if( z > .5 || z < -.5)
        *incl = acos( z);
    else
    {
        *incl = asine( sqrt( x * x + y * y));
        if( z < 0.)                   /* added 11 Sep 2006 to fix some */
            *incl = DPI - *incl;        /* retrograde orbit errors       */
    }

    y = -sin( eta) * sin_asc_node_minus_pi;
    x = sin_incl * cos( eta) - cos_incl * sin( eta) * cos_asc_node_minus_pi;
    *arg_per += atan2( y, x);
}

/// 添加校验码
int add_tle_checksum_data( char *buff)
{
    int count = 69, rval = 0;

    if( (*buff != '1' && *buff != '2') || buff[1] != ' ')
        return( 0);    /* not a .TLE */
    while( --count)
    {
        if( *buff < ' ' || *buff > 'Z')
            return( 0);             /* wups!  those shouldn't happen in .TLEs */
        if( *buff > '0' && *buff <= '9')
            rval += *buff - '0';
        else if( *buff == '-')
            ++rval;
        ++buff;
    }
    if( *buff != 10 && buff[1] != 10 && buff[2] != 10)
        return( 0);                              /* _still_ not a .TLE */
    *buff++ = (char)( '0' + (rval % 10));
    *buff++ = 13;
    *buff++ = 10;
    *buff++ = '\0';
    return( 1);
}

/// 将数值写入到文本中
void put_sci( char *obuff, double ival)
{
    if( !ival)
        memcpy( obuff, " 00000-0", 7);
    else
    {
        int oval = 0, exponent = 0;

        if( ival > 0.)
            *obuff++ = ' ';
        else
        {
            *obuff++ = '-';
            ival = -ival;
        }

        while( ival > 1.)    /* start exponent search in floats,  to evade */
        {                 /* an integer overflow */
            ival /= 10.;
            ++exponent;
        }                 /* then do it in ints,  to evade roundoffs */
        while( oval > 99999 || oval < 10000)
        {
            oval = (int)( ival * 100000. + .5);
            if( oval > 99999)
            {
                ival /= 10;
                ++exponent;
            }
            if( oval < 10000)
            {
                ival *= 10;
                --exponent;
            }
        }
        sprintf( obuff, "%5d", oval);
        if( exponent > 0)
        {
            obuff[5] = '+';
            obuff[6] = (char)( '0' + exponent);
        }
        else
        {
            obuff[5] = '-';
            obuff[6] = (char)( '0' - exponent);
        }
    }
}

/// 将TLE数据转换成字符串
void write_elements_in_tle_format( char buff[2][73], const tle_t *tle)
{
    long year;
    double day_of_year, deriv_mean_motion;

    year = (int)( tle->epoch - J1900) / 365 + 1;
    do
    {
        double start_of_year;

        year--;
        start_of_year = J1900 + (double)year * 365. + (double)((year - 1) / 4);
        day_of_year = tle->epoch - start_of_year;
    }
    while( day_of_year < 1.);
    sprintf( buff[0],
            /*                                     xndt2o    xndd6o   bstar  eph bull */
            "1 %05d%c %-8s %02ld%12.8lf -.00000000 +00000-0 +00000-0 %c %4dZ\n",
            tle->norad_number, tle->classification, tle->intl_desig,
            year % 100L, day_of_year,
            tle->ephemeris_type, tle->bulletin_number);
    if( buff[0][20] == ' ')       /* fill in leading zeroes for day of year */
        buff[0][20] = '0';
    if( buff[0][21] == ' ')
        buff[0][21] = '0';
    deriv_mean_motion = tle->xndt2o * MINUTES_PER_DAY_SQUARED / D2PI;
    if( deriv_mean_motion >= 0)
        buff[0][33] = ' ';
    deriv_mean_motion = fabs( deriv_mean_motion * 100000000.) + .5;
    sprintf( buff[0]+35, "%08ld", (long)deriv_mean_motion);
    buff[0][43] = ' ';
    put_sci( buff[0]+44, tle->xndd6o * MINUTES_PER_DAY_CUBED / D2PI);
    put_sci( buff[0]+53, tle->bstar / AE);
    add_tle_checksum_data( buff[0]);
    // buff += strlen( buff);
    sprintf( buff[1], "2 %05d %8.4lf %8.4lf %07ld %8.4lf %8.4lf %11.8lf%5dZ\n",
            tle->norad_number,
            double(iauAnp( tle->xincl) * DR2D),
            double(iauAnp( tle->xnodeo) * DR2D),
            (long)( tle->eo * 10000000. + .5),
            double(iauAnp( tle->omegao) * DR2D),
            double(iauAnp( tle->xmo) * DR2D),
            tle->xno,
            tle->revolution_number);

    add_tle_checksum_data( buff[1]);
}

/// Harris-Priester大气密度模型
double Density_HP(const CVector &vecRSun, const CVector& r_tod )
{
    // Constants

    static const double upper_limit =   1000.0;           // Upper height limit [km]
    static const double lower_limit =    100.0;           // Lower height limit [km]
    static const double ra_lag      = 0.523599;           // Right ascension lag [rad]
    static const int    n_prm       =        3;           // Harris-Priester parameter
    // 2(6) low(high) inclination 低倾角卫星用2高倾角用6

    // Harris-Priester atmospheric density model parameters
    // Height [km], minimum density, maximum density [gm/km^3]

    static const int    N_Coef = 50;
    static const double Data_h[N_Coef]= {
        100.0, 120.0, 130.0, 140.0, 150.0, 160.0, 170.0, 180.0, 190.0, 200.0,
        210.0, 220.0, 230.0, 240.0, 250.0, 260.0, 270.0, 280.0, 290.0, 300.0,
        320.0, 340.0, 360.0, 380.0, 400.0, 420.0, 440.0, 460.0, 480.0, 500.0,
        520.0, 540.0, 560.0, 580.0, 600.0, 620.0, 640.0, 660.0, 680.0, 700.0,
        720.0, 740.0, 760.0, 780.0, 800.0, 840.0, 880.0, 920.0, 960.0,1000.0};
    static const double Data_c_min[N_Coef] = {
        4.974e+05, 2.490e+04, 8.377e+03, 3.899e+03, 2.122e+03, 1.263e+03,
        8.008e+02, 5.283e+02, 3.617e+02, 2.557e+02, 1.839e+02, 1.341e+02,
        9.949e+01, 7.488e+01, 5.709e+01, 4.403e+01, 3.430e+01, 2.697e+01,
        2.139e+01, 1.708e+01, 1.099e+01, 7.214e+00, 4.824e+00, 3.274e+00,
        2.249e+00, 1.558e+00, 1.091e+00, 7.701e-01, 5.474e-01, 3.916e-01,
        2.819e-01, 2.042e-01, 1.488e-01, 1.092e-01, 8.070e-02, 6.012e-02,
        4.519e-02, 3.430e-02, 2.632e-02, 2.043e-02, 1.607e-02, 1.281e-02,
        1.036e-02, 8.496e-03, 7.069e-03, 4.680e-03, 3.200e-03, 2.210e-03,
        1.560e-03, 1.150e-03                                            };
    static const double Data_c_max[N_Coef] = {
        4.974e+05, 2.490e+04, 8.710e+03, 4.059e+03, 2.215e+03, 1.344e+03,
        8.758e+02, 6.010e+02, 4.297e+02, 3.162e+02, 2.396e+02, 1.853e+02,
        1.455e+02, 1.157e+02, 9.308e+01, 7.555e+01, 6.182e+01, 5.095e+01,
        4.226e+01, 3.526e+01, 2.511e+01, 1.819e+01, 1.337e+01, 9.955e+00,
        7.492e+00, 5.684e+00, 4.355e+00, 3.362e+00, 2.612e+00, 2.042e+00,
        1.605e+00, 1.267e+00, 1.005e+00, 7.997e-01, 6.390e-01, 5.123e-01,
        4.121e-01, 3.325e-01, 2.691e-01, 2.185e-01, 1.779e-01, 1.452e-01,
        1.190e-01, 9.776e-02, 8.059e-02, 5.741e-02, 4.210e-02, 3.130e-02,
        2.360e-02, 1.810e-02                                            };


    // Variables

    int    i, ih;                              // Height section variables
    double height;                             // Earth flattening
    double dec_Sun, ra_Sun, c_dec;             // Sun declination, right asc.
    double c_psi2;                             // Harris-Priester modification
    double density, h_min, h_max, d_min, d_max;// Height, density parameters
    CVector u(3);                          // Apex of diurnal bulge


    // 卫星高度
    CVector vTemp;
    GisMath::XYZ2LBH(r_tod,vTemp);
    height = vTemp(2)/1000.0;              //  [km]


    // 如果超过范围则直接返回 0.0
    if ( height >= upper_limit || height <= lower_limit )
    {
        return 0.0;
    }


    // Sun right ascension, declination

    ra_Sun  = atan2( vecRSun(1), vecRSun(0) );
    dec_Sun = atan2( vecRSun(2), sqrt( pow(vecRSun(0),2)+pow(vecRSun(1),2) ) );


    // Unit vector u towards the apex of the diurnal bulge
    // in inertial geocentric coordinates

    c_dec = cos(dec_Sun);
    u(0) = c_dec * cos(ra_Sun + ra_lag);
    u(1) = c_dec * sin(ra_Sun + ra_lag);
    u(2) = sin(dec_Sun);


    // Cosine of half angle between satellite position vector and
    // apex of diurnal bulge

    c_psi2 = 0.5 + 0.5 * CVecMat::Dot(r_tod,u)/r_tod.Length();


    // Height index search and exponential density interpolation

    ih = 0;                           // section index reset
    for ( i=0; i<N_Coef-1; ++i)       // loop over N_Coef height regimes
    {
        if ( height >= Data_h[i] && height < Data_h[i+1] )
        {
            ih = i;                       // ih identifies height section
            break;
        }
    }

    h_min = ( Data_h[ih] - Data_h[ih+1] )/log( Data_c_min[ih+1]/Data_c_min[ih] );
    h_max = ( Data_h[ih] - Data_h[ih+1] )/log( Data_c_max[ih+1]/Data_c_max[ih] );

    d_min = Data_c_min[ih] * exp( (Data_h[ih]-height)/h_min );
    d_max = Data_c_max[ih] * exp( (Data_h[ih]-height)/h_max );

    /// 计算密度
    density = d_min + (d_max-d_min)*pow(c_psi2,n_prm);


    return (density * 1.0e-12);       // [kg/m^3]

}

/// 美国1976标准大气模型
double SA76(double dH)
{
    static const double hpi[8] = {0,11000,20000,32000,47000,51000,71000,84852};
    static const double tti[7] = {288.15,216.65,216.65,228.65,270.65,270.65,214.65};
    static const double eli[7] = {-0.0065,0.0,0.001,0.0028,0.0,-0.0028,-0.002};
    static const double pi[8] = {0.101325e+6, 0.226321e+5, 0.547488e+4, 0.868018e+3,
                                 0.110906e+3, 0.669387e+2, 0.395641e+1, 0.363382e+0};
    static const double feoc[25][5] = {
        {85.,   -.508515e+01, -.733375e-01, -.825405e-03,  .479193e-04},
        {90.,   -.546648e+01, -.779976e-01, -.106615e-03,  .561803e-05},
        {100.,  -.625150e+01, -.784445e-01,  .619256e-04,  .168843e-04},
        {110.,  -.701287e+01, -.721407e-01,  .568454e-03,  .241761e-04},
        {120.,  -.765326e+01, -.535188e-01,  .129374e-02, -.296655e-04},
        {130.,  -.808874e+01, -.365437e-01,  .403772e-03, -.289208e-05},
        {140.,  -.841669e+01, -.293359e-01,  .317010e-03, -.442685e-05},
        {150.,  -.868227e+01, -.243238e-01,  .184205e-03, -.144722e-05},
        {160.,  -.890904e+01, -.210738e-01,  .140788e-03, -.137461e-05},
        {170.,  -.910707e+01, -.186705e-01,  .995495e-04, -.677454e-06},
        {180.,  -.928450e+01, -.168827e-01,  .792259e-04, -.593218e-06},
        {190.,  -.944600e+01, -.154761e-01,  .614293e-04, -.381119e-06},
        {200.,  -.959500e+01, -.143619e-01,  .499958e-04, -.249568e-06},
        {220.,  -.986423e+01, -.126615e-01,  .350217e-04, -.154281e-06},
        {240.,  -.101047e+02, -.114458e-01,  .257648e-04, -.925137e-07},
        {260.,  -.103240e+02, -.105262e-01,  .202140e-04, -.774691e-07},
        {280.,  -.105271e+02, -.981064e-02,  .155659e-04, -.650883e-07},
        {300.,  -.107176e+02, -.926611e-02,  .116606e-04, -.249222e-07},
        {400.,  -.115525e+02, -.768166e-02,  .418393e-05, -.388723e-08},
        {500.,  -.122827e+02, -.696149e-02,  .301776e-05,  .447753e-08},
        {600.,  -.129442e+02, -.622361e-02,  .436102e-05,  .998741e-08},
        {700.,  -.135130e+02, -.505179e-02,  .735724e-05, -.124098e-10},
        {800.,  -.139446e+02, -.358071e-02,  .735352e-05, -.104955e-07},
        {900.,  -.142397e+02, -.242487e-02,  .420487e-05, -.833682e-08},
        {1001., -.000000e+00, -.000000e+00,  .000000e+00, -.000000e+00}};

    if( dH > 1000.e3 )
    {
        // 1000km以上的大气密度为0
        return 0;
    }
    if(dH<85.0e3)
    {
        /// 如果高度小于0,按照100m高度的大气计算
        if(dH<0)
        {
            dH=100.;
        }

        static const double g0=9.80665;
        static const double rr=287.053;
        static const double r0=6356766.;
        double hp,dhp,tt,a,b,p;
        int k,j=7;

        hp=dH/(r0+dH)*r0;
        for(k=0;k<=6;++k)
        {
            if(hpi[k]<=hp && hp<hpi[k+1]) j=k;
        }
        dhp=hp-hpi[j];
        tt=tti[j]+eli[j]*dhp;
        if(eli[j]!=0)
        {
            a=1+dhp*eli[j]/tti[j];
            b=-g0/(rr*eli[j]);
            p=pi[j]*pow(a,b);
        }
        else
        {
            a=-g0*dhp/(rr*tti[j]);
            p=pi[j]*exp(a);
        }
        return (p/(rr*tt));
    }
    else
    {
        double hk,dif,diff,difff,dsr;
        int k,n=24;

        hk=dH*.001;
        for(k=0;k<=23;++k)
        {
            if(hk>=feoc[k][0] && hk<feoc[k+1][0])
            {
                n=k;
            }
        }

        /// 如果等于 1001km
        if(fabs(hk-1001.)<1e-09)
        {
            n=24;
        }

        dif=hk-feoc[n][0];
        diff=dif*dif;
        difff=diff*dif;
        dsr=feoc[n][1]+feoc[n][2]*dif+feoc[n][3]*diff+feoc[n][4]*difff;
        return (exp(log(10.)*dsr));
    }
}

double EquationOfEquinoxes(double dpsi, double deps, double eps, double om)
{
    return(dpsi*cos(eps+deps) + DAS2R*(0.00264*sin(om) - 0.000009*sin(om+om)));
}
