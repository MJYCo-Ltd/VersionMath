#include <cmath>
#include <limits>
#include <iostream>
using namespace std;
#include "Kepler.h"
#include "Vector.h"
#include "Matrix.h"
#include "VecMat.h"
#include "sofam.h"
#include "MathCommonAlgorithm.h"
#include "CommonAlgorithm.h"
#include "CoorSys.h"

using namespace Satellite;
using namespace Aerospace;
using namespace Math;

static const double eps_mach(numeric_limits<double>::epsilon());
#define centralize_angle(x) (fmod( (x) + DPI * 10., D2PI))
#define J0 (DJ00 - 2000. * 365.25)
#define SQRT_2 1.41421356

CKepler::CKepler()
{

}

CKepler::~CKepler()
{

}

/// 解开普勒方程
CVector CKepler::State(const double &dGM, const CVector &vKep, double dT)
{
    /// 变量
    double  a,e,i,Omega,omega,M,M0,n;
    double  E,cosE,sinE, fac, R,V;
    CVector  r,v;
    CMatrix  PQW;

    if(vKep.Size() < 6)
    {
        return(r);
    }

    /// 开普勒六根数
    a = vKep(0);  Omega = vKep(3);
    e = vKep(1);  omega = vKep(4);
    i = vKep(2);  M0    = vKep(5);

    /// 计算平近点角位置
    if (0.0 == dT)
    {
        M = M0;
    }
    else
    {
        n = sqrt(dGM/(a*a*a));
        M = M0 + n*dT;
    };

    /// 计算偏近点角
    E  = EccAnom(M,e);

    cosE = cos(E);
    sinE = sin(E);

    /// 轨道平面坐标系下
    fac = sqrt ( (1.0-e)*(1.0+e) );

    R = a*(1.0-e*cosE);  /// 距离
    V = sqrt(dGM*a)/R;    /// 速度

    r = CVector( a*(cosE-e), a*fac*sinE , 0.0 );
    v = CVector( -V*sinE   , +V*fac*cosE, 0.0 );

    /// 轨道面到惯性系旋转矩阵
    PQW = CVecMat::R_z(-Omega) * CVecMat::R_x(-i) * CVecMat::R_z(-omega);

    r = PQW*r;
    v = PQW*v;

    return(CVecMat::Stack(r,v));
}

/// 根据位置速度，计算经典根数
CVector CKepler::Elements(const double &dGM, const CVector &vY)
{
    CVector  r,v,h;
    double  H, u, R;
    double  eCosE, eSinE, e2, E, nu;
    double  a,e,i,Omega,omega,M;

    if(vY.Size() < 6)
    {
        return(r);
    }

    /// 读取 位置 速度
    r = vY.slice(0,2);
    v = vY.slice(3,5);

    h = CVecMat::Cross(r,v);
    H = CVecMat::Norm(h);

    Omega = atan2 ( h(0), -h(1) );                     /// 升交点赤经
    Omega = Modulo(Omega,D2PI);
    i     = atan2 ( sqrt(h(0)*h(0)+h(1)*h(1)), h(2) ); /// 轨道倾角
    u     = atan2 ( r(2)*H, -r(0)*h(1)+r(1)*h(0) );    // Arg. of latitude

    R  = CVecMat::Norm(r);                         /// 地心距离

    a = 1.0 / (2.0/R-CVecMat::Dot(v,v)/dGM);        /// 长半轴

    eCosE = 1.0-R/a;                                   // e*cos(E)
    eSinE = CVecMat::Dot(r,v)/sqrt(dGM*a);          // e*sin(E)

    e2 = eCosE*eCosE +eSinE*eSinE;
    e  = sqrt(e2);                                     /// 偏心率
    E  = atan2(eSinE,eCosE);                           /// 偏近点角

    M  = Modulo(E-eSinE,D2PI);                         /// 平近点角

    nu = atan2(sqrt(1.0-e2)*eSinE, eCosE-e2);          /// 真近点角

    omega = Modulo(u-nu,D2PI);                         /// 近地点幅角

    /// 返回开普勒六根数
    return(CVector(a,e,i,Omega,omega,M));
}

/// 根据两组位置 计算经典根数
CVector CKepler::Elements(const double &dGM, const double &dMjd1, const double &dMjd2,
                                  const CVector &vPos1, const CVector &vPos2)
{
    /// Variables

    double  tau, eta, p;
    double  n, nu, E, u;
    double  s_a, s_b, s_0, fac, sinhH;
    double  cos_dnu, sin_dnu, ecos_nu, esin_nu;
    double  a, e, i, Omega, omega, M;
    CVector  e_a(3), r_0(3), e_0(3), W(3);

    s_a = CVecMat::Norm(vPos1);  e_a = vPos1/s_a;
    s_b = CVecMat::Norm(vPos2);

    ///fac是r_b在r_a方向上的投影长度，r_0垂直于向量r_a
    fac = CVecMat::Dot(vPos2,e_a); r_0 = vPos2-fac*e_a;
    s_0 = CVecMat::Norm(r_0);  e_0 = r_0/s_0;

    /// 高斯矢量W
    W     = CVecMat::Cross(e_a,e_0);

    /// Long. ascend. node 升交点赤经
    Omega = atan2 ( W(0), -W(1) );
    Omega = Modulo(Omega,D2PI);

    /// Inclination   倾角
    i     = atan2 ( sqrt(W(0)*W(0)+W(1)*W(1)), W(2) );

    if (i==0.0)
    {
        /// u是纬度辐角（u=w+i）,当i=0时，轨道平面与赤道平面重合
        u = atan2 ( vPos1(1), vPos1(0) );
    }
    else
    {
        u = atan2 ( +e_a(2) , -e_a(0)*W(1)+e_a(1)*W(0) );
    }

    /// Semilatus rectum

    tau = sqrt(dGM) * 86400.0*fabs(dMjd2-dMjd1);
    eta = FindEta ( vPos1, vPos2, tau );
    p   = pow ( s_a*s_0*eta/tau, 2 );                 /// 半通径

    /// Eccentricity, true anomaly and argument of perihelion

    cos_dnu = fac / s_b;
    sin_dnu = s_0 / s_b;

    ecos_nu = p / s_a - 1.0;
    esin_nu = ( ecos_nu * cos_dnu - (p/s_b-1.0) ) / sin_dnu;

    e  = sqrt ( ecos_nu*ecos_nu + esin_nu*esin_nu ); /// 偏心率
    nu = atan2(esin_nu,ecos_nu);                     /// 真近点角

    omega = Modulo(u-nu,D2PI);                       /// 近地点辐角

    /// 半长轴和平动角速度

    a = p/(1.0-e*e);
    n = sqrt ( dGM / fabs(a*a*a) );

    /// 平近点角

    if (e<1.0)
    {
        E = atan2 ( sqrt((1.0-e)*(1.0+e)) * esin_nu,  ecos_nu + e*e );
        M = Modulo ( E - e*sin(E), D2PI );
    }
    else
    {
        sinhH = sqrt((e-1.0)*(e+1.0)) * esin_nu / ( e + e * ecos_nu );
        M = e * sinhH - log ( sinhH + sqrt(1.0+sinhH*sinhH) );
    }

    /// 返回开普勒六根数
    return(CVector(a,e,i,Omega,omega,M));
}

/// 转换成经典根数
CVector CKepler::ClassicalElements(const double& dGM, const CVector& vY)
{
    double h0, n0, tval;
    double ecc2;
    CVector h(3),e(3),vecR,vecV;
    double ecc;
    int i;
    double  da(0),de,di,dOmega,domega,dM;

    /// 如果没有位置、速度则退出
    if(vY.Size() < 6)
    {
        return(vecR);
    }

    /// 获取位置、速度
    vecR = vY.slice(0,2);
    vecV = vY.slice(3,5);

    const double dist = vecR.Length();
    const double v2 = CVecMat::Dot(vecV,vecV);
    double inv_major_axis = 2. / dist - v2 / dGM;
    double dNInvMajo;

    h = CVecMat::Cross(vecR,vecV);

    n0 = h(0) * h(0) + h(1) * h(1);
    h0 = n0 + h(2) * h(2);

    n0 = sqrt( n0);
    h0 = sqrt( h0);

    /// 升交点赤经
    dOmega = Modulo(atan2( h(0), -h(1)),D2PI);

    /// 轨道倾角
    {
        di = asine( n0 / h0);

        if( h(2) < 0.)                   /* 逆行轨道 */
        {
            di = DPI - di;
        }
    }

    e = CVecMat::Cross(vecV,h);

    for( i = 0; i < 3; ++i)
    {
        e(i) = e(i) / dGM - vY(i) / dist;
    }

    tval = CVecMat::Dot( e, h) / h0;     /* "flatten" e vector into the rv */

    for( i = 0; i < 3; ++i)                  /* plane to avoid roundoff; see   */
    {
        e(i) -= h(i) * tval;                 /* above comments                 */
    }

    ecc2 = CVecMat::Dot( e, e);

    double dSqrtE;
    dSqrtE = sqrt( fabs( 1. - ecc2));

    /// 轨道离心率
    ecc = de = sqrt( ecc2);

    if( fabs(ecc)<1e-3)            /* for purely circular orbits,  e is */
    {                              /* arbitrary in the orbit plane; choose */
        for( i = 0; i < 3; ++i)    /* r normalized                         */
        {
            e(i) = vY(i) / dist;
        }
    }
    else                           /* ..(0)nd if it's not circular,  */
    {
        for( i = 0; i < 3; ++i)     /* normalize e:                  */
        {
            e(i) /= ecc;
        }
    }

    if( ecc < .9)
    {
        dNInvMajo = (1. - ecc) / inv_major_axis;
    }
    else        /* at eccentricities near one,  the above suffers  */
    {           /* a loss of precision problem,  and we switch to: */
        const double gm_over_h0 = dGM / h0;
        const double perihelion_speed = gm_over_h0 +
                sqrt( gm_over_h0 * gm_over_h0 - inv_major_axis * dGM);

        dNInvMajo = h0 / perihelion_speed;
        inv_major_axis = (1. - ecc) / dNInvMajo;
    }

    /// 轨道长半轴
    if( inv_major_axis)
    {
        da = 1. / inv_major_axis;
    }

    CVector vTemp = CVecMat::Cross(h,e);

    /// 近地点幅角
    {
        const double cos_arg_per = (h(0) * e(1) - h(1) * e(0)) / n0;

        if( cos_arg_per < .7 && cos_arg_per > -.7)
        {
            domega = acos( cos_arg_per);
        }
        else
        {
            const double sin_arg_per =
                    (e(0) * h(0) * h(2) + e(1) * h(1) * h(2) - e(2) * n0 * n0)
                    / (n0 * h0);

            domega = fabs(asin( sin_arg_per));
            if( cos_arg_per < 0.)
            {
                domega = DPI - domega;
            }
        }

        if( e(2) < 0.)
        {
            domega = D2PI - domega;
        }

        domega = Modulo(domega,D2PI);
    }

    /// 平近点角
    if( inv_major_axis && dSqrtE)
    {
        const double r_cos_true_anom = CVecMat::Dot(vecR, e);
        const double r_sin_true_anom = CVecMat::Dot(vecR, vTemp) / h0;
        const double cos_E = r_cos_true_anom * inv_major_axis + ecc;
        const double sin_E = r_sin_true_anom * inv_major_axis / dSqrtE;


        if( inv_major_axis > 0.)          /* 椭圆 */
        {
            const double ecc_anom = atan2( sin_E, cos_E);

            dM = ecc_anom * (1 - ecc)
                    - ecc * ecc_anom * remaining_terms( -ecc_anom * ecc_anom);

        }
        else                             /* 双曲线 */
        {
            /// G++ 用
            /// const double ecc_anom = atanh( sin_E / cos_E);
            /// VS 用
            double dTemp = sin_E/cos_E;
            const double ecc_anom = 0.5*log((1+dTemp)/(1-dTemp));

            dM = ecc_anom * (1 - ecc)
                    - ecc * ecc_anom * remaining_terms( ecc_anom * ecc_anom);
        }
    }
    else              /* 抛物线 */
    {
        dM = (3. / SQRT_2) / (dNInvMajo * sqrt( dNInvMajo / dGM));
    }

    dM = Modulo(dM,D2PI);

    return(CVector(da,de,di,dOmega,domega,dM));
}

/// 将经典六根数转化成 TLE数据
///
bool CKepler::Classical2TLE(const CVector &vKep, const double &dJDEpoch,
                                int nNorad, char str[2][73])
{
    /// 判断数据是否满足
    if(vKep.Size() < 6 || vKep(1)>.99)
    {
        return(false);
    }

    tle_t tle;
    /// 计算卫星周期
    const double t0 =
        vKep(0) * sqrt( vKep(0) / 398600.4415e+9);
    double incl = vKep(2);
    double asc_node = vKep(3);
    double arg_per = vKep(4);
    double mean_anomaly = vKep(5);

    tle.epoch = dJDEpoch;
          /* The elements are in J2000,  but TLEs are given  */
          /* in epoch of date: */
    convert_elements( 2000., (dJDEpoch - J0) / DJY,
                     &incl, &asc_node, &arg_per);     /* conv_ele.cpp */
    tle.xincl = centralize_angle( incl);
    tle.xnodeo = centralize_angle( asc_node);
    tle.omegao = centralize_angle( arg_per);
    tle.xmo = centralize_angle( mean_anomaly);
    tle.xno = D2PI*t0/DAYSEC;
    tle.eo = vKep(1);
          /* Address these three values later: */
    tle.xndt2o = tle.xndd6o = tle.bstar = 0.;

    /// 自己添加的信息
    tle.norad_number = nNorad;
    tle.intl_desig[0] = '\0';
    tle.classification = 'U';
    tle.ephemeris_type = '0';

    write_elements_in_tle_format(str,&tle);

    return(true);
}

/// 根据根数计算两者的运动关系
CVector CKepler::RIC(const CVector &vKepChaser,
                             const CVector &vKepTarget, const double &dGM)
{
    CVector tmp(6);
    /// 判断数据是否有效
    if(vKepChaser.Size() < 6 || vKepTarget.Size() < 6)
    {
        return(tmp);
    }

    CVector vChaser = State(dGM,vKepChaser);
    CVector vTarget = State(dGM,vKepTarget);

    /// 如果计算失败则返回
    if(!vChaser || !vTarget)
    {
        return(tmp);
    }

    CVector rc=vChaser.slice(0,2),
                vc=vChaser.slice(3,5),
                rt=vTarget.slice(0,2),
                vt=vTarget.slice(3,5),
            RelPos,RelVel;

    CMatrix Coti = CCoorSys::ECI2VVLH(vKepTarget);
    RelPos = Coti*(rc - rt);
    //Vector w(0,-RefSat.n(),0);
    // 这里的轨道角速度应该用位置速度来计算，与使用平均轨道
    // 相比，相对速度有差别，差别大小与目标轨道偏心率有关
    CVector w = CVecMat::Cross(rt,vt)/rt.Length()/rt.Length();
    w = Coti*w;
    RelVel = Coti*(vc - vt) - CVecMat::Cross(w,RelPos);

    /// 赋值操作
    tmp(0) = RelPos(0);
    tmp(1) = RelPos(1);
    tmp(2) = RelPos(2);
    tmp(3) = RelVel(0);
    tmp(4) = RelVel(1);
    tmp(5) = RelVel(2);

    return(tmp);
}

CVector CKepler::CIR(const CVector &vKepTarget, const CVector& vPV, const double& dGM)
{
    CVector tmp;
    /// 判断数据是否有效
    if(vKepTarget.Size() < 6 || vPV.Size() <6)
    {
        return(tmp);
    }

    /// 获取ECI 到 轨道面的旋转矩阵
    CMatrix Coti = CCoorSys::ECI2VVLH(vKepTarget);
    CVector rv;
    rv = State(dGM,vKepTarget);
    if(!rv)
    {
        return(rv);
    }
    //Vector w(0,-RefSat.n(),0);
    // 这里的轨道角速度应该用位置速度来计算，与使用平均轨道
    // 相比，相对速度有差别，差别大小与目标轨道偏心率有关

    CVector rt(rv.slice(0,2));
    CVector vt(rv.slice(3,5));

    CVector w = CVecMat::Cross(rt,vt)/rt.Length()/rt.Length();
    w = Coti*w;
    CMatrix Itoc(CVecMat::Transp(Coti));
    CVector cPos = rt + Itoc*vPV.slice(0,2); //目标在惯性坐标系下的位置矢量
    CVector cVel = vt + Itoc*(vPV.slice(3,5) + CVecMat::Cross(w,vPV.slice(0,2))); //目标在惯性坐标系下的速度矢量
    CVector endRv(6);

    endRv(0) = cPos(0);
    endRv(1) = cPos(1);
    endRv(2) = cPos(2);
    endRv(3) = cVel(0);
    endRv(4) = cVel(1);
    endRv(5) = cVel(2);

    return ClassicalElements(dGM,endRv);   //由目标惯性系位置速度计算其轨道根数
}

/// 将瞬时根数转换成平根数
CVector CKepler::Instant2Mean(const CVector &rIKep)
{
    double ST[6],AS2;
    CVector ME(6);
    if(rIKep.Size() < 6)
    {
        return(ME);
    }

    int count=0;
    double aeps,ksieps,etaeps,ieps,oeps,lmdeps;
    double a2,ksi2,eta2,i2,o2,lamda2;
    double ksi = rIKep(1)*cos(rIKep(4));
    double eta = rIKep(1)*sin(rIKep(4));
    double lamda = Modulo(rIKep(4) + rIKep(5),D2PI);
    double mksi,meta,mlamda;

    while(true)
    {
        Short2(rIKep,ST,AS2);
        a2 = ME(0) + ST[0] + AS2;
        i2 = ME(2) + ST[1];
        o2 = ME(3) + ST[2];
        ksi2 = ME(1)*cos(ME(4)) + ST[3];
        eta2 = ME(1)*sin(ME(4)) + ST[4];
        lamda2 = ME(4) + ME(5) + ST[5];
        aeps = a2 - rIKep(0);
        ieps = i2 - rIKep(2);
        oeps = Modulo(o2 - rIKep(3),D2PI);
        ksieps = ksi2 - ksi;
        etaeps = eta2 - eta;
        lmdeps = Modulo(lamda2 - lamda,D2PI);

        if( fabs(aeps)<1e-6 && fabs(ksieps)<1e-10 && fabs(etaeps)<1e-10
         && fabs(ieps)<1e-8 && fabs(oeps)<1e-8 && fabs(lmdeps)<1e-8 )
        {
            break;
        }

        /// 最多循环 50次
        if( ++count>50 )
        {
            break;
        }

        ME(0) = rIKep(0) - ST[0] - AS2;
        ME(2) = rIKep(2) - ST[1];
        ME(3) = rIKep(3) - ST[2];
        mksi = rIKep(1)*cos(rIKep(4)) - ST[3];
        meta = rIKep(1)*sin(rIKep(4)) - ST[4];
        mlamda = rIKep(4) + rIKep(5) - ST[5];
        ME(1) = sqrt(mksi*mksi+meta*meta);
        ME(4) = Modulo(atan2(meta,mksi)+D2PI,D2PI);
        ME(5) = Modulo(mlamda-ME(4)+2*D2PI,D2PI);
    }
    return (ME);
}

/// 将平根数转换成瞬时根数
CVector CKepler::Mean2Instant(const CVector &rMKep)
{
    CVector InstKepler(6);
    if(rMKep.Size() < 6)
    {
        return(InstKepler);
    }

    double ShortTerm[6],as2;//ShortTerm[6]:短周期项,as2:半长轴的二阶短周期项
    Short2(rMKep,ShortTerm,as2);

    InstKepler(0) = rMKep(0) + ShortTerm[0] + as2;
    InstKepler(2) = rMKep(2) + ShortTerm[1];
    InstKepler(3) = rMKep(3) + ShortTerm[2];
    double ksi = rMKep(1)*cos(rMKep(4)) + ShortTerm[3];
    double eta = rMKep(1)*sin(rMKep(4)) + ShortTerm[4];
    double lamda = rMKep(4) + rMKep(5) + ShortTerm[5];
    InstKepler(1) = sqrt(ksi*ksi + eta*eta);
    InstKepler(4) = Modulo(atan2(eta,ksi) + D2PI,D2PI);
    InstKepler(5) = Modulo(lamda - InstKepler(4) + 2*D2PI,D2PI);
    return (InstKepler);
}

/// 计算真近点角
double CKepler::TrueAnom(double de, double dM)
{
    double E = EccAnom(dM,de);
    if(de<=1.0)
    {
        return (Modulo(atan2(sqrt(1-de*de)*sin(E),cos(E)-de)+D2PI,D2PI));
    }
    else
    {
        return (Modulo(atan2(-sqrt(de*de-1)*sinh(E),cosh(E)-de)+D2PI,D2PI));
    }
}

/// 真近点角转平近点角
double CKepler::MeanAnom(double da, double de, double dT)
{
    double r = da*(1-de*de)/(1+de*cos(dT));
    // 对于椭圆和双曲线，下式在形式上基本相同，只是在根号下差一个符号，只要对1-e^2取绝对值就可以了
    double E;
    if(de<1.0)
    {
        E = Modulo(atan2(r*sin(dT)/da/sqrt(1-de*de),r/da*cos(dT)+de)+D2PI,D2PI);
        return(E - de*sin(E));
    }
    else if(de>1.0)
    {
        double x = (-r*sin(dT)/da/sqrt(de*de-1)) / (r/da*cos(dT)+de);
        E = 0.5*log((1+x)/(1-x));
        return(de*sinh(E) - E);
    }

    return(dT);
}

/// 计算偏近点角
double CKepler::EccAnom(double dM, double de)
{
    // Constants
    const int maxit = 15;
    const double eps = 100.0*eps_mach;

    // Variables

    int    i=0;
    double E, f;

    // Starting value

    dM = Modulo(dM, D2PI);
    if (de<0.8)
    {
        E=dM;
    }
    else
    {
        E=DPI;
    }

    /// 牛顿迭代
    do
    {
        f = E - de*sin(E) - dM;
        E = E - f / ( 1.0 - de*cos(E) );
        ++i;
        if (i==maxit)
        {
            cerr << " convergence problems in EccAnom" << endl;
            break;
        }
    }
    while (fabs(f) > eps);

    return(E);
}

/// 根据轨道的长半轴计算卫星的运行周期
double CKepler::T(double da, double dGM)
{
    return(D2PI/sqrt(dGM/pow(da,3.)));
}

/// 计算轨道的远地点
double CKepler::Apogee(double da, double de)
{
    return(da*(1+de));
}

/// 计算轨道的近地点
double CKepler::Perigee(double da, double de)
{
    return(da*(1-de));
}
