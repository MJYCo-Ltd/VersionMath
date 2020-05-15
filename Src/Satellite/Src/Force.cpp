#include <cmath>
using namespace std;
#include "VecMat.h"
#include "Force.h"
#include "CoorSys.h"
#include "sofam.h"
#include "MathCommonAlgorithm.h"
#include "CommonAlgorithm.h"
#include "GisMath.h"

using namespace Physical;
using namespace Math;
using namespace Aerospace;
CForce::CForce()
{

}

CForce::~CForce()
{

}

/// 计算地球非球形引力
CVector CForce::AccelHarmonic (const CVector& vecR, const CMatrix& matT,
                                       double dGM, double dR, const CMatrix& matCS,
                                       int nN_Max, int nM_Max )
{

  // Local variables

  int     n,m;                           // Loop counters
  double  r_sqr, rho, Fac;               // Auxiliary quantities
  double  x0,y0,z0;                      // Normalized coordinates
  double  ax,ay,az;                      // Acceleration vector
  double  C,S;                           // Gravitational coefficients
  CVector  r_bf(3);                  // Body-fixed position
  CVector  a_bf(3);                  // Body-fixed acceleration

  n = matCS.Row();
  m = matCS.Col();

  nN_Max = nN_Max > n ? n : nN_Max;
  nM_Max = nM_Max > m ? m : nM_Max;

  CMatrix  V(nN_Max+2,nN_Max+2);       // Harmonic functions
  CMatrix  W(nN_Max+2,nN_Max+2);       // work array (0..n_max+1,0..n_max+1)


  /// 地固系下位置

  r_bf = matT * vecR;

  // Auxiliary quantities

  /// 到地心的距离
  r_sqr =  CVecMat::Dot(r_bf,r_bf);
  rho   =  dR*dR / r_sqr;

  x0 = dR * r_bf(0) / r_sqr;          // Normalized
  y0 = dR * r_bf(1) / r_sqr;          // coordinates
  z0 = dR * r_bf(2) / r_sqr;

  //
  // Evaluate harmonic functions
  //   V_nm = (R_ref/r)^(n+1) * P_nm(sin(phi)) * cos(m*lambda)
  // and
  //   W_nm = (R_ref/r)^(n+1) * P_nm(sin(phi)) * sin(m*lambda)
  // up to degree and order n_max+1
  //

  // Calculate zonal terms V(n,0); set W(n,0)=0.0

  V(0,0) = dR / sqrt(r_sqr);
  W(0,0) = 0.0;

  V(1,0) = z0 * V(0,0);
  W(1,0) = 0.0;

  for (n=2; n<=nN_Max+1; ++n)
  {
    V(n,0) = ( (2*n-1) * z0 * V(n-1,0) - (n-1) * rho * V(n-2,0) ) / n;
    W(n,0) = 0.0;
  };

  // Calculate tesseral and sectorial terms

  for (m=1; m<=nM_Max+1; ++m)
  {

    // Calculate V(m,m) .. V(n_max+1,m)

    V(m,m) = (2*m-1) * ( x0*V(m-1,m-1) - y0*W(m-1,m-1) );
    W(m,m) = (2*m-1) * ( x0*W(m-1,m-1) + y0*V(m-1,m-1) );

    if (m<=nN_Max)
    {
      V(m+1,m) = (2*m+1) * z0 * V(m,m);
      W(m+1,m) = (2*m+1) * z0 * W(m,m);
    };

    for (n=m+2; n<=nN_Max+1; ++n)
    {
      V(n,m) = ( (2*n-1)*z0*V(n-1,m) - (n+m-1)*rho*V(n-2,m) ) / (n-m);
      W(n,m) = ( (2*n-1)*z0*W(n-1,m) - (n+m-1)*rho*W(n-2,m) ) / (n-m);
    };

  };

  //
  // Calculate accelerations ax,ay,az
  //

  ax = ay = az = 0.0;

  for (m=0; m<=nM_Max; ++m)
  {
      for (n=m; n<=nN_Max ; ++n)
      {
          if (m==0)
          {
              C = matCS(n,0);   // = C_n,0
              ax -=       C * V(n+1,1);
              ay -=       C * W(n+1,1);
              az -= (n+1)*C * V(n+1,0);
          }
          else
          {
              C = matCS(n,m);   // = C_n,m
              S = matCS(m-1,n); // = S_n,m
              Fac = 0.5 * (n-m+1) * (n-m+2);
              ax +=   + 0.5 * ( - C * V(n+1,m+1) - S * W(n+1,m+1) )
                      + Fac * ( + C * V(n+1,m-1) + S * W(n+1,m-1) );
              ay +=   + 0.5 * ( - C * W(n+1,m+1) + S * V(n+1,m+1) )
                      + Fac * ( - C * W(n+1,m-1) + S * V(n+1,m-1) );
              az += (n-m+1) * ( - C * V(n+1,m)   - S * W(n+1,m)   );
          };
      }
  }

  // Body-fixed acceleration

  a_bf = (dGM/(dR*dR)) * CVector(ax,ay,az);

  // Inertial acceleration

  return  CVecMat::Transp(matT)*a_bf;

}

/// 带谐项引力
CVector CForce::ZonalHarmonic(const CVector &vecR, const CMatrix &matT,
                                      double dGM, double dR, const CMatrix &matCS,
                                      int nN_Max)
{
    static CVector k(0,0,1);
    static double PL[71];

    CVector SatPos = matT * vecR;
    SatPos /= dR;
    double r = SatPos.Length();
    /// 计算 sin
    double zr = SatPos.GetZ()/r;

    if(nN_Max<2)
    {
        //即不考虑地球非球形引力摄动
        return (CVector(0,0,0));
    }

    /// 计算勒让德函数
    Legendre_sphPl(nN_Max,zr,PL);

    double WP,RL,dPL,pVpr=0,pVpsf=0;
    double dTemp = 1.0 - zr*zr;
    for(int L=2;L<=nN_Max;++L)
    {
        /// 计算勒让德函数的导数
        WP = sqrt((2.0*L+1.0)/(2.0*L-1.0));
        dPL = L/dTemp*(WP*PL[L-1]-zr*PL[L]);

        /// 计算系数
        RL = pow(1/r,L+3);
        pVpr += matCS(L,0)*RL*((L+1.0)*PL[L]+zr*dPL);
        pVpsf += matCS(L,0)*RL*dPL*r;

    }

    /// 将无纲量数据变成 m/s^2的数据
    double miurrr = dGM/dR/dR;
    pVpr  *= miurrr;
    pVpsf *= miurrr;

    return (CVecMat::Transp(matT)*(k*pVpsf - pVpr*SatPos));
}

/// 田谐项引力
CVector CForce::TesseralHarmonic(const CVector &vecR, const CMatrix &matT,
                                         double dGM, double dR, const CMatrix &matCS,
                                         int nN_Max)
{

    static CVector k(0,0,1);
    static double PLM[71][71]={{0},{0}},SX[71]={0},CX[71]={0};

    CVector SatPos = matT * vecR;

    double dL,dB,dH;
    GisMath::XYZ2LBH(SatPos(0),SatPos(1),SatPos(2),dL,dB,dH);

    SatPos /= dR;
    double RL,dPlm,WCS,WSC,pVpr=0,pVpsf=0,pVpl=0;
    double r = SatPos.Length();
    double zr = SatPos.GetZ()/r;

    /// 计算勒让德函数
    Legendre_sphPlm(nN_Max,zr,PLM);

    double dCosL = cos(dL);
    double dSinL = sin(dL);

    double gxy = sqrt(SatPos(0)*SatPos(0)+SatPos(1)*SatPos(1));

    CVector g(-dSinL,dCosL,0.0);

    /// 计算Cos mX Sin mX
    SmxCmx(nN_Max, dSinL, dCosL, SX, CX);

    double dTemp = 1.0 / sqrt(1.0 - zr*zr);

    double dRL;
    for(int L=2;L<=nN_Max;++L)
    {
        RL = pow(1/r,L+3);
        dRL = pow(1/r,L+1);
        for(int M=1;M<=L;++M)
        {
            /// 计算勒让德函数的导数
            dPlm = dTemp *(sqrt((double)(L+M+1)*(L-M))*PLM[L-1][M]- dTemp*M*zr*PLM[L][M]);

            WCS = matCS(L,M)*CX[M] + matCS(M-1,L)*SX[M];
            WSC = matCS(L,M)*SX[M] - matCS(M-1,L)*CX[M];

            pVpr -= RL*((L+1.0)*PLM[L][M]+zr*dPlm)*WCS;
            pVpsf += RL*dPlm*r*WCS;
            pVpl -= dRL*M/gxy*PLM[L][M]*WSC;
        }
    }

    double miurrr = dGM/dR/dR;
    pVpr *= miurrr;
    pVpsf *= miurrr;
    pVpl *= miurrr;

    return (CVecMat::Transp(matT)*(pVpr*SatPos + pVpsf*k + pVpl*g));
}

/// 太阳光压
CVector CForce::AccelSolrad(const CVector &vecR, const CVector &vecRSun,
                                    double dArea, double dMass, double dCR,
                                    double dP0, double dAU)
{

    CVector d = vecR - vecRSun;

    /// 太阳相对于卫星的位置

    /**
        *             《卫星轨道模型方法和应用》
        *   (德)门斯布吕克(Oliver Montenbruck),(德)吉尔(Eberhard Gill)著
        *            王家松        祝开建        胡小工      译
        *    国防工业出版社,2012.4第1版第1次印刷
        *    75页
        */
    return  dCR*(dArea/dMass)*dP0*(dAU*dAU) * d / pow(CVecMat::Norm(d),3);
}

/// 大气阻力
CVector CForce::AccelDrag(const CVector &vecRSun, const CVector &vecR,
                                  const CVector &vecV, const CMatrix &matT, double dMJD,
                                  double dArea, double dMass, double dCD, AtmosModel typeAtmos)
{
    /// 地球自转角速度
    const CVector omega ( 0.0, 0.0, 7.29212e-5 );


    // Variables

    double v_abs, dens;
    CVector r_tod(3), v_tod(3);
    CVector v_rel(3), a_tod(3);
    CMatrix T_trp(3,3);


    /// 计算 ECF 到 ECI 矩阵

    T_trp = CVecMat::Transp(matT);


    /// 将卫星位置速度转换到地固系下

    r_tod = matT * vecR;
    v_tod = matT * vecV;


    /// 卫星相对于大气运动的速度

    v_rel = v_tod - CVecMat::Cross(omega,r_tod);
    v_abs = CVecMat::Norm(v_rel);


    /// 大气密度模型 Harris-Priester

    // dens = Density_HP (vecRSun,r_tod);
    dens = CAtmosphere::GetDensity(dMJD,r_tod,vecRSun,typeAtmos);

    /// 计算受力情况

    a_tod = -0.5*dCD*(dArea/dMass)*dens*v_abs*v_rel;

    return T_trp * a_tod;
}

/// 第三体引力
CVector CForce::AccelPointMass(const CVector &vecR, const CVector &vecPoint, double dGM)
{
    /// 计算相对于卫星的位置
    CVector d = vecR - vecPoint;

    // Acceleration

    return  (-dGM) * ( d/pow(CVecMat::Norm(d),3) + vecPoint/pow(CVecMat::Norm(vecPoint),3) );
}

/// 后牛顿效应
CVector CForce::AccelPostNewton(const CVector &vecR, const CVector &vecV, double dGM)
{
    double c2 = DLIGHT*DLIGHT;
    double r1 = vecR.Length();
    double r3 = r1*r1*r1;
    double c2r3 = c2*r3;
    double v2 = vecV.Length();
    double r4v2 = 4.0 * dGM / r1 - v2 * v2;
    double rv4 = 4.0 * CVecMat::Dot(vecR,vecV);
    return (vecR*r4v2 + vecV*rv4) * dGM / c2r3;
}
