#include <Math/VecMat.h>
#include <Math/Intpol.h>
#include <Satellite/TimeSys.h>
#include <Satellite/OpticalAlg.h>
#include <VersionMathCommon.h>
#include "sofa.h"

using namespace Math;
using namespace Satellite;
using namespace Numerical;

COpticalAlg::COpticalAlg()
{

}

COpticalAlg::~COpticalAlg()
{

}

/// 计算天体阴影
double COpticalAlg::Shadow(const CVector &vecPos, const CVector &vecSun,
                               const CVector &vecCbPos, const double &dRcb)
{
    CVector vecRS = vecPos - vecSun;
    CVector vecRC = vecPos - vecCbPos;

    double rm = vecRC.Length(), drm = vecRS.Length();
    double thetaES = acos(CVecMat::Dot(vecRC,vecRS)/rm/drm );  // 天体-卫星-太阳夹角
    double as = asin(R_Sun/drm);  // 卫星上看太阳的视半角
    double ae = asin(dRcb/rm);    // 卫星上看天体的视半角

    double ae2 = ae*ae;
    double as2 = as*as;
    double as2Inv = 1.0/as2;

    if( thetaES <= ae-as)
    {
        return 0.0; // 本影区，日全食
    }
    else if( thetaES < ae+as )
    {
        if(as-thetaES>=ae)
        {
            return 1.0-ae2*as2Inv;  // 天体比太阳视半径小，日环食
        }
        else
        {
            double a = acos( (ae2+thetaES*thetaES-as2)/(2.0*ae*thetaES) );
            double b = acos( (as2+thetaES*thetaES-ae2)/(2.0*as*thetaES) );
            // 太阳可视面积与太阳圆面积的比值，日偏食
            return 1.0 - ( as2*(b-sin(b)*cos(b)) + ae2*(a-sin(a)*cos(a)))*as2Inv/DPI;
        }
    }
    else
    {
        return 1.0; // 阳照区
    }
}

/// 对流层对光(电磁波)的折射
double COpticalAlg::TroposphericRef(double dPa, double dfH, double dT, double dElev)
{
    double dTC = dT+T0;                                  // Temperature [C]
    double eh = 6.10*dfH*exp(17.15*dTC/(234.7+dTC));     // Partial water pressure
    double Ns = 77.64*dPa/dT + 3.734e5*eh/(dT*dT);       // Refractivity
    return(Ns*1.0e-6/tan(dElev));                        // Tropospheric refraction
}

/// 下行链路的光行迭代(卫星到地面站)
double COpticalAlg::DownlegLightTime(const vector<CVector> &vSatPos, const CVector &vecStaion,
                                         int nPos, double dMJd, double dStep, CVector &vecSat)
{
    /// 判断是否有效
    if(int(vSatPos.size()) <= nPos)
    {
        return(-1);
    }

    CVector vecTem;
    double tau_down = 0.0,dStarJD,dRho;
    double d[3],dX[3],dY[3],dZ[3];
    int nSize = vSatPos.size() - nPos;
    int nMax = nSize > 3 ? 3 : nSize;

    dStarJD = dMJd + nPos*dStep/DAYSEC;
    for(int j=0; j<nMax; ++j)
    {
        vecTem = vSatPos[j];
        d[j] = dStarJD + j*dStep/DAYSEC;
        dX[j] = vecTem(0);
        dY[j] = vecTem(1);
        dZ[j] = vecTem(2);
    }

    vecSat.Resize(3);
    /// 经过两次迭代才能准确 (经验值)
    for(int i=0; i<2; ++i)
    {
        vecSat(0) = Intpol::ItNeville(nMax,d,dX,dStarJD-tau_down/DAYSEC);
        vecSat(1) = Intpol::ItNeville(nMax,d,dY,dStarJD-tau_down/DAYSEC);
        vecSat(2) = Intpol::ItNeville(nMax,d,dZ,dStarJD-tau_down/DAYSEC);


        /// 计算两者的距离
        dRho = CVecMat::Norm(vecSat - vecStaion);

        /// 光传播的时间
        tau_down = dRho/DLIGHT;
    }

    return(tau_down);
}

/// 上行链路的光行迭代(地面站到卫星)
double COpticalAlg::UplegLightTime(const vector<CVector> &vSatPos, const CVector &vecStaion,
                                       int nPos, double dMJd, double dStep, CVector &vecSat)
{
    /// 判断是否有效
    if(int(vSatPos.size()) <= nPos)
    {
        return(-1);
    }

    vecSat = vSatPos[nPos];
    double tau_up = 0.0,dStarJD,dRho,dEartAngle;

    dStarJD = dMJd + nPos*dStep/DAYSEC;

    CVector vecTem = vecStaion;

    /// 经过两次迭代才能准确 (经验值)
    for(int i=0; i<2; ++i)
    {
        CTimeSys tmpTimeSys(dStarJD-tau_up/DAYSEC);
        dEartAngle = iauGmst06(DJM0,tmpTimeSys.GetUT1(),DJM0,tmpTimeSys.GetTT());
        vecTem = CVecMat::Transp(CVecMat::R_z(dEartAngle))*vecStaion;
        /// 计算两者的距离
        dRho = CVecMat::Norm(vecTem - vecStaion);

        /// 光传播的时间
        tau_up = dRho/DLIGHT;
    }

    return(tau_up);
}
