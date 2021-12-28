#include <VersionMathCommon.h>
#include <Satellite/J2.h>
#include <Satellite/Kepler.h>
#include "Sofa/sofa.h"

using namespace Satellite;

CJ2::CJ2(double dA, double dE, double dI, double dRAAN, double dMA, double dAP)
    :m_dA(dA),m_dE(dE),m_dI(dI),m_dRAAN(dRAAN),m_dMA(dMA),m_dAP(dAP)
{
    m_dN=sqrt(GM_Earth/(m_dA*m_dA*m_dA));
    m_dFac = sqrt(1.0-m_dE*m_dE);
    m_dV = sqrt(GM_Earth*m_dA);
    m_vPV.Resize(6);


    dRAAN = cos(m_dI);
    dE = 4*m_dA*m_dA*m_dFac*m_dFac*m_dFac*m_dFac;
    dA = 3*DJ2*m_dN*R_Earth*R_Earth/dE;
    m_dOmega = -2*dA*dRAAN;
    m_dW = dA*(5*dRAAN*dRAAN-1);
    m_domega = m_dN + dA*(3*dRAAN*dRAAN-1)*m_dFac;
}

const CVector &CJ2::CalPV(double dT)
{
    double dM,dE;
    /// 计算平近点角位置
    if (0.0 == dT)
    {
        dM = m_dMA;
    }
    else
    {
        dM = m_dMA + m_domega*dT;
    };

    /// 计算偏近点角
    dE  = CKepler::EccAnom(dM,m_dE);

    double dCosE = cos(dE);
    double dSinE = sin(dE);


    double dR = m_dA*(1.0-m_dE*dCosE);  /// 距离
    double dV = m_dV/dR;    /// 速度


    m_vPV(0) = m_dA*(dCosE-m_dE);
    m_vPV(1) = m_dA*m_dFac*dSinE;
    m_vPV(2) = m_vPV(5) = 0.;

    m_vPV(3) = -dV*dSinE;
    m_vPV(4) = dV*m_dFac*dCosE;

    double omega = m_dAP + m_dW*dT;
    double Omega = m_dRAAN + m_dOmega * dT;
    double dTempMat[3][3];
    iauIr(dTempMat);
    iauRz(-omega,dTempMat);
    iauRx(-m_dI,dTempMat);
    iauRz(-Omega,dTempMat);
    /// 轨道面到惯性系旋转矩阵
    CMatrix  PQW(dTempMat,3,3);
    PQW.Translate(m_vPV,m_vPV);

    return(m_vPV);
}
