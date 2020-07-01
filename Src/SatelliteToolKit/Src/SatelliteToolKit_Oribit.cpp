#include "VecMat.h"
#include "sofa.h"
#include "Date.h"
#include "CoorSys.h"
#include "Kepler.h"
#include "SGP4.h"
#include "SatelliteToolKit.h"
#include "../Inc/SatelliteToolKitCommon.h"

using namespace Math;
using namespace Satellite;
using namespace Aerospace;

/// 根据开普勒六根数计算轨道
bool TwoBody(const BJTime &stStartTime, const BJTime &stEndTime,
             unsigned int nStep, const Kepler &stKepler, SatellitePos &stSatPos)
{

    double dMJDStart,dMJDEnd;

    if(!JudgeTimeValid(stStartTime,stEndTime,dMJDStart,dMJDEnd))
    {
        return(false);
    }

    /// 构建历元时间
    CDate dateEpoch(stKepler.stEpoch.nYear,stKepler.stEpoch.uMonth,stKepler.stEpoch.uDay,
                    stKepler.stEpoch.uHour,stKepler.stEpoch.uMinute,
                    static_cast<double>(stKepler.stEpoch.uSecond)
                    +static_cast<double>(stKepler.stEpoch.uMSecond)*0.001);

    /// 开普勒六根数
    CVector vKepler(stKepler.dA,stKepler.dE,stKepler.dI,
                    stKepler.dRAAN,stKepler.dW,stKepler.dMA);

    double dMJDEpoch = dateEpoch.GetMJD(),
           dMJDStep = nStep*SECDAY,
           dMJDCal = dMJDStart;



    CVector vECI(6),vECF(6);

    PV tmpPV;

    /// 清空原有的数据
    stSatPos.stStart = stStartTime;
    stSatPos.stEnd = stEndTime;
    stSatPos.uStep = nStep;
    stSatPos.vJ2000.clear();
    stSatPos.vECF.clear();
    stSatPos.vTimes.clear();

    double dSpace = dMJDEnd - dMJDStart;
    long long nTimeLenMs = dSpace * DAYSEC * (long long) 1000;

    long long nStepMs = nStep * 1000;
    long long nCount = nTimeLenMs / nStepMs;

    if(nTimeLenMs % nStepMs != 0)
    {
        stSatPos.vJ2000.resize(nCount+2);
        stSatPos.vECF.resize(nCount+2);
        stSatPos.vTimes.resize(nCount+2);
    }
    else
    {
        stSatPos.vJ2000.resize(nCount+1);
        stSatPos.vECF.resize(nCount+1);
        stSatPos.vTimes.resize(nCount+1);
    }

    double dCalSec=0;
    for(int i=0;i<=nCount;++i,dMJDCal = dMJDStart + i *dMJDStep,dCalSec=nStep*i)
    {
        vECI = CKepler::State(GM_Earth,vKepler,dCalSec);
        tmpPV.stP.dX = vECI(0);
        tmpPV.stP.dY = vECI(1);
        tmpPV.stP.dZ = vECI(2);
        tmpPV.stV.dX = vECI(3);
        tmpPV.stV.dY = vECI(4);
        tmpPV.stV.dZ = vECI(5);

        CCoorSys::ECI2ECF(dMJDCal,vECI,vECF);

        stSatPos.vJ2000.push_back(tmpPV);
        stSatPos.vTimes.push_back(dMJDCal);

        tmpPV.stP.dX = vECF(0);
        tmpPV.stP.dY = vECF(1);
        tmpPV.stP.dZ = vECF(2);
        tmpPV.stV.dX = vECF(3);
        tmpPV.stV.dY = vECF(4);
        tmpPV.stV.dZ = vECF(5);
        stSatPos.vECF.push_back(tmpPV);
    }

    if(nTimeLenMs % nStepMs != 0)
    {
        vECI = CKepler::State(GM_Earth,vKepler,(dMJDEnd-dMJDEpoch)*DAYSEC);
        tmpPV.stP.dX = vECI(0);
        tmpPV.stP.dY = vECI(1);
        tmpPV.stP.dZ = vECI(2);
        tmpPV.stV.dX = vECI(3);
        tmpPV.stV.dY = vECI(4);
        tmpPV.stV.dZ = vECI(5);

        CCoorSys::ECI2ECF(dMJDEnd,vECI,vECF);

        stSatPos.vJ2000.push_back(tmpPV);
        stSatPos.vTimes.push_back(dMJDEnd);

        tmpPV.stP.dX = vECF(0);
        tmpPV.stP.dY = vECF(1);
        tmpPV.stP.dZ = vECF(2);
        tmpPV.stV.dX = vECF(3);
        tmpPV.stV.dY = vECF(4);
        tmpPV.stV.dZ = vECF(5);

        stSatPos.vECF.push_back(tmpPV);
    }

    return(true);
}

/// 根据两行星历计算轨道
bool SGP4(const BJTime &stStartTime, const BJTime &stEndTime, unsigned int nStep,
          const string& sLine1, const string& sLine2, SatellitePos &stSatPos)
{
    double dMJDStart,dMJDEnd;

    /// 判断时间是否有效
    if(!JudgeTimeValid(stStartTime,stEndTime,dMJDStart,dMJDEnd))
    {
        return(false);
    }

    CSGP4 calSGP4(sLine1,sLine2);
    if(!calSGP4)
    {
        return(false);
    }

    double dMJDStep = nStep*SECDAY,
           dMJDCal = dMJDStart;



    CVector vTEME(6),vECI(6),vECF(6);

    PV tmpPV;

    /// 清空原有的数据
    stSatPos.stStart = stStartTime;
    stSatPos.stEnd = stEndTime;
    stSatPos.uStep = nStep;
    stSatPos.vJ2000.clear();
    stSatPos.vECF.clear();
    stSatPos.vTimes.clear();

    double dSpace = dMJDEnd - dMJDStart;
    long long nTimeLenMs = dSpace * DAYSEC * (long long) 1000;

    long long nStepMs = nStep * 1000;
    long long nCount = nTimeLenMs / nStepMs;

    if(nTimeLenMs % nStepMs != 0)
    {
        stSatPos.vJ2000.resize(nCount+2);
        stSatPos.vECF.resize(nCount+2);
        stSatPos.vTimes.resize(nCount+2);
    }
    else
    {
        stSatPos.vJ2000.resize(nCount+1);
        stSatPos.vECF.resize(nCount+1);
        stSatPos.vTimes.resize(nCount+1);
    }

    for(int i=0;i<=nCount;++i,dMJDCal = dMJDStart + i *dMJDStep)
    {
        vTEME = calSGP4.CalPV(dMJDCal);

        CCoorSys::TEME2ECI(dMJDCal,vTEME,vECI);
        tmpPV.stP.dX = vECI(0);
        tmpPV.stP.dY = vECI(1);
        tmpPV.stP.dZ = vECI(2);
        tmpPV.stV.dX = vECI(3);
        tmpPV.stV.dY = vECI(4);
        tmpPV.stV.dZ = vECI(5);

        CCoorSys::TEME2ECF(dMJDCal,vTEME,vECF);

        stSatPos.vJ2000[i] = tmpPV;
        stSatPos.vTimes[i] = dMJDCal;

        tmpPV.stP.dX = vECF(0);
        tmpPV.stP.dY = vECF(1);
        tmpPV.stP.dZ = vECF(2);
        tmpPV.stV.dX = vECF(3);
        tmpPV.stV.dY = vECF(4);
        tmpPV.stV.dZ = vECF(5);

        stSatPos.vECF[i] = tmpPV;
    }

    if(nTimeLenMs % nStepMs != 0)
    {
        vTEME = calSGP4.CalPV(dMJDEnd);

        CCoorSys::TEME2ECI(dMJDEnd,vTEME,vECI);
        tmpPV.stP.dX = vECI(0);
        tmpPV.stP.dY = vECI(1);
        tmpPV.stP.dZ = vECI(2);
        tmpPV.stV.dX = vECI(3);
        tmpPV.stV.dY = vECI(4);
        tmpPV.stV.dZ = vECI(5);

        CCoorSys::TEME2ECF(dMJDEnd,vTEME,vECF);

        stSatPos.vJ2000[nCount+1] = tmpPV;
        stSatPos.vTimes[nCount+1] = dMJDEnd;

        tmpPV.stP.dX = vECF(0);
        tmpPV.stP.dY = vECF(1);
        tmpPV.stP.dZ = vECF(2);
        tmpPV.stV.dX = vECF(3);
        tmpPV.stV.dY = vECF(4);
        tmpPV.stV.dZ = vECF(5);

        stSatPos.vECF[nCount+1] = tmpPV;
    }

    return(true);
}

/// 生成卫星星座
vector<Satellite_Element> CreateConstellatory(Satellite_Element satTemplet,
                                              int nPlanes,
                                              int nNumSats)
{
    vector<Satellite_Element> vSatElement;
    vSatElement.reserve(nPlanes*nNumSats);

    double dRAANSpace = D2PI/nPlanes;
    double dMSpace = D2PI / nNumSats;

    Satellite_Element tmpSatellite=satTemplet;
    for(int i=0;i<nPlanes;++i)
    {
        for(int j=0; j<nNumSats; ++j)
        {
            tmpSatellite.stKepler.dRAAN = satTemplet.stKepler.dRAAN + i*dRAANSpace;
            tmpSatellite.stKepler.dMA = satTemplet.stKepler.dMA + j*dMSpace;

            while(tmpSatellite.stKepler.dRAAN > DPI)
            {
                tmpSatellite.stKepler.dRAAN -= D2PI;
            }

            while(tmpSatellite.stKepler.dRAAN < -DPI)
            {
                tmpSatellite.stKepler.dRAAN += D2PI;
            }

            while(tmpSatellite.stKepler.dMA > D2PI)
            {
                tmpSatellite.stKepler.dMA -= D2PI;
            }

            vSatElement.push_back(tmpSatellite);
        }
    }

    return(vSatElement);
}
