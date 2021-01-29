#include <VersionMathCommon.h>
#include <Math/VecMat.h>
#include <GisMath/GisMath.h>
#include <Satellite/Date.h>
#include <Satellite/CoorSys.h>
#include <Satellite/Kepler.h>
#include <Satellite/SGP4.h>
#include <SatelliteToolKit/SatelliteToolKit.h>
#include "../Inc/SatelliteToolKitCommon.h"

extern double iauAnp(double a);

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

    /// 清空原有的数据
    stSatPos.stStart = stStartTime;
    stSatPos.stEnd = stEndTime;
    stSatPos.uStepMs = nStep;
    stSatPos.vJ2000.clear();
    stSatPos.vECF.clear();
    stSatPos.vTimes.clear();

    double dSpace = dMJDEnd - dMJDStart;
    long long nTimeLenMs = dSpace * DAYSEC * (long long)1000;

    long long nStepMs = nStep * 1000;
    long long nCount = nTimeLenMs / nStepMs;

    if(nTimeLenMs % nStepMs != 0)
    {
        stSatPos.vJ2000.resize(nCount+2);
        stSatPos.vECF.resize(nCount+2);
        stSatPos.vTimes.resize(nCount+2);
        stSatPos.vLLA.resize(nCount+2);
    }
    else
    {
        stSatPos.vJ2000.resize(nCount+1);
        stSatPos.vECF.resize(nCount+1);
        stSatPos.vTimes.resize(nCount+1);
        stSatPos.vLLA.resize(nCount+1);
    }

    double dCalSec=0;
    for(int i=0;i<=nCount;++i,dMJDCal = dMJDStart + i *dMJDStep,dCalSec=nStep*i)
    {
        vECI = CKepler::State(GM_Earth,vKepler,dCalSec);
        stSatPos.vJ2000[i].stP.dX = vECI(0);
        stSatPos.vJ2000[i].stP.dY = vECI(1);
        stSatPos.vJ2000[i].stP.dZ = vECI(2);
        stSatPos.vJ2000[i].stV.dX = vECI(3);
        stSatPos.vJ2000[i].stV.dY = vECI(4);
        stSatPos.vJ2000[i].stV.dZ = vECI(5);

        CCoorSys::ECI2ECF(dMJDCal,vECI,vECF);

        stSatPos.vTimes[i] = dMJDCal;

        stSatPos.vECF[i].stP.dX = vECF(0);
        stSatPos.vECF[i].stP.dY = vECF(1);
        stSatPos.vECF[i].stP.dZ = vECF(2);
        stSatPos.vECF[i].stV.dX = vECF(3);
        stSatPos.vECF[i].stV.dY = vECF(4);
        stSatPos.vECF[i].stV.dZ = vECF(5);

        GisMath::XYZ2LBH(vECF(0),vECF(1),vECF(2),
                         stSatPos.vLLA[i].dX,stSatPos.vLLA[i].dY,stSatPos.vLLA[i].dZ);
        stSatPos.vLLA[i].dX *= DR2D;
        stSatPos.vLLA[i].dY *= DR2D;
    }

    if(nTimeLenMs % nStepMs != 0)
    {
        vECI = CKepler::State(GM_Earth,vKepler,(dMJDEnd-dMJDEpoch)*DAYSEC);
        stSatPos.vJ2000[nCount+1].stP.dX = vECI(0);
        stSatPos.vJ2000[nCount+1].stP.dY = vECI(1);
        stSatPos.vJ2000[nCount+1].stP.dZ = vECI(2);
        stSatPos.vJ2000[nCount+1].stV.dX = vECI(3);
        stSatPos.vJ2000[nCount+1].stV.dY = vECI(4);
        stSatPos.vJ2000[nCount+1].stV.dZ = vECI(5);

        CCoorSys::ECI2ECF(dMJDEnd,vECI,vECF);

        stSatPos.vTimes[nCount+1] = dMJDEnd;

        stSatPos.vECF[nCount+1].stP.dX = vECF(0);
        stSatPos.vECF[nCount+1].stP.dY = vECF(1);
        stSatPos.vECF[nCount+1].stP.dZ = vECF(2);
        stSatPos.vECF[nCount+1].stV.dX = vECF(3);
        stSatPos.vECF[nCount+1].stV.dY = vECF(4);
        stSatPos.vECF[nCount+1].stV.dZ = vECF(5);

        GisMath::XYZ2LBH(vECF(0),vECF(1),vECF(2),
                         stSatPos.vLLA[nCount+1].dX,stSatPos.vLLA[nCount+1].dY,stSatPos.vLLA[nCount+1].dZ);
        stSatPos.vLLA[nCount+1].dX *= DR2D;
        stSatPos.vLLA[nCount+1].dY *= DR2D;
    }

    return(true);
}

/// 根据两行星历计算轨道
bool SGP4(const BJTime &stStartTime, const BJTime &stEndTime, unsigned int nStepMs,
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

    double dMJDStep = (nStepMs/1000.0)*SECDAY,
           dMJDCal = dMJDStart;



    CVector vTEME(6),vECI(6),vECF(6);

    /// 清空原有的数据
    stSatPos.stStart = stStartTime;
    stSatPos.stEnd = stEndTime;
    stSatPos.uStepMs = nStepMs;
    stSatPos.vJ2000.clear();
    stSatPos.vECF.clear();
    stSatPos.vTimes.clear();

    double dSpace = dMJDEnd - dMJDStart;
    long long nTimeLenMs = dSpace * DAYSEC * (long long) 1000.0 + 0.5;

    //long long nStepMs = nStep * 1000;
    long long nCount = nTimeLenMs / nStepMs;

    if(nTimeLenMs % nStepMs != 0)
    {
        stSatPos.vJ2000.resize(nCount+2);
        stSatPos.vECF.resize(nCount+2);
        stSatPos.vTimes.resize(nCount+2);
        stSatPos.vLLA.resize(nCount+2);
    }
    else
    {
        stSatPos.vJ2000.resize(nCount+1);
        stSatPos.vECF.resize(nCount+1);
        stSatPos.vTimes.resize(nCount+1);
        stSatPos.vLLA.resize(nCount+1);
    }

    for(int i=0;i<=nCount;++i,dMJDCal = dMJDStart + i *dMJDStep)
    {
        vTEME = calSGP4.CalPV(dMJDCal);

        CCoorSys::TEME2ECI(dMJDCal,vTEME,vECI);
        stSatPos.vJ2000[i].stP.dX = vECI(0);
        stSatPos.vJ2000[i].stP.dY = vECI(1);
        stSatPos.vJ2000[i].stP.dZ = vECI(2);
        stSatPos.vJ2000[i].stV.dX = vECI(3);
        stSatPos.vJ2000[i].stV.dY = vECI(4);
        stSatPos.vJ2000[i].stV.dZ = vECI(5);

        CCoorSys::TEME2ECF(dMJDCal,vTEME,vECF);

        stSatPos.vTimes[i] = dMJDCal;

        stSatPos.vECF[i].stP.dX = vECF(0);
        stSatPos.vECF[i].stP.dY = vECF(1);
        stSatPos.vECF[i].stP.dZ = vECF(2);
        stSatPos.vECF[i].stV.dX = vECF(3);
        stSatPos.vECF[i].stV.dY = vECF(4);
        stSatPos.vECF[i].stV.dZ = vECF(5);

        GisMath::XYZ2LBH(vECF(0),vECF(1),vECF(2),
                         stSatPos.vLLA[i].dX,stSatPos.vLLA[i].dY,stSatPos.vLLA[i].dZ);
        stSatPos.vLLA[i].dX *= DR2D;
        stSatPos.vLLA[i].dY *= DR2D;
    }

    if(nTimeLenMs % nStepMs != 0)
    {
        vTEME = calSGP4.CalPV(dMJDEnd);

        CCoorSys::TEME2ECI(dMJDEnd,vTEME,vECI);
        stSatPos.vJ2000[nCount+1].stP.dX = vECI(0);
        stSatPos.vJ2000[nCount+1].stP.dY = vECI(1);
        stSatPos.vJ2000[nCount+1].stP.dZ = vECI(2);
        stSatPos.vJ2000[nCount+1].stV.dX = vECI(3);
        stSatPos.vJ2000[nCount+1].stV.dY = vECI(4);
        stSatPos.vJ2000[nCount+1].stV.dZ = vECI(5);

        CCoorSys::TEME2ECF(dMJDEnd,vTEME,vECF);

        stSatPos.vTimes[nCount+1] = dMJDEnd;

        stSatPos.vECF[nCount+1].stP.dX = vECF(0);
        stSatPos.vECF[nCount+1].stP.dY = vECF(1);
        stSatPos.vECF[nCount+1].stP.dZ = vECF(2);
        stSatPos.vECF[nCount+1].stV.dX = vECF(3);
        stSatPos.vECF[nCount+1].stV.dY = vECF(4);
        stSatPos.vECF[nCount+1].stV.dZ = vECF(5);

        GisMath::XYZ2LBH(vECF(0),vECF(1),vECF(2),
                         stSatPos.vLLA[nCount+1].dX,stSatPos.vLLA[nCount+1].dY,stSatPos.vLLA[nCount+1].dZ);
        stSatPos.vLLA[nCount+1].dX *= DR2D;
        stSatPos.vLLA[nCount+1].dY *= DR2D;
    }

    return(true);
}

/// 生成卫星星座
vector<Satellite_Element> CreateConstellatory(Satellite_Element satTemplet,
                                              int nPlanes,
                                              int nNumSats,
                                              int nFactor,
                                              double dRAANDelt)
{
    static vector<Satellite_Element> vSatElement;
    if(nFactor<0 || (nFactor+1 > nPlanes))
    {
        return(vSatElement);
    }
    double dDelt = nFactor * D2PI/(nPlanes*nNumSats);

    return(CreateConstellatoryBase(satTemplet,nPlanes,nNumSats,dDelt,dRAANDelt));
}

/// 生成卫星星座
vector<Satellite_Element> CreateConstellatoryBase(Satellite_Element satTemplet,
                                              int nPlanes,
                                              int nNumSats,
                                              double dDelt,
                                              double dRAANDelt)
{
    vector<Satellite_Element> vSatElement;

    double dRAANSpace = dRAANDelt/nPlanes;
    double dMSpace = D2PI / nNumSats;

    Math::CVector kepler(6);
    Satellite::CSGP4 tmpSGP4(nullptr,nullptr);
    switch(satTemplet.elemType)
    {
    case SAT_TLE:
        tmpSGP4.SetTLE(satTemplet.stTLE.sLine1,satTemplet.stTLE.sLine2);
        if(!tmpSGP4) return(vSatElement);
        kepler = tmpSGP4.ClassicalElements();
        break;
    case SAT_PV:
    {
        Math::CVector vPV(6);
        vPV(0) = satTemplet.stSatPV.stPV.stP.dX;
        vPV(1) = satTemplet.stSatPV.stPV.stP.dY;
        vPV(2) = satTemplet.stSatPV.stPV.stP.dZ;
        vPV(3) = satTemplet.stSatPV.stPV.stV.dX;
        vPV(4) = satTemplet.stSatPV.stPV.stV.dY;
        vPV(5) = satTemplet.stSatPV.stPV.stV.dZ;
        kepler = Satellite::CKepler::ClassicalElements(GM_Earth,vPV);
    }
        break;
    case SAT_TWOBODY:
        kepler(0) = satTemplet.stKepler.dA;
        kepler(1) = satTemplet.stKepler.dE;
        kepler(2) = satTemplet.stKepler.dI;
        kepler(3) = satTemplet.stKepler.dRAAN;
        kepler(4) = satTemplet.stKepler.dW;
        kepler(5) = satTemplet.stKepler.dMA;
        break;
    }

    vSatElement.reserve(nPlanes*nNumSats);
    Math::CVector tempKepler(kepler);

    Satellite_Element tmpSatellite = satTemplet;
    char tleOut[2][73];
    for(int i=0;i<nPlanes;++i)
    {
        for(int j=0; j<nNumSats; ++j)
        {
            tempKepler(3) = kepler(3) + i*dRAANSpace;
            tempKepler(5) = kepler(5) + j*dMSpace + i*dDelt;

            tempKepler(3) = iauAnp(tempKepler(3));
            tempKepler(5) = iauAnp(tempKepler(5));

            switch (tmpSatellite.elemType)
            {
            case SAT_TLE:
                Satellite::CKepler::Classical2TLE(tempKepler,tmpSGP4.GetTLEEpoch()+DJM0,99,tleOut);
                tmpSatellite.stTLE.sLine1=std::string(tleOut[0],73);
                tmpSatellite.stTLE.sLine1=std::string(tleOut[1],73);
                break;
            case SAT_PV:
            {
                Math::CVector vPV = Satellite::CKepler::State(GM_Earth,tempKepler);
                tmpSatellite.stSatPV.stPV.stP.dX=vPV(0);
                tmpSatellite.stSatPV.stPV.stP.dY=vPV(1);
                tmpSatellite.stSatPV.stPV.stP.dZ=vPV(2);
                tmpSatellite.stSatPV.stPV.stV.dX=vPV(3);
                tmpSatellite.stSatPV.stPV.stV.dY=vPV(4);
                tmpSatellite.stSatPV.stPV.stV.dZ=vPV(5);
            }
                break;
            case SAT_TWOBODY:
            {
                tmpSatellite.stKepler.dA=tempKepler(0);
                tmpSatellite.stKepler.dE=tempKepler(1);
                tmpSatellite.stKepler.dI=tempKepler(2);
                tmpSatellite.stKepler.dRAAN=tempKepler(3);
                tmpSatellite.stKepler.dW=tempKepler(4);
                tmpSatellite.stKepler.dMA=tempKepler(5);
            }
                break;
            }
            vSatElement.push_back(tmpSatellite);
        }
    }

    return(vSatElement);
}
