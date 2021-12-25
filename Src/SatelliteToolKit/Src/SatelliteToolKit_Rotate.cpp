#include <cmath>
#include <Math/VecMat.h>
#include <Math/Quaternion.h>
#include <Math/YPRAngle.h>
#include <Satellite/SGP4.h>
#include <Satellite/Kepler.h>
#include <Satellite/Date.h>
#include <Satellite/CoorSys.h>
#include <SatelliteToolKit/SatelliteToolKit.h>
#include <VersionMathCommon.h>
#include "../Inc/SatelliteToolKitCommon.h"

/// 计算姿态
bool CalPRY(const PV& satPV,const Pos& rPos,RotateType eRoteType,Pos& satPRY)
{
    if(InsertEarth(satPV.stP,rPos))
    {
        return(false);
    }

    CMatrix tmpMatrix = CalSatMatrix(satPV);
    CVector globalPos(rPos.dX - satPV.stP.dX,rPos.dY - satPV.stP.dY,rPos.dZ - satPV.stP.dZ);
    CVector localPos = tmpMatrix * globalPos;

    CVector vZ(0,0,1);
    localPos.Normalize();

    /// 根据叉乘 计算法向量
    CVector vAxis = CVecMat::Cross(vZ,localPos);

    /// 根据点乘 计算旋转角度
    double dDot = CVecMat::Dot(vZ,localPos);
    double dAngle =acos(dDot);
    CQuaternion quatation(vAxis,dAngle);
    CMatrix rotateMatrix = quatation.GetMatrix();

    YPR_Rotate tmpRotate;
    switch (eRoteType)
    {
    case Rota_PRY:
    case Rota_YXZ:
        CYPRAngle::CalTransform(rotateMatrix,PRY,tmpRotate);
        break;
    case Rota_PYR:
    case Rota_YZX:
        CYPRAngle::CalTransform(rotateMatrix,PYR,tmpRotate);
        break;
    case Rota_RYP:
    case Rota_XZY:
        CYPRAngle::CalTransform(rotateMatrix,RYP,tmpRotate);
        break;
    case Rota_RPY:
    case Rota_XYZ:
        CYPRAngle::CalTransform(rotateMatrix,RPY,tmpRotate);
        break;
    case Rota_YPR:
    case Rota_ZYX:
        CYPRAngle::CalTransform(rotateMatrix,YPR,tmpRotate);
        break;
    case Rota_YRP:
    case Rota_ZXY:
        CYPRAngle::CalTransform(rotateMatrix,YRP,tmpRotate);
        break;

    }

    satPRY.dX = tmpRotate.dX;
    satPRY.dY = tmpRotate.dY;
    satPRY.dZ = tmpRotate.dZ;

    return(true);
}

/// 计算相交时刻的卫星轨道
Kepler InsertOribit(const Satellite_Element& originOribit,double dAngle,const BJTime& stInsertTime)
{
    Kepler tmpKepler;
    PV satPV;

    /// 构建历元时间
    Aerospace::CDate dateNow(stInsertTime.nYear,stInsertTime.uMonth,stInsertTime.uDay,
                    stInsertTime.uHour,stInsertTime.uMinute,
                    static_cast<double>(stInsertTime.uSecond)
                    +static_cast<double>(stInsertTime.uMSecond)*0.001);

    double dMJD = dateNow.GetMJD();
    Math::CVector vECI;


    switch(originOribit.elemType)
    {
    case SAT_TLE:
    {
        Satellite::CSGP4 calSGP4(originOribit.stTLE.sLine1,originOribit.stTLE.sLine2);
        if(!calSGP4)
        {
            return(tmpKepler);
        }
        CCoorSys::TEME2J2000(dMJD,calSGP4.CalPV(dMJD),vECI);
        satPV.stP.dX = vECI(0);
        satPV.stP.dY = vECI(1);
        satPV.stP.dZ = vECI(2);

        satPV.stV.dX = vECI(3);
        satPV.stV.dY = vECI(4);
        satPV.stV.dZ = vECI(5);
    }
        break;
    case SAT_TWOBODY:
    {
        Aerospace::CDate dateEpoch(originOribit.stKepler.stEpoch.nYear,originOribit.stKepler.stEpoch.uMonth,
                                   originOribit.stKepler.stEpoch.uDay,originOribit.stKepler.stEpoch.uHour,
                                   originOribit.stKepler.stEpoch.uMinute,
                        static_cast<double>(originOribit.stKepler.stEpoch.uSecond)
                        +static_cast<double>(originOribit.stKepler.stEpoch.uMSecond)*0.001);
        CVector vKepler(originOribit.stKepler.dA,originOribit.stKepler.dE,originOribit.stKepler.dI,
                        originOribit.stKepler.dRAAN,originOribit.stKepler.dW,originOribit.stKepler.dMA);
        vECI = Satellite::CKepler::State(GM_Earth,vKepler,(dMJD-dateEpoch.GetMJD())*DAYSEC);
        satPV.stP.dX = vECI(0);
        satPV.stP.dY = vECI(1);
        satPV.stP.dZ = vECI(2);

        satPV.stV.dX = vECI(3);
        satPV.stV.dY = vECI(4);
        satPV.stV.dZ = vECI(5);
    }
        break;
    case SAT_PV:
    {
        Aerospace::CDate dateEpoch(originOribit.stSatPV.stEpoch.nYear,originOribit.stSatPV.stEpoch.uMonth,
                                   originOribit.stSatPV.stEpoch.uDay,originOribit.stSatPV.stEpoch.uHour,
                                   originOribit.stSatPV.stEpoch.uMinute,
                        static_cast<double>(originOribit.stSatPV.stEpoch.uSecond)
                        +static_cast<double>(originOribit.stSatPV.stEpoch.uMSecond)*0.001);

        CVector vPV(originOribit.stSatPV.stPV.stP.dX,originOribit.stSatPV.stPV.stP.dY,
                    originOribit.stSatPV.stPV.stP.dZ,originOribit.stSatPV.stPV.stV.dX,
                    originOribit.stSatPV.stPV.stV.dY,originOribit.stSatPV.stPV.stV.dZ);
        vECI = Satellite::CKepler::State(GM_Earth,Satellite::CKepler::ClassicalElements(
                                             GM_Earth, vPV),(dMJD-dateEpoch.GetMJD())*DAYSEC);
        satPV.stP.dX = vECI(0);
        satPV.stP.dY = vECI(1);
        satPV.stP.dZ = vECI(2);

        satPV.stV.dX = vECI(3);
        satPV.stV.dY = vECI(4);
        satPV.stV.dZ = vECI(5);
    }
        break;
    }

    Math::CVector vLocal = Math::CVector(1,0,0) * Math::CVecMat::R_z(dAngle);
    /// 计算位置
    Math::CVector vGlobal = vLocal * CalSatMatrix(satPV);
    vGlobal.Normalize();
    Math::CVector vNewV = vECI.slice(3,5).Length() * vGlobal;
    vECI(3) = vNewV(0);
    vECI(4) = vNewV(1);
    vECI(5) = vNewV(2);

    CVector vKepler = Satellite::CKepler::ClassicalElements(GM_Earth, vECI);
    tmpKepler.stEpoch = stInsertTime;
    tmpKepler.dA = vKepler(0);
    tmpKepler.dE = vKepler(1);
    tmpKepler.dI = vKepler(2);
    tmpKepler.dRAAN = vKepler(3);
    tmpKepler.dW = vKepler(4);
    tmpKepler.dMA = vKepler(5);

    return(tmpKepler);
}
