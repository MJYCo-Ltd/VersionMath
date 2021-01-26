#include <iomanip>
#include <float.h>

#include <VersionMathCommon.h>
#include <Math/VecMat.h>
#include <GisMath/GisMath.h>
#include <Satellite/TimeSys.h>
#include <Satellite/CoorSys.h>
#include <Satellite/SGP4.h>
#include <Satellite/Kepler.h>
#include <SatelliteToolKit/SatelliteToolKit.h>
#include "../Inc/SatelliteToolKitCommon.h"

using namespace Aerospace;
using namespace Math;
using namespace Satellite;

double CalSolarAltitude(const BJTime& BJTime,const Pos&  stPos)
{
    CDate bjTime(BJTime.nYear,BJTime.uMonth,BJTime.uDay
                 ,BJTime.uHour,BJTime.uMinute,static_cast<double>(BJTime.uSecond)+BJTime.uMSecond*1e-3);

    return(CalSolarAltitude(bjTime.GetMJD(),stPos));
}

/// 计算满足要求的太阳高度角的时间段
vector<Period> SolarAltitude(const BJTime& stStartTime, const BJTime& stEndTime,
                             float fIncidentMinAngle, float fIncidentMaxAngle, const Pos& stPos)
{
    vector<Period> vResult;

    double dMJDStart,dMJDEnd;

    /// 时间不对直接返回
    if(!JudgeTimeValid(stStartTime,stEndTime,dMJDStart,dMJDEnd))
    {
        return(vResult);
    }

    /// 时间步长
    double dStep = 60. * SECDAY;

    CVector vECISun,vECFSun,vLocal,vStation(stPos.dX,stPos.dY,stPos.dZ);
    CVector vGlobal;

    vector<TimeElev> vElev;
    TimeElev tmpTimeElev;

    /// 遍历所有时间
    for(double dMJDTime = dMJDStart; dMJDTime < dMJDEnd; dMJDTime += dStep)
    {
        tmpTimeElev.dElev = CalSolarAltitude(dMJDTime,stPos);
        tmpTimeElev.dMJDTime = dMJDTime;

        /// 将角度保存到数组中
        vElev.push_back(tmpTimeElev);
    }

    /// 判断数据是否有效
    JudgeDataTrend(vElev,fIncidentMinAngle);

    return(vResult);
}

struct bestPoint
{
    double dEnterTime;
    double dOutTime;
    bool    bIsTowPos;
};

struct minAngle
{
    double dTime;
    double dAngle;
};

typedef bool (*pFun)(const PV& satECFPV, const Pos& station3DPos, const Pos& satPRY,
                     double dVisibleHAngle, double dVisibleVAngle, RotateType eRotate);

typedef bool (*pJudgeFun)(const double& dVisibleHAngle, const double& dVisibleVAngle, const double& dAngle, const CVector& vInsert);
/// 地面站
vector<Period> VisiblePeriod(const BJTime& stStartTime,
                             const BJTime& stEndTime,
                             const Satellite_Element& stSatElement,
                             const Pos&   stPos,
                             RotateType emRotate,
                             const Pos& satPRY,
                             double  dHAngle,
                             double  dVAngle,
                             ShapeType emType)
{
    pFun CallFun;
    pJudgeFun CallJudge;
    switch (emType)
    {
    case eEllipse:
        CallFun = EllipseVisible;
        CallJudge = EllipseVisible;
        break;
    case eRectangle:
        CallFun = RectangleVisible;
        CallJudge = RectangleVisible;
        break;
    }

    vector<Period> vResult;

    double dMJDStart,dMJDEnd;

    /// 判断开始时间和结束时间
    if(!JudgeTimeValid(stStartTime,stEndTime,dMJDStart,dMJDEnd))
    {
        return(vResult);
    }

    CSGP4 tmpSGP4(stSatElement.stTLE.sLine1,stSatElement.stTLE.sLine2);
    if(!tmpSGP4)
    {
        return(vResult);
    }

    CVector vKepler = tmpSGP4.ClassicalElements();

    if(!vKepler)
    {
        return(vResult);
    }

    double dT = CKepler::T(vKepler(0), GM_Earth) * SECDAY;

    CVector vGround3D(3);
    GisMath::LBH2XYZ(stPos.dX*DD2R,stPos.dY*DD2R,0,vGround3D(0),vGround3D(1),vGround3D(2));

    /// 计算步长
    double dStep = dT/360.;

    double dMJDStep = dStep,
            dMJDCal = dMJDStart;

    double dSpace = dMJDEnd - dMJDStart;
    long long nTimeLenMs = dSpace * DAYSEC * 1000;

    long long nStepMs = dStep * 1000 * DAYSEC;
    long long nCount = nTimeLenMs / nStepMs;

    CVector vTEME(6),vSatPos(3),vTEMEGround(3),vCenter(0,0,1);

    vector<bestPoint> vBestPos;
    vector<minAngle>    vMinAngle;
    int nReserveCount = dSpace / dT * 1.5;
    vBestPos.reserve(nReserveCount);
    vMinAngle.reserve(nReserveCount);
    minAngle tmpMinAngle;

    double dCosAngle,dDot,dLength,dMinAngle(DBL_MAX);
    bool   bEntered(false),bCalRotate(false),bBiger(false),bCanPush(false);

    CMatrix matRotate;

    bestPoint tmpBestPos;

    /// 如果姿态角没有设置
    if(fabs(satPRY.dX) < DF_ZERO && fabs(satPRY.dY) < DF_ZERO
            && fabs(satPRY.dZ) < DF_ZERO)
    {
        bCalRotate = false;
    }
    else
    {
        matRotate = CalRotateMatrix(satPRY,emRotate);
        bCalRotate = true;
    }

    for(int i=0;i<=nCount;++i,dMJDCal = dMJDStart + i *dMJDStep)
    {
        vTEME = tmpSGP4.CalPV(dMJDCal);
        /// 只考虑地球自转
        vSatPos = vTEME.slice(0,2);
        vTEMEGround = vGround3D * CCoorSys::TEME2ECF(dMJDCal);

        /// 排除不可见的点
        if(InsertEarth(vSatPos,vTEMEGround))
        {
            if(bEntered)
            {
                tmpBestPos.bIsTowPos = false;
                vBestPos.push_back(tmpBestPos);
                bEntered =false;
            }

            if(bCanPush)
            {
                vMinAngle.push_back(tmpMinAngle);
                bCanPush = false;
                dMinAngle = DBL_MAX;
            }
        }
        else
        {
            /// 目标相对于卫星的位置
            vSatPos = vTEMEGround - vSatPos;

            /// 计算旋转后的地面目标
            if(bCalRotate)
            {
                vSatPos = matRotate * CalSatMatrix(vTEME)*vSatPos;
            }
            else
            {
                vSatPos = CalSatMatrix(vTEME) * vSatPos;
            }

            dDot = CVecMat::Dot(vSatPos,vCenter);

            /// 如果与面平行
            if(fabs(dDot) < DF_ZERO || dDot < 0)
            {
                if(bEntered)
                {
                    tmpBestPos.bIsTowPos = false;
                    vBestPos.push_back(tmpBestPos);
                    bEntered =false;
                }

                if(bCanPush)
                {
                    vMinAngle.push_back(tmpMinAngle);
                    bCanPush = false;
                    dMinAngle = DBL_MAX;
                }
            }
            else
            {
                dLength = vSatPos.Length();
                dCosAngle = acos(dDot/dLength);

                /// 寻找最小值
                if(dMinAngle > dCosAngle)
                {
                    dMinAngle = dCosAngle;
                    tmpMinAngle.dTime = dMJDCal;
                    tmpMinAngle.dAngle = dMinAngle;
                    bCanPush = true;
                }
                else
                {
                    bBiger=true;
                }

                dLength = 1.0 / vSatPos.GetZ();
                vSatPos *= dLength;
                /// 找到小值
                if(CallJudge(dHAngle,dVAngle,dCosAngle,vSatPos))
                {
                    if(!bEntered)
                    {
                        tmpBestPos.dEnterTime = dMJDCal;
                        tmpBestPos.bIsTowPos = false;
                        bEntered = true;
                    }
                    else
                    {
                        tmpBestPos.dOutTime = dMJDCal;
                        tmpBestPos.bIsTowPos = true;
                    }

                    if(nCount == i)
                    {
                        vBestPos.push_back(tmpBestPos);
                    }
                }
                else if(bEntered)
                {
                    vBestPos.push_back(tmpBestPos);

                    tmpBestPos.bIsTowPos = false;
                    bEntered = false;
                }
                else if(bCanPush && bBiger)
                {
                    vMinAngle.push_back(tmpMinAngle);
                    bCanPush = false;
                    bBiger = false;
                }
            }
        }
    }

    /// 计算详细相交时刻
    CVector vECF(6);
    PV enterPV;
    Pos ground3D={vGround3D(0),vGround3D(1),vGround3D(2)};

    double dMinStep = SECDAY * 1e-3,dEnterStep;


    double dSecond,dInsertTime;
    CVector vMinAngleSatPos(3);

    /// 遍历所有的
    int nIndexMinAngle(0),nSizeMinAngle(vMinAngle.size());
    for(;nIndexMinAngle<nSizeMinAngle;++nIndexMinAngle)
    {
        minAngle& one = vMinAngle[nIndexMinAngle];
        /// 不可见时间点
        dInsertTime = one.dTime;
        dMinAngle = dInsertTime;
        dEnterStep = dMJDStep*0.5;
        dMJDCal = dInsertTime-dEnterStep;

        for(dEnterStep*=0.5;
            dEnterStep>dMinStep;
            dEnterStep*=0.5)
        {
            vTEME = tmpSGP4.CalPV(dMJDCal);
            /// 只考虑地球自转
            vSatPos = vTEME.slice(0,2);
            vTEMEGround = vGround3D * CCoorSys::TEME2ECF(dMJDCal);

            if(!InsertEarth(vSatPos,vTEMEGround))
            {
                /// 目标相对于卫星的位置
                vSatPos = vTEMEGround - vSatPos;

                /// 计算旋转后的地面目标
                if(bCalRotate)
                {
                    vSatPos = matRotate * CalSatMatrix(vTEME)*vSatPos;
                }
                else
                {
                    vSatPos = CalSatMatrix(vTEME) * vSatPos;
                }

                dDot = CVecMat::Dot(vSatPos,vCenter);
                if(fabs(dDot) < DF_ZERO || dDot < 0)
                {
                    dMJDCal += dEnterStep;
                }
                else
                {
                    dLength = vSatPos.Length();
                    dCosAngle = acos(dDot/dLength);
                    if(dCosAngle > one.dAngle)
                    {
                        dMJDCal += dEnterStep;
                    }
                    else
                    {
                        vMinAngleSatPos = vSatPos;
                        dLength = 1.0 / vSatPos.GetZ();
                        vMinAngleSatPos *= dLength;
                        dMinAngle = dMJDCal;
                        one.dAngle = dCosAngle;
                        dMJDCal -= dEnterStep;
                    }
                }
            }
            else
            {
                dMJDCal += dEnterStep;
            }
        }

        dInsertTime = one.dTime;
        dEnterStep = dMJDStep*0.5;
        dMJDCal = dInsertTime+dEnterStep;

        for(dEnterStep*=0.5;
            dEnterStep>dMinStep;
            dEnterStep*=0.5)
        {
            vTEME = tmpSGP4.CalPV(dMJDCal);
            /// 只考虑地球自转
            vSatPos = vTEME.slice(0,2);
            vTEMEGround = vGround3D * CCoorSys::TEME2ECF(dMJDCal);

            if(!InsertEarth(vSatPos,vTEMEGround))
            {
                /// 目标相对于卫星的位置
                vSatPos = vTEMEGround - vSatPos;

                /// 计算旋转后的地面目标
                if(bCalRotate)
                {
                    vSatPos = matRotate * CalSatMatrix(vTEME)*vSatPos;
                }
                else
                {
                    vSatPos = CalSatMatrix(vTEME) * vSatPos;
                }

                dDot = CVecMat::Dot(vSatPos,vCenter);
                if(fabs(dDot) < DF_ZERO || dDot < 0)
                {
                    dMJDCal -= dEnterStep;
                }
                else
                {
                    dLength = vSatPos.Length();
                    dCosAngle = acos(dDot/dLength);
                    if(dCosAngle > one.dAngle)
                    {
                        dMJDCal -= dEnterStep;
                    }
                    else
                    {
                        vMinAngleSatPos = vSatPos;
                        dLength = 1.0 / vSatPos.GetZ();
                        vMinAngleSatPos *= dLength;

                        dMinAngle = dMJDCal;
                        one.dAngle = dCosAngle;
                        dMJDCal += dEnterStep;
                    }
                }
            }
            else
            {
                dMJDCal -= dEnterStep;
            }
        }

        /// 最小角度值差别大于1毫秒
        if(fabs(dMinAngle - one.dTime)>dMinStep && CallJudge(dHAngle,dVAngle,one.dAngle,vMinAngleSatPos))
        {
            tmpBestPos.dEnterTime = dMinAngle;
            tmpBestPos.bIsTowPos = false;
            vBestPos.push_back(tmpBestPos);
        }

    }


    Period tmpPeriod;
    int nYear,nMonth,nDay,nHour,nMinute;

    CDate tmpDate(UTC);
    int nResultCount=vBestPos.size();
    vResult.resize(nResultCount);
    int nIndex(0);

    --nResultCount;
    for(nIndexMinAngle=0,nSizeMinAngle=vBestPos.size();nIndexMinAngle<nSizeMinAngle;++nIndexMinAngle)
    {
        bestPoint& one = vBestPos[nIndexMinAngle];
        if(one.bIsTowPos)
        {
            dEnterStep = dMJDStep*0.5;
            dInsertTime = one.dEnterTime;
            dMJDCal = dInsertTime - dEnterStep;

            for(dEnterStep*=0.5;
                dEnterStep>dMinStep;
                dEnterStep*=0.5)
            {
                vTEME = tmpSGP4.CalPV(dMJDCal);
                CCoorSys::TEME2ECF(dMJDCal,vTEME,vECF);
                enterPV.stP.dX = vECF(0);
                enterPV.stP.dY = vECF(1);
                enterPV.stP.dZ = vECF(2);
                enterPV.stV.dX = vECF(3);
                enterPV.stV.dY = vECF(4);
                enterPV.stV.dZ = vECF(5);

                if(CallFun(enterPV,ground3D,satPRY,dHAngle,dVAngle,emRotate))
                {
                    dInsertTime = dMJDCal;
                    dMJDCal -= dEnterStep;
                }
                else
                {
                    dMJDCal += dEnterStep;
                }
            }

            if(nIndex == nResultCount)
            {
                if(dInsertTime > dMJDEnd)
                {
                    vResult.erase(vResult.end()-1);
                    break;
                }
            }

            /// 生成输出数据
            tmpDate.SetJD(DJM0+dInsertTime);
            tmpDate.GetDate(nYear,nMonth,nDay,nHour,nMinute,dSecond);
            BJTime& rBJTime = vResult[nIndex].stStart;
            rBJTime.nYear = nYear;
            rBJTime.uMonth = nMonth;
            rBJTime.uDay = nDay;
            rBJTime.uHour = nHour;
            rBJTime.uMinute = nMinute;
            rBJTime.uSecond = dSecond;
            rBJTime.uMSecond=(dSecond-rBJTime.uSecond)*1000;


            vResult[nIndex].vLLa.resize(4);
            /// 构建卫星本体到全局的转换矩阵
            CVector Pos(enterPV.stP.dX,enterPV.stP.dY,enterPV.stP.dZ);

            CMatrix tmpMatrix = CalSatRoteMatrix(enterPV,satPRY,emRotate);

            CVector vLocal(0,tan(dVAngle),1);
            CVector vGlobal(3);

            /// 将局部的坐标转成全局坐标
            vGlobal = vLocal * tmpMatrix;

            if(GisMath::CalLineInterEarth(Pos,vGlobal,vLocal))
            {
                GisMath::XYZ2LBH(vLocal,vGlobal);
                vResult[nIndex].vLLa[0].dX = vGlobal(0)*DR2D;
                vResult[nIndex].vLLa[0].dY = vGlobal(1)*DR2D;
                vResult[nIndex].vLLa[0].dZ = vGlobal(2);
            }

            vLocal.Set(0,-tan(dVAngle),1);
            /// 将全局的坐标转成局部坐标
            vGlobal = vLocal * tmpMatrix;

            if(GisMath::CalLineInterEarth(Pos,vGlobal,vLocal))
            {
                GisMath::XYZ2LBH(vLocal,vGlobal);
                vResult[nIndex].vLLa[3].dX = vGlobal(0)*DR2D;
                vResult[nIndex].vLLa[3].dY = vGlobal(1)*DR2D;
                vResult[nIndex].vLLa[3].dZ = vGlobal(2);
            }


            dSecond = dInsertTime;


            dEnterStep = dMJDStep*0.5;
            dInsertTime = one.dOutTime;
            dMJDCal = dInsertTime + dEnterStep;

            for(dEnterStep*=0.5;
                dEnterStep>dMinStep;
                dEnterStep*=0.5)
            {
                vTEME = tmpSGP4.CalPV(dMJDCal);
                CCoorSys::TEME2ECF(dMJDCal,vTEME,vECF);
                enterPV.stP.dX = vECF(0);
                enterPV.stP.dY = vECF(1);
                enterPV.stP.dZ = vECF(2);
                enterPV.stV.dX = vECF(3);
                enterPV.stV.dY = vECF(4);
                enterPV.stV.dZ = vECF(5);

                if(CallFun(enterPV,ground3D,satPRY,dHAngle,dVAngle,emRotate))
                {
                    dInsertTime = dMJDCal;
                    dMJDCal += dEnterStep;
                }
                else
                {
                    dMJDCal -= dEnterStep;
                }
            }

            if(nIndex == nResultCount)
            {
                if(dInsertTime > dMJDEnd)
                {
                    dInsertTime = dMJDEnd;
                }
            }

            vResult[nIndex].dDurationTime = (dInsertTime-dSecond)*DAYSEC;
            /// 生成输出数据
            tmpDate.SetJD(DJM0+dInsertTime);
            tmpDate.GetDate(nYear,nMonth,nDay,nHour,nMinute,dSecond);
            BJTime& rEndBJTime = vResult[nIndex].stEnd;
            rEndBJTime.nYear = nYear;
            rEndBJTime.uMonth = nMonth;
            rEndBJTime.uDay = nDay;
            rEndBJTime.uHour = nHour;
            rEndBJTime.uMinute = nMinute;
            rEndBJTime.uSecond = dSecond;
            rEndBJTime.uMSecond=(dSecond-rEndBJTime.uSecond)*1000;

            /// 构建卫星本体到全局的转换矩阵
            Pos.Set(enterPV.stP.dX,enterPV.stP.dY,enterPV.stP.dZ);

            tmpMatrix = CalSatRoteMatrix(enterPV,satPRY,emRotate);

            vLocal.Set(0,tan(dVAngle),1);

            /// 将局部的坐标转成全局坐标
            vGlobal = vLocal * tmpMatrix;

            if(GisMath::CalLineInterEarth(Pos,vGlobal,vLocal))
            {
                GisMath::XYZ2LBH(vLocal,vGlobal);
                vResult[nIndex].vLLa[1].dX = vGlobal(0)*DR2D;
                vResult[nIndex].vLLa[1].dY = vGlobal(1)*DR2D;
                vResult[nIndex].vLLa[1].dZ = vGlobal(2);
            }

            vLocal.Set(0,-tan(dVAngle),1);
            /// 将全局的坐标转成局部坐标
            vGlobal = vLocal * tmpMatrix;

            if(GisMath::CalLineInterEarth(Pos,vGlobal,vLocal))
            {
                GisMath::XYZ2LBH(vLocal,vGlobal);
                vResult[nIndex].vLLa[2].dX = vGlobal(0)*DR2D;
                vResult[nIndex].vLLa[2].dY = vGlobal(1)*DR2D;
                vResult[nIndex].vLLa[2].dZ = vGlobal(2);
            }
        }
        else
        {
            dEnterStep = dMJDStep*0.5;
            dInsertTime = one.dEnterTime;
            dMJDCal = dInsertTime-dEnterStep;

            for(dEnterStep*=0.5;
                dEnterStep>dMinStep;
                dEnterStep*=0.5)
            {
                vTEME = tmpSGP4.CalPV(dMJDCal);
                CCoorSys::TEME2ECF(dMJDCal,vTEME,vECF);
                enterPV.stP.dX = vECF(0);
                enterPV.stP.dY = vECF(1);
                enterPV.stP.dZ = vECF(2);
                enterPV.stV.dX = vECF(3);
                enterPV.stV.dY = vECF(4);
                enterPV.stV.dZ = vECF(5);

                if(CallFun(enterPV,ground3D,satPRY,dHAngle,dVAngle,emRotate))
                {
                    dInsertTime = dMJDCal;
                    dMJDCal -= dEnterStep;
                }
                else
                {
                    dMJDCal += dEnterStep;
                }
            }

            if(nIndex == nResultCount)
            {
                if(dInsertTime > dMJDEnd)
                {
                    vResult.erase(vResult.end()-1);
                    break;
                }
            }

            /// 生成输出数据
            tmpDate.SetJD(DJM0+dInsertTime);
            tmpDate.GetDate(nYear,nMonth,nDay,nHour,nMinute,dSecond);
            BJTime& rBJTime = vResult[nIndex].stStart;
            rBJTime.nYear = nYear;
            rBJTime.uMonth = nMonth;
            rBJTime.uDay = nDay;
            rBJTime.uHour = nHour;
            rBJTime.uMinute = nMinute;
            rBJTime.uSecond = dSecond;
            rBJTime.uMSecond=(dSecond-rBJTime.uSecond)*1000;

            vResult[nIndex].vLLa.resize(4);
            /// 构建卫星本体到全局的转换矩阵
            CVector Pos(enterPV.stP.dX,enterPV.stP.dY,enterPV.stP.dZ);

            CMatrix tmpMatrix = CalSatRoteMatrix(enterPV,satPRY,emRotate);

            CVector vLocal(0,tan(dVAngle),1);
            CVector vGlobal(3);

            /// 将局部的坐标转成全局坐标
            vGlobal = vLocal * tmpMatrix;

            if(GisMath::CalLineInterEarth(Pos,vGlobal,vLocal))
            {
                GisMath::XYZ2LBH(vLocal,vGlobal);
                vResult[nIndex].vLLa[0].dX = vGlobal(0)*DR2D;
                vResult[nIndex].vLLa[0].dY = vGlobal(1)*DR2D;
                vResult[nIndex].vLLa[0].dZ = vGlobal(2);
            }

            vLocal.Set(0,-tan(dVAngle),1);
            /// 将全局的坐标转成局部坐标
            vGlobal = vLocal * tmpMatrix;

            if(GisMath::CalLineInterEarth(Pos,vGlobal,vLocal))
            {
                GisMath::XYZ2LBH(vLocal,vGlobal);
                vResult[nIndex].vLLa[3].dX = vGlobal(0)*DR2D;
                vResult[nIndex].vLLa[3].dY = vGlobal(1)*DR2D;
                vResult[nIndex].vLLa[3].dZ = vGlobal(2);
            }

            dSecond = dInsertTime;


            dEnterStep = dMJDStep*0.5;
            dInsertTime = one.dEnterTime;
            dMJDCal = dInsertTime+dEnterStep;

            for(dEnterStep*=0.5;
                dEnterStep>dMinStep;
                dEnterStep*=0.5)
            {
                vTEME = tmpSGP4.CalPV(dMJDCal);
                CCoorSys::TEME2ECF(dMJDCal,vTEME,vECF);
                enterPV.stP.dX = vECF(0);
                enterPV.stP.dY = vECF(1);
                enterPV.stP.dZ = vECF(2);
                enterPV.stV.dX = vECF(3);
                enterPV.stV.dY = vECF(4);
                enterPV.stV.dZ = vECF(5);

                if(CallFun(enterPV,ground3D,satPRY,dHAngle,dVAngle,emRotate))
                {
                    dInsertTime = dMJDCal;
                    dMJDCal += dEnterStep;
                }
                else
                {
                    dMJDCal -= dEnterStep;
                }
            }

            if(nIndex == nResultCount)
            {
                if(dInsertTime > dMJDEnd)
                {
                    dInsertTime = dMJDEnd;
                }
            }

            vResult[nIndex].dDurationTime = (dInsertTime-dSecond)*DAYSEC;
            /// 生成输出数据
            tmpDate.SetJD(DJM0+dInsertTime);
            tmpDate.GetDate(nYear,nMonth,nDay,nHour,nMinute,dSecond);
            BJTime& rEndBJTime = vResult[nIndex].stEnd;
            rEndBJTime.nYear = nYear;
            rEndBJTime.uMonth = nMonth;
            rEndBJTime.uDay = nDay;
            rEndBJTime.uHour = nHour;
            rEndBJTime.uMinute = nMinute;
            rEndBJTime.uSecond = dSecond;
            rEndBJTime.uMSecond=(dSecond-rEndBJTime.uSecond)*1000;

            /// 构建卫星本体到全局的转换矩阵
            Pos.Set(enterPV.stP.dX,enterPV.stP.dY,enterPV.stP.dZ);

            tmpMatrix = CalSatRoteMatrix(enterPV,satPRY,emRotate);

            vLocal.Set(0,tan(dVAngle),1);

            /// 将局部的坐标转成全局坐标
            vGlobal = vLocal * tmpMatrix;

            if(GisMath::CalLineInterEarth(Pos,vGlobal,vLocal))
            {
                GisMath::XYZ2LBH(vLocal,vGlobal);
                vResult[nIndex].vLLa[1].dX = vGlobal(0)*DR2D;
                vResult[nIndex].vLLa[1].dY = vGlobal(1)*DR2D;
                vResult[nIndex].vLLa[1].dZ = vGlobal(2);
            }

            vLocal.Set(0,-tan(dVAngle),1);
            /// 将全局的坐标转成局部坐标
            vGlobal = vLocal * tmpMatrix;

            if(GisMath::CalLineInterEarth(Pos,vGlobal,vLocal))
            {
                GisMath::XYZ2LBH(vLocal,vGlobal);
                vResult[nIndex].vLLa[2].dX = vGlobal(0)*DR2D;
                vResult[nIndex].vLLa[2].dY = vGlobal(1)*DR2D;
                vResult[nIndex].vLLa[2].dZ = vGlobal(2);
            }
        }

        ++nIndex;
    }
    return(vResult);
}

vector<Period> VisiblePeriod(const SatellitePos& stSatPosEye,
                             const SatellitePos& stSatPosTarget,
                             float fVisibleAngle)
{
    vector<Period> vResult;

    return(vResult);
}

vector<Period> IntersectCircle(const SatellitePos& stSatPos,const Pos& stPos,float fAngle)
{
    vector<Period> vResult;

    return(vResult);
}

vector<Period> IntersectRectangle(const SatellitePos& stSatPos,const Pos& stPos,
                                  float fHAngle, float fVAngle)
{
    vector<Period> vResult;

    return(vResult);
}
