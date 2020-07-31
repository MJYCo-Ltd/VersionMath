#include <Date.h>
#include <VecMat.h>
#include <sofam.h>
#include <YPRAngle.h>
#include <JPLEphemeris.h>
#include <TimeSys.h>
#include <CoorSys.h>
#include <GisMath.h>

#include "../Inc/SatelliteToolKitCommon.h"

using namespace Aerospace;
using namespace Math;

/// 判断时间是否有效
bool JudgeTimeValid(const BJTime &stStartTime, const BJTime &stEndTime,
                    double& dMJDStart, double& dMJDEnd)
{
    CDate dateStart(stStartTime.nYear,stStartTime.uMonth,stStartTime.uDay,
               stStartTime.uHour,stStartTime.uMinute,
               static_cast<double>(stStartTime.uSecond)
               +static_cast<double>(stStartTime.uMSecond)*0.001);

    CDate dateEnd(stEndTime.nYear,stEndTime.uMonth,stEndTime.uDay,
               stEndTime.uHour,stEndTime.uMinute,
               static_cast<double>(stEndTime.uSecond)
               +static_cast<double>(stEndTime.uMSecond)*0.001);

    dMJDStart = dateStart.GetMJD();
    dMJDEnd = dateEnd.GetMJD();

    if(dMJDStart > dMJDEnd)
    {
        return(false);
    }

    return(true);
}

/// 判断数据是否有效
int  JudgeDataTrend(const vector<TimeElev>& vElev,float dMin)
{
    return(0);
}

/// 判断点是与载荷平面相交
bool JudgeIsInsert(const PV& satECFPV, const Pos& station3DPos, const Pos& satPRY,
                   RotateType eRotate, CVector& vInsert, double& dAngle)
{
    if(InsertEarth(satECFPV.stP,station3DPos))
    {
        return(false);
    }

    /// 构建卫星本体到全局的转换矩阵
    CVector Pos(satECFPV.stP.dX,satECFPV.stP.dY,satECFPV.stP.dZ);

    CMatrix tmpMatrix = CalSatRoteMatrix(satECFPV,satPRY,eRotate);

    CVector StationPos(station3DPos.dX,station3DPos.dY,station3DPos.dZ);
    CVector vGlobal = StationPos - Pos;

    /// 将全局的坐标转成局部坐标
    vInsert = tmpMatrix * vGlobal;

    CVector Center(0,0,1);
    double dDot = CVecMat::Dot(vInsert,Center);

    /// 如果与面平行
    if(fabs(dDot) < DF_ZERO || dDot < 0)
    {
        return(false);
    }
    else
    {
        dAngle = acos(dDot/vInsert.Length());
        /// 计算与z=1平面相交的点
        double dScale = 1.0 / vInsert.GetZ();
        vInsert *= dScale;
        return(true);
    }
}

/// 判断两点是否在地球之外
bool InsertEarth(const Pos &rPos1, const Pos &rPos2)
{
    CVector vR1(rPos1.dX,rPos1.dY,rPos1.dZ);
    CVector vR2(rPos2.dX,rPos2.dY,rPos2.dZ);
    return InsertEarth(vR1,vR2);
}

bool InsertEarth(const CVector &vR1, const CVector &vR2)
{
    double dR1 = vR1.Length();

    double dDot = CVecMat::Dot(vR1,vR2);

    ///两者夹角
    double dB   = acos(dDot/dR1/vR2.Length());
    double dA   = acos(R_Earth/dR1);

    if(dA > dB)
    {
        return(false);
    }

    CVector vR3 = vR1 - vR2;

    CVector vCross = CVecMat::Cross(vR1,vR2);
    double dS = vCross.Length();
    double dC = vR3.Length();

    double dR = dS / dC;

    return(dR < R_Earth);
}
/// 计算旋转矩阵
CMatrix CalSatRoteMatrix(const PV& satPV, const Pos& satPRY, RotateType eRotate)
{
    CMatrix tmpMatrix = CalSatMatrix(satPV);

    /// 如果姿态角没有设置
    if(fabs(satPRY.dX) < DF_ZERO && fabs(satPRY.dY) < DF_ZERO
      && fabs(satPRY.dZ) < DF_ZERO)
    {
        return(tmpMatrix);
    }


    CMatrix matRotate = CalRotateMatrix(satPRY,eRotate);


    return(matRotate*tmpMatrix);
}

/// 计算太阳高度角
double CalSolarAltitude(double dMJD,const Pos& gruondStationGeo)
{
    double dAltitude(-1);
    /// 初始化JPL星历数据
    CJPLEphemeris* pJPL = CJPLEphemeris::GetInstance();

    if(nullptr == pJPL || !pJPL->IsInit())
    {
        return(dAltitude);
    }

    CTimeSys timeSys(dMJD);
    CVector vECISun = pJPL->GetSunPos(timeSys.GetTT());

    CVector vECFSun;
    CCoorSys::ECI2ECF(dMJD,vECISun,vECFSun);

    CVector vStation3D(3);
    GisMath::LBH2XYZ(gruondStationGeo.dX,gruondStationGeo.dY,gruondStationGeo.dZ,
                      vStation3D(0),vStation3D(1),vStation3D(2));

    /// 计算太阳相对于地面站位置
    CVector vGlobal = vECFSun - vStation3D;

    CVector vLocal(3);

    /// 转成局部坐标
    GisMath::GLOBAL2LOCAL(gruondStationGeo.dX,gruondStationGeo.dY,vGlobal,vLocal);

    double dAzim,dElev;
    CVecMat::AzEl(vLocal,dAzim,dElev);

    return(dElev);
}

/// 计算叉乘
double CalCross(const STK_Point &rP1, const STK_Point &rP2, const STK_Point &rP)
{
    return((rP2.dX - rP1.dX) * (rP.dY - rP1.dY) -(rP.dX - rP1.dX) * (rP2.dY - rP1.dY));
}

/// 计算卫星旋转矩阵
CMatrix CalSatMatrix(const PV &satPV)
{
    /// 构建卫星本体到全局的转换矩阵
    CVector Pos(satPV.stP.dX,satPV.stP.dY,satPV.stP.dZ);
    CVector rZ = -Pos;
    CVector V(satPV.stV.dX,satPV.stV.dY,satPV.stV.dZ);

    CVector rY = CVecMat::Cross(rZ , V);
    CVector rX = CVecMat::Cross(rY,rZ);

    rX.Normalize();
    rY.Normalize();
    rZ.Normalize();

    CMatrix tmpMatrix(3,3);
    tmpMatrix.SetRow(0,rX);
    tmpMatrix.SetRow(1,rY);
    tmpMatrix.SetRow(2,rZ);
    return(tmpMatrix);
}

CMatrix CalSatMatrix(const CVector& vPV)
{
    /// 构建卫星本体到全局的转换矩阵
    CVector Pos(vPV.slice(0,2));
    CVector rZ = -Pos;
    CVector V(vPV.slice(3,5));

    CVector rY = CVecMat::Cross(rZ , V);
    CVector rX = CVecMat::Cross(rY,rZ);

    rX.Normalize();
    rY.Normalize();
    rZ.Normalize();

    CMatrix tmpMatrix(3,3);
    tmpMatrix.SetRow(0,rX);
    tmpMatrix.SetRow(1,rY);
    tmpMatrix.SetRow(2,rZ);
    return(tmpMatrix);
}

CMatrix CalRotateMatrix(const Pos &satPRY, RotateType eRotate)
{
    switch (eRotate)
    {
    case Rota_PRY:
    case Rota_YXZ:
        return(CYPRAngle::CreateMatrix(satPRY.dX,satPRY.dY,satPRY.dZ,PRY));
    case Rota_PYR:
    case Rota_YZX:
        return(CYPRAngle::CreateMatrix(satPRY.dX,satPRY.dY,satPRY.dZ,PYR));
    case Rota_RYP:
    case Rota_XZY:
        return(CYPRAngle::CreateMatrix(satPRY.dX,satPRY.dY,satPRY.dZ,RYP));
    case Rota_RPY:
    case Rota_XYZ:
        return(CYPRAngle::CreateMatrix(satPRY.dX,satPRY.dY,satPRY.dZ,RPY));
    case Rota_YPR:
    case Rota_ZYX:
        return(CYPRAngle::CreateMatrix(satPRY.dX,satPRY.dY,satPRY.dZ,YPR));
    case Rota_YRP:
    case Rota_ZXY:
        return(CYPRAngle::CreateMatrix(satPRY.dX,satPRY.dY,satPRY.dZ,YRP));
    }
	// 默认返回PRY
	return(CYPRAngle::CreateMatrix(satPRY.dX, satPRY.dY, satPRY.dZ, PRY));
}

/// 判断载荷和地球是否相交
bool CanIntersection(const Pos &satPos, const CVector &rDir, const double &dMaxAngle)
{
    CVector vCenter(0,0,1);

    /// 两者夹角
    double dθ = acos(CVecMat::Dot(rDir,vCenter));

    CVector vSatPos(satPos.dX,satPos.dY,satPos.dZ);

    double dSatθ =  asin(R_Earth/vSatPos.Length());


    /// 判断是否相交
    return((dSatθ + dMaxAngle) >= dθ);
}

bool EllipseVisible(const double &dVisibleHAngle, const double &dVisibleVAngle, const double &dAngle, const CVector &vInsert)
{
    /// 如果是圆
    if(fabs(dVisibleVAngle - dVisibleHAngle) < 1e-10)
    {
        if(dVisibleVAngle >= dAngle)
        {
            return(true);
        }
        else
        {
            return (false);
        }
    }

    /// 如果是椭圆
    double dMinAngle = dVisibleVAngle > dVisibleHAngle ? dVisibleHAngle : dVisibleVAngle;

    if(dMinAngle > dAngle)
    {
        return(true);
    }
    else
    {
        double dY= tan(dVisibleVAngle),dX=tan(dVisibleHAngle);
        dY = dY * dY;
        dX = dX * dX;

        /// x*x/dX/dX + y*y/dY/dY > 1则在椭圆外
        double dValue = vInsert.GetX()*vInsert.GetX()/dX + vInsert.GetY() * vInsert.GetY()/dY;
        return(dValue <= 1);
    }
}

bool RectangleVisible(const double &dVisibleHAngle, const double &dVisibleVAngle, const double &dAngle, const CVector &vInsert)
{
    double dMinAngle = dVisibleVAngle > dVisibleHAngle ? dVisibleHAngle : dVisibleVAngle;

    /// 如果最小的角都大于对地的角
    if(dMinAngle > dAngle)
    {
        return(true);
    }
    else
    {
        double dY= tan(dVisibleVAngle),dX=tan(dVisibleHAngle);
        STK_Point point1={-dX,-dY},
                  point2={dX,-dY},
                  point3={dX,dY},
                  point4={-dX,dY},
                  pointP={vInsert.GetX(),vInsert.GetY()};
        return CalCross(point1,point2,pointP) * CalCross(point3,point4,pointP) >= 0
                && CalCross(point2,point3,pointP) * CalCross(point4,point1,pointP) >= 0;
    }
}
