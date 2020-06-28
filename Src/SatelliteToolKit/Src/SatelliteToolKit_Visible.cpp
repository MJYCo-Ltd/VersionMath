#include <VecMat.h>
#include <GisMath.h>
#include <CoorSys.h>
#include "SatelliteToolKit.h"
#include "../Inc/SatelliteToolKitCommon.h"

using namespace Math;
using namespace Aerospace;

/// 地面站对卫星的可见
bool IsVisible(const Pos& satECFPos, const Pos& station3D, double dVisibleAngle)
{
    double dL,dB,dH;
    GisMath::XYZ2LBH(station3D.dX, station3D.dY, station3D.dZ, dL, dB, dH);

    CVector vGlobal(satECFPos.dX-station3D.dX,satECFPos.dY-station3D.dY,satECFPos.dZ-station3D.dZ)
            ,vLocal(3);

    GisMath::GLOBAL2LOCAL(dL,dB,vGlobal,vLocal);

    double dAzim,dElev;
    CVecMat::AzEl(vLocal,dAzim,dElev);

    if(dElev > dVisibleAngle)
    {
        return(true);
    }

    return(false);
}

/// 判断两个星是否可见
bool IsVisible(const PV& satECIPV, const PV& satoTherECIPV, double dConeAngle)
{
    if(InsertEarth(satECIPV.stP,satoTherECIPV.stP))
    {
        return(false);
    }

    CMatrix tmpMatrix = CalSatMatrix(satECIPV);
    /// 构建卫星本体到全局的转换矩阵
    CVector Pos(satECIPV.stP.dX,satECIPV.stP.dY,satECIPV.stP.dZ);

    CVector StationPos(satoTherECIPV.stP.dX,satoTherECIPV.stP.dY,satoTherECIPV.stP.dZ);
    CVector vGlobal = StationPos - Pos;
    /// 将全局的坐标转成局部坐标
    CVector vInsert = tmpMatrix * vGlobal;
    /// 将全局的坐标转成局部坐标
    vInsert = tmpMatrix * vGlobal;

    CVector Center(0,0,1);
    double dDot = CVecMat::Dot(vInsert,Center);
    double dAngle = acos(dDot/vInsert.Length());

    return(dConeAngle >= dAngle);
}

/// 卫星对地面目标的可见
bool EllipseVisible(const PV& satECFPV, const Pos& station3DPos, const Pos& satPRY,
               double dVisibleHAngle, double dVisibleVAngle, RotateType eRotate)
{
    double dAngle;
    CVector vInsert;

    if(JudgeIsInsert(satECFPV,station3DPos,satPRY,eRotate,vInsert,dAngle))
    {
        return(EllipseVisible(dVisibleHAngle,dVisibleVAngle,dAngle,vInsert));
    }

    return(false);
}

bool RectangleVisible(const PV& satECFPV, const Pos& station3DPos, const Pos& satPRY,
                      double dVisibleHAngle, double dVisibleVAngle, RotateType eRotate)
{
    double dAngle;
    CVector vInsert;

    if(JudgeIsInsert(satECFPV,station3DPos,satPRY,eRotate,vInsert,dAngle))
    {
        return(RectangleVisible(dVisibleHAngle,dVisibleVAngle,dAngle,vInsert));
    }
    return(false);
}
