#include "JPLEphemeris.h"
#include "SatelliteToolKit.h"
#include "../Inc/SatelliteToolKitCommon.h"
#include "sofam.h"

const double C_SMaxAngle(DPI*0.5*0.99);
using namespace Aerospace;


vector<Pos> Intersect(const PV& satPV,
                      const CMatrix& rotateMatrix,
                      const vector<Pos>& vAngle)
{
    vector<Pos> vResult;

    /// 将位置点旋转到卫星本体坐标系下
    vector<Pos> vSatPos;
    Pos tmpPos;
    CVector satPos(3),localPos(3);
    int nSize = vAngle.size();
    for(int nIndex=0; nIndex < nSize; ++nIndex)
    {
        const Pos& one = vAngle[nIndex];
        localPos.Set(one.dX,one.dY,one.dZ);

        satPos = localPos * rotateMatrix;

        tmpPos.dX = satPos.GetX();
        tmpPos.dY = satPos.GetY();
        tmpPos.dZ = satPos.GetZ();

        vSatPos.push_back(tmpPos);
    }



    return(vResult);
}

/// 计算圆与地球相交的区域
vector<Pos> IntersectCircle(const PV& satPV,
                            const Pos& satPRY,
                            RotateType eRotate,
                            double dHAngle,
                            double dVAngle)
{
    vector<Pos> vResult;

    /// 计算旋转矩阵
    CMatrix rotateMatrix = CalRotateMatrix(satPRY,eRotate);
    CVector vDir(0.,0.,1.);
    vDir = vDir * rotateMatrix;

    double dMaxAngle = dHAngle > dVAngle ? dHAngle : dVAngle;
    if(!CanIntersection(satPV.stP,vDir,dMaxAngle))
    {
        return(vResult);
    }

    //// x=a cosθ　 y=b sinθ
    double dA = tan(dHAngle);
    double dB = tan(dVAngle);

    vector<Pos> vAllCircle;
    Pos tmpPos;
    tmpPos.dZ = 1;

    for(double dStart=0; dStart < D2PI; dStart += DD2R)
    {
        tmpPos.dX = dA*cos(dStart);
        tmpPos.dY = dB*cos(dStart);
        vAllCircle.push_back(tmpPos);
    }

    return(Intersect(satPV,rotateMatrix,vAllCircle));
}


vector<Pos> IntersectRectangle(const PV& satPV,
                               const Pos& satPRY,
                               RotateType eRotate,
                               double dHAngle,
                               double dVAngle)
{
    vector<Pos> vResult;

    /// 计算旋转矩阵
    CMatrix rotateMatrix = CalRotateMatrix(satPRY,eRotate);
    CVector vDir(0.,0.,1.);
    vDir = vDir * rotateMatrix;

    double dMaxAngle(0);
    if(dHAngle > C_SMaxAngle || dVAngle > C_SMaxAngle)
    {
        dMaxAngle = DPI * 0.5;
    }
    else
    {
        double dH = tan(dHAngle);
        double dV = tan(dVAngle);

        dMaxAngle = atan(sqrt(dH*dH + dV*dV));
    }

    /// 判断两者是否相交
    if(!CanIntersection(satPV.stP,vDir,dMaxAngle))
    {
        return(vResult);
    }

    return(vResult);
}
