#include "JPLEphemeris.h"
#include "SatelliteToolKit.h"
#include "../Inc/SatelliteToolKitCommon.h"
#include "sofam.h"

const double C_SMaxAngle(DPI*0.5*0.99);
using namespace Aerospace;

/// 判断射线与椭球相交
bool CalLineInterEllipsoid(const Pos& pt,const Pos& stDir,Pos& rInsertPos)
{
    /// 椭球参数
    double a  = R_Earth;
    double b  = R_Earth2;
    CVector vDir(stDir.dX,stDir.dY,stDir.dZ);
    vDir.Normalize();

    double a2 = a*a;
    double b2 = b*b;
    double dx = vDir.GetX();
    double dy = vDir.GetY();
    double dz = vDir.GetZ();
    double ex = pt.dX;
    double ey = pt.dY;
    double ez = pt.dZ;

    double A   = b2*dx*dx + b2*dy*dy + a2*dz*dz;
    double B   = 2 * (b2*dx*ex + b2*dy*ey+ a2*dz*ez);
    double C   = b2*ex*ex + b2*ey*ey + a2*ez*ez - a2*b2;

    /// 求该视线与椭球的两个交点
    double delta = B*B-4*A*C;
    if(delta<0) return false; //没有交点

    double deltaS = sqrt(delta);
    double t1 = (-B + deltaS)/(2*A);
    double t2 = (-B - deltaS)/(2*A);

    /// 如果t1,t2都小于0，则认为视线与椭球没有交点
    if(t1<0 && t2<0)
    {
        return(false);
    }

    /// 取t1, t2中绝对值较小的那个
    double t = fabs(t1)<fabs(t2)? t1:t2;

    ///  根据直线参数方程：X = Dt + E，可得直线与椭球的交点坐标
    rInsertPos.dX = t * dx + ex;
    rInsertPos.dY = t * dy + ey;
    rInsertPos.dZ = t * dz + ez;

    return(true);
}

vector<Pos> Intersect(const PV& satPV,
                      const CMatrix& rotateMatrix,
                      const vector<Pos>& vAngle)
{
    vector<Pos> vResult;

    /// 将位置点旋转到卫星本体坐标系下
    vector<Pos> vSatPos;
    Pos tmpPos;
    CVector satPos(3),localPos(3);
    for(auto one : vAngle)
    {
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
