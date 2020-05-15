#include <VecMat.h>
#include <Quaternion.h>
#include <YPRAngle.h>
#include "SatelliteToolKit.h"
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
