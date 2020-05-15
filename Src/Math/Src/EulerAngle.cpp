#include "Matrix.h"
#include "Quaternion.h"
#include "EulerAngle.h"

using namespace Math;

CEulerAngle::CEulerAngle(const CMatrix& rMatrix,EULERROTATE type):
    m_emOriginalType(type),m_matRote(rMatrix)
{
}

CEulerAngle::CEulerAngle(const Euler_Rotate &stRotate, EULERROTATE type):
    m_emOriginalType(type),m_stRotate(stRotate)
{
    m_matRote = CreateMatrix(m_stRotate.dA,m_stRotate.dB,m_stRotate.dC,m_emOriginalType);
}

CEulerAngle::CEulerAngle(double dA, double dB, double dC, EULERROTATE type):
    m_emOriginalType(type)
{
    m_stRotate.dA = dA;
    m_stRotate.dB = dB;
    m_stRotate.dC = dC;

    m_matRote = CreateMatrix(m_stRotate.dA,m_stRotate.dB,m_stRotate.dC,m_emOriginalType);
}

CEulerAngle::~CEulerAngle()
{
}

/// 通过角A、角B、角C角度值，以及旋转类型，计算旋转矩阵
CMatrix CEulerAngle::CreateMatrix(double dA, double dB, double dC, EULERROTATE type)
{
    switch(type)
    {
    case E121:
        return((CQuaternion(CVector(1,0,0),dC)*CQuaternion(CVector(0,1,0),dB)*
                CQuaternion(CVector(1,0,0),dA)).GetMatrix());
    case E123:
        return((CQuaternion(CVector(0,0,1),dC)*CQuaternion(CVector(0,1,0),dB)*
                CQuaternion(CVector(1,0,0),dA)).GetMatrix());
    case E131:
        return((CQuaternion(CVector(1,0,0),dC)*CQuaternion(CVector(0,0,1),dB)*
                CQuaternion(CVector(1,0,0),dA)).GetMatrix());
    case E132:
        return((CQuaternion(CVector(0,1,0),dC)*CQuaternion(CVector(0,0,1),dB)*
                CQuaternion(CVector(1,0,0),dA)).GetMatrix());
    case E212:
        return((CQuaternion(CVector(0,1,0),dC)*CQuaternion(CVector(1,0,0),dB)*
                CQuaternion(0,1,0,dA)).GetMatrix());
    case E213:
        return((CQuaternion(CVector(0,0,1),dC)*CQuaternion(CVector(1,0,0),dB)*
                CQuaternion(CVector(0,1,0),dA)).GetMatrix());
    case E231:
        return((CQuaternion(CVector(1,0,0),dC)*CQuaternion(CVector(0,0,1),dB)*
                CQuaternion(CVector(0,1,0),dA)).GetMatrix());
    case E232:
        return((CQuaternion(CVector(0,1,0),dC)*CQuaternion(CVector(0,0,1),dB)*
                CQuaternion(CVector(0,1,0),dA)).GetMatrix());
    case E312:
        return((CQuaternion(CVector(0,1,0),dC)*CQuaternion(CVector(1,0,0),dB)*
                CQuaternion(CVector(0,0,1),dA)).GetMatrix());
    case E313:
        return((CQuaternion(CVector(0,0,1),dC)*CQuaternion(CVector(1,0,0),dB)*
                CQuaternion(CVector(0,0,1),dA)).GetMatrix());
    case E321:
        return((CQuaternion(CVector(1,0,0),dC)*CQuaternion(CVector(0,1,0),dB)*
                CQuaternion(CVector(0,0,1),dA)).GetMatrix());
    case E323:
        return((CQuaternion(CVector(0,0,1),dC)*CQuaternion(CVector(0,1,0),dB)*
                CQuaternion(CVector(0,0,1),dA)).GetMatrix());
    default:
        break;
    }

    CMatrix tmpMatrix;
    return(tmpMatrix);
}
