#include <Math/Matrix.h>
#include <Math/Quaternion.h>
#include <Math/EulerAngle.h>
#include <sofa.h>

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
    double dTemp[3][3];
    iauIr(dTemp);
    switch(type)
    {
    case E121:
        iauRx(dA,dTemp);
        iauRy(dB,dTemp);
        iauRx(dC,dTemp);
        break;
    case E123:
        iauRx(dA,dTemp);
        iauRy(dB,dTemp);
        iauRz(dC,dTemp);
        break;
    case E131:
        iauRx(dA,dTemp);
        iauRz(dB,dTemp);
        iauRx(dC,dTemp);
        break;
    case E132:
        iauRx(dA,dTemp);
        iauRz(dB,dTemp);
        iauRy(dC,dTemp);
        break;
    case E212:
        iauRy(dA,dTemp);
        iauRx(dB,dTemp);
        iauRy(dC,dTemp);
        break;
    case E213:
        iauRy(dA,dTemp);
        iauRx(dB,dTemp);
        iauRz(dC,dTemp);
        break;
    case E231:
        iauRy(dA,dTemp);
        iauRz(dB,dTemp);
        iauRx(dC,dTemp);
        break;
    case E232:
        iauRy(dA,dTemp);
        iauRz(dB,dTemp);
        iauRy(dC,dTemp);
        break;
    case E312:
        iauRz(dA,dTemp);
        iauRx(dB,dTemp);
        iauRy(dC,dTemp);
        break;
    case E313:
        iauRz(dA,dTemp);
        iauRx(dB,dTemp);
        iauRz(dC,dTemp);
        break;
    case E321:
        iauRz(dA,dTemp);
        iauRy(dB,dTemp);
        iauRx(dC,dTemp);
        break;
    case E323:
        iauRz(dA,dTemp);
        iauRy(dB,dTemp);
        iauRz(dC,dTemp);
        break;
    default:
        break;
    }

    CMatrix tmpMatrix(dTemp,3,3);
    return(tmpMatrix);
}
