#include <cmath>
#include <limits>
#include <VersionMathCommon.h>
#include <Math/VecMat.h>
#include <Math/YPRAngle.h>
#include <Math/Quaternion.h>
#include <sofa.h>

using namespace std;

/************************************************************************/
/*         1   0     0           +cosy 0 -siny        +cosz +sinz 0     */
/*  R(x) = 0 +cosx +sinx  R(y) =  0    1   0   R(z) = -sinz +cosz 0     */
/*         0 -sinx +cosx         +siny 0 +cosy          0    0    1     */
/************************************************************************/
/*********************************************************************************
               +cosy        0     -siny                +cosz     +sinz         0
    R(x,y) = +sinx*siny  +cosx +sinx*cosy  R(x,z) = -cosx*sinz  +cosx*cosz  +sinx
             +cosx*siny  -sinx +cosx*cosy           +sinx*sinz  -sinx*cosz  +cosx
//////////////////////////////////////////////////////////////////////////////////////
            +cosy   +sinx*siny  -cosx*siny          +cosy*cosz  +cosy*sinz  -siny
    R(y,x) =  0       +cosx       +sinx    R(y,z) =   -sinz        +cosz      0
            +siny  -sinx*cosy   +cosx*cosy          +siny*cosz  +siny*sinz  +cosy
//////////////////////////////////////////////////////////////////////////////////////
             +cosz +cosx*sinz  +sinx*sinz           +cosy*cosz   +sinz  -siny*cosz
    R(z,x) = -sinz +cosx*cosz  +sinx*cosz  R(z,y) = -cosy*sinz   +cosz  +siny*sinz
              0      -sinx        +cosx               +siny        0      +cosy
************************************************************************************/
/***************************************************************************************************************************************************************
                    +cosy*cosz                 +cosy*sinz               -siny                   +cosy*cosz                 +sinz         -siny*cosz
    R(x,y,z) = +sinx*siny*cosz-cosx*sinz  +sinx*siny*sinz+cosx*cosz   +sinx*cosy   R(x,z,y) =  -cosx*cosy*sinz+sinx*siny  +cosx*cosz +cosx*siny*sinz+sinx*cosy
               +cosx*siny*cosz+sinx*sinz  +cosx*siny*sinz-sinx*cosz   +cosx*cosy               +sinx*cosy*sinz+cosx*siny  -sinx*cosz -sinx*siny*sinz+cosx*cosy
****************************************************************************************************************************************************************
               +cosy*cosz-sinx*siny*sinz    +cosy*sinz+sinx*siny*cosz  -cosx*siny               +cosy*cosz  +cosx*cosy*sinz+sinx*siny +sinx*cosy*sinz-cosx*siny
    R(y,x,z) =      -cosx*sinz                      +cosx*cosz           +sinx     R(y,z,x) =      -sinz       +cosx*cosz                   +sinx*cosz
               +siny*cosz+sinx*cosy*sinz    +siny*sinz-sinx*cosy*cosz  +cosx*cosy               +siny*cosz  +cosx*siny*cosz-sinx*cosy +sinx*siny*sinz+cosx*cosy
****************************************************************************************************************************************************************
               +cosy*cosz+sinx*siny*sinz   +cosx*sinz  -siny*cosz+sinx*cosy*sinz               +cosy*cosz   +cosx*sinz+sinx*siny*cosz +sinx*sinz-cosx*siny*cosz
    R(z,x,y) = -cosy*sinz+cosx*siny*cosz   +cosx*cosz  +siny*sinz+sinx*cosy*cosz   R(z,y,x) =  -cosy*sinz   +cosx*cosz-sinx*siny*sinz +sinx*cosz+cosx*siny*sinz
                     +cosx*siny               -sinx            +cosx*cosy                         +siny      -sinx*cosy         +cosx*cosy
*****************************************************************************************************************************************************************/
/************************************************************************
  asin 取值范围 -Pi/2 到 Pi/2
  acos 取值范围  0    到 Pi
  atan 取值范围 -Pi/2 到 Pi/2
  atan2取值范围 -Pi   到 Pi
*************************************************************************/
static const double DBL_ZERO = numeric_limits<double>::min();
static const double PI_2     = DPI*0.5;
using namespace Math;
CYPRAngle::CYPRAngle(const CMatrix& rMatrix, YPRROTATE type):
    m_emOriginalType(type),m_matRote(rMatrix)
{
    CalTransform(m_matRote,m_emOriginalType,m_stRotate);
}

CYPRAngle::CYPRAngle(double dRoll, double dPitch, double dYaw, YPRROTATE type):m_emOriginalType(type)
{
    m_stRotate.dRoll = dRoll;
    m_stRotate.dPitch = dPitch;
    m_stRotate.dYaw = dYaw;

    m_matRote = CreateMatrix(dRoll,dPitch,dYaw,m_emOriginalType);
}

CYPRAngle::CYPRAngle(const YPR_Rotate &stRotate, YPRROTATE type):m_emOriginalType(type),m_stRotate(stRotate)
{
    m_matRote = CreateMatrix(stRotate.dRoll,stRotate.dPitch,stRotate.dYaw,m_emOriginalType);
}

CYPRAngle::~CYPRAngle()
{

}

/// 计算旋转角度
bool CYPRAngle::CalTransform(YPRROTATE OrigType, YPRROTATE PurType, const YPR_Rotate &OrigRotate, YPR_Rotate &PurRotate)
{
    //如果坐标类型相同，则直接返回
    if(OrigType == PurType)
    {
        PurRotate = OrigRotate;
        return(true);
    }
    //如果坐标类型不同，则返回计算完成的数据
    else
    {
        CMatrix matrxRotate = CreateMatrix(OrigRotate.dX,OrigRotate.dY,OrigRotate.dZ,OrigType);
        CalTransform(matrxRotate,PurType,PurRotate);

        return (true);
    }
}

/// 根据矩阵计算旋转角度
void CYPRAngle::CalTransform(const CMatrix &rotateMatrix, YPRROTATE type, YPR_Rotate &stRotate)
{
    switch(type)
    {
    /**********************************************************************************
                      +cosy*cosz                 +cosy*sinz               -siny
      R(x,y,z) = +sinx*siny*cosz-cosx*sinz  +sinx*siny*sinz+cosx*cosz   +sinx*cosy
                 +cosx*siny*cosz+sinx*sinz  +cosx*siny*sinz-sinx*cosz   +cosx*cosy
     **********************************************************************************/
    case XYZ:
    case RPY:
        if(fabs(rotateMatrix(0,0))<DBL_ZERO && fabs(rotateMatrix(0,1))<DBL_ZERO)
        {
            stRotate.dX = DPI;

            if(rotateMatrix(0,2) < 0.0) stRotate.dY = PI_2;
            else                        stRotate.dY = -PI_2;

            stRotate.dZ = atan2(rotateMatrix(1,0),-rotateMatrix(1,1));
        }
        else
        {
            stRotate.dX = atan2(rotateMatrix(1,2),rotateMatrix(2,2));
            stRotate.dY = atan2(-rotateMatrix(0,2),sqrt(rotateMatrix(1,2)*rotateMatrix(1,2)+rotateMatrix(2,2)*rotateMatrix(2,2)));
            stRotate.dZ = atan2(rotateMatrix(0,1),rotateMatrix(0,0));
        }
        break;
        /****************************************************************************
                   +cosy*cosz                 +sinz         -siny*cosz
      R(x,z,y) =  -cosx*cosy*sinz+sinx*siny  +cosx*cosz +cosx*siny*sinz+sinx*cosy
                  +sinx*cosy*sinz+cosx*siny  -sinx*cosz -sinx*siny*sinz+cosx*cosy
         ******************************************************************************/
    case XZY:
    case RYP:
        if(fabs(rotateMatrix(0,0))<DBL_ZERO && fabs(rotateMatrix(0,2))<DBL_ZERO)
        {
            stRotate.dX = DPI;
            stRotate.dY = atan2(-rotateMatrix(0,2),-rotateMatrix(2,2));
            if(rotateMatrix(0,1) < 0.0) stRotate.dZ = -PI_2;
            else                        stRotate.dZ = PI_2;
        }
        else
        {
            stRotate.dX = atan2(-rotateMatrix(2,1),rotateMatrix(1,1));
            stRotate.dY = atan2(-rotateMatrix(0,2),rotateMatrix(0,0));
            stRotate.dZ = atan2(rotateMatrix(0,1),sqrt(rotateMatrix(0,0)*rotateMatrix(0,0)+rotateMatrix(0,2)*rotateMatrix(0,2)));
        }
        break;
        /**********************************************************************************
                 +cosy*cosz-sinx*siny*sinz    +cosy*sinz+sinx*siny*cosz  -cosx*siny
      R(y,x,z) =      -cosx*sinz                      +cosx*cosz           +sinx
                 +siny*cosz+sinx*cosy*sinz    +siny*sinz-sinx*cosy*cosz  +cosx*cosy
         **********************************************************************************/
    case YXZ:
    case PRY:
        if(fabs(rotateMatrix(1,0))<DBL_ZERO && fabs(rotateMatrix(1,1))<DBL_ZERO)
        {
            if(rotateMatrix(1,2) < 0.0) stRotate.dX = -PI_2;
            else                        stRotate.dX = PI_2;

            stRotate.dY = DPI;
            stRotate.dZ = atan2(-rotateMatrix(0,1),-rotateMatrix(0,0));
        }
        else
        {
            stRotate.dX = atan2(rotateMatrix(1,2),sqrt(rotateMatrix(1,0)*rotateMatrix(1,0)+rotateMatrix(1,1)*rotateMatrix(1,1)));
            stRotate.dY = atan2(-rotateMatrix(0,2),rotateMatrix(2,2));
            stRotate.dZ = atan2(-rotateMatrix(1,0),rotateMatrix(1,1));
        }
        break;
        /******************************************************************************
                   +cosy*cosz  +cosx*cosy*sinz+sinx*siny +sinx*cosy*sinz-cosx*siny
      R(y,z,x) =      -sinz       +cosx*cosz                   +sinx*cosz
                   +siny*cosz  +cosx*siny*cosz-sinx*cosy +sinx*siny*sinz+cosx*cosy
         ******************************************************************************/
    case YZX:
    case PYR:
        if(fabs(rotateMatrix(0,0))<DBL_ZERO && fabs(rotateMatrix(2,0))<DBL_ZERO)
        {
            stRotate.dX = DPI;
            stRotate.dY = atan2(-rotateMatrix(0,2),-rotateMatrix(2,2));
            if(rotateMatrix(1,0) < 0.0) stRotate.dZ = PI_2;
            else                        stRotate.dZ = -PI_2;
        }
        else
        {
            stRotate.dX = atan2(rotateMatrix(1,2),rotateMatrix(1,1));
            stRotate.dY = atan2(rotateMatrix(2,0),rotateMatrix(0,0));
            stRotate.dZ = atan2(-rotateMatrix(1,0),sqrt(rotateMatrix(0,0)*rotateMatrix(0,0)+rotateMatrix(2,0)*rotateMatrix(2,0)));
        }
        break;
        /********************************************************************************
                 +cosy*cosz+sinx*siny*sinz   +cosx*sinz  -siny*cosz+sinx*cosy*sinz
      R(z,x,y) = -cosy*sinz+cosx*siny*cosz   +cosx*cosz  +siny*sinz+sinx*cosy*cosz
                       +cosx*siny               -sinx            +cosx*cosy
         ********************************************************************************/
    case ZXY:
    case YRP:
        if(fabs(rotateMatrix(2,0))<DBL_ZERO && fabs(rotateMatrix(2,2))<DBL_ZERO)
        {
            if(rotateMatrix(2,1) < 0.0) stRotate.dX = PI_2;
            else                        stRotate.dX = -PI_2;

            stRotate.dY = DPI;
            stRotate.dZ = atan2(rotateMatrix(1,0),-rotateMatrix(0,0));
        }
        else
        {
            stRotate.dX = atan2(-rotateMatrix(2,1),sqrt(rotateMatrix(0,1)*rotateMatrix(0,1)+rotateMatrix(1,1)*rotateMatrix(1,1)));
            stRotate.dY = atan2(rotateMatrix(2,0),rotateMatrix(2,2));
            stRotate.dZ = atan2(rotateMatrix(0,1),rotateMatrix(1,1));
        }
        break;
        /********************************************************************************
                  +cosy*cosz   +cosx*sinz+sinx*siny*cosz +sinx*sinz-cosx*siny*cosz
      R(z,y,x) =  -cosy*sinz   +cosx*cosz-sinx*siny*sinz +sinx*cosz+cosx*siny*sinz
                     +siny      -sinx*cosy         +cosx*cosy
         ********************************************************************************/
    case ZYX:
    case YPR:
        if(fabs(rotateMatrix(0,0))<DBL_ZERO && fabs(rotateMatrix(1,0))<DBL_ZERO)
        {
            stRotate.dX = DPI;

            if(rotateMatrix(0,2) < 0.0) stRotate.dY = -PI_2;
            else                        stRotate.dY = PI_2;

            stRotate.dZ = atan2(-rotateMatrix(0,1),-rotateMatrix(1,1));
        }
        else
        {
            stRotate.dX = atan2(-rotateMatrix(2,1),rotateMatrix(2,2));
            stRotate.dY = atan2(rotateMatrix(2,0),sqrt(rotateMatrix(2,1)*rotateMatrix(2,1)+rotateMatrix(2,2)*rotateMatrix(2,2)));
            stRotate.dZ = atan2(-rotateMatrix(1,0),rotateMatrix(0,0));
        }
        break;
    default:
        return;
    }
}

/// 根据指定的旋转角度计算旋转矩阵
CMatrix CYPRAngle::CreateMatrix(double dRoll, double dPitch, double dYaw, YPRROTATE type)
{
    double dTemp[3][3];
    iauIr(dTemp);
    switch(type)
    {
    case XYZ:
    case RPY:
        iauRz(dYaw,dTemp);
        iauRy(dPitch,dTemp);
        iauRx(dRoll,dTemp);
        break;
    case XZY:
    case RYP:
        iauRy(dPitch,dTemp);
        iauRz(dYaw,dTemp);
        iauRx(dRoll,dTemp);
        break;
    case YXZ:
    case PRY:
        iauRz(dYaw,dTemp);
        iauRx(dRoll,dTemp);
        iauRy(dPitch,dTemp);
        break;
    case YZX:
    case PYR:
        iauRx(dRoll,dTemp);
        iauRz(dYaw,dTemp);
        iauRy(dPitch,dTemp);
        break;
    case ZXY:
    case YRP:
        iauRy(dPitch,dTemp);
        iauRx(dRoll,dTemp);
        iauRz(dYaw,dTemp);
        break;
    case ZYX:
    case YPR:
        iauRx(dRoll,dTemp);
        iauRy(dPitch,dTemp);
        iauRz(dYaw,dTemp);
        break;
    default:
        break;
    }
    return(CMatrix(dTemp,3,3));
}

/// 根据指定的旋转顺序计算 各角度
YPR_Rotate CYPRAngle::GetRotate(YPRROTATE type) const
{
    if(type == m_emOriginalType)
    {
        return(m_stRotate);
    }
    else
    {
        YPR_Rotate tmp;
        CalTransform(m_matRote,type,tmp);
        return(tmp);
    }
}
