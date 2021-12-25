#include <Math/MathCommonAlgorithm.h>
#include <Math/VecMat.h>
#include <Satellite/CoorSys.h>
#include <Satellite/TimeSys.h>
#include <Satellite/IRESInfo.h>
#include <Satellite/Kepler.h>
#include "sofa.h"
#include "sofam.h"
#include "CommonAlgorithm.h"

using namespace Math;
using namespace Aerospace;
using namespace Satellite;

/// 将TEME的矢量旋转到J2000坐标系下
bool CCoorSys::TEME2J2000(const double &dMJD, const CVector &vTeme, CVector &vJ2000)
{
    /// 判断两个矩阵元素是否相同
    if(!CVecMat::IsValid(vTeme,vJ2000))
    {
        return(false);
    }

    /// 进行矩阵运算
    if(!CVecMat::CalReault(TEME2J2000(dMJD),vTeme,vJ2000))
    {
        return(false);
    }

    return(true);
}

/// 将J2000的矢量旋转到TEME坐标系下
bool CCoorSys::J20002TEME(const double &dMJD, const CVector &vJ2000, CVector &vTeme)
{
    /// 判断两个矩阵元素是否相同
    if(!CVecMat::IsValid(vJ2000,vTeme))
    {
        return(false);
    }

    /// 进行矩阵运算
    if(!CVecMat::CalReault(J20002TEME(dMJD),vJ2000,vTeme))
    {
        return(false);
    }

    return(true);
}

/// TEME到J2000的旋转矩阵
CMatrix CCoorSys::TEME2J2000(const double &dMJD)
{
    /// 正交矩阵的转置与矩阵的逆相同
    return(CVecMat::Transp(J20002TEME(dMJD)));
}

CMatrix CCoorSys::J20002GCRF()
{
    static CMatrix tmpMatrix(3,3);
    static bool bIsInit(false);
    if(!bIsInit)
    {
        tmpMatrix(0,0) = tmpMatrix(1,1) = tmpMatrix(2,2) = 1.0;

        tmpMatrix(0,1) = tmpMatrix(1,0) = -0.000274e-8;
        tmpMatrix(0,2) = tmpMatrix(2,0) = -9.790527e-8;
        tmpMatrix(1,2) = tmpMatrix(2,1) = -1.147780e-8;

        tmpMatrix(0,1) *= -1;
        tmpMatrix(2,0) *= -1;
        tmpMatrix(2,1) *= -1;
        bIsInit = true;
    }
    return(tmpMatrix);
}

CMatrix CCoorSys::GCRF2J2000()
{
    static CMatrix tmpMatrix(3,3);
    static bool bIsInit(false);
    if(!bIsInit)
    {
        tmpMatrix = CVecMat::Transp(J20002GCRF());
        bIsInit = true;
    }

    return(tmpMatrix);
}

bool CCoorSys::TEME2ECF(const double &dMJD, const CVector &vTeme, CVector &vEcf)
{
    /// 判断两个矩阵元素是否相同
    if(!CVecMat::IsValid(vTeme,vEcf))
    {
        return(false);
    }

    /// 进行矩阵运算
    if(!CVecMat::CalReault(TEME2ECF(dMJD),vTeme,vEcf))
    {
        return(false);
    }

    return(true);
}

bool CCoorSys::ECF2TEME(const double &dMJD, const CVector &vEcf, CVector &vTeme)
{
    /// 判断两个矩阵元素是否相同
    if(!CVecMat::IsValid(vEcf,vTeme))
    {
        return(false);
    }

    /// 进行矩阵运算
    if(!CVecMat::CalReault(ECF2TEME(dMJD),vEcf,vTeme))
    {
        return(false);
    }

    return(true);
}

/// TEME2ECF 的旋转矩阵
CMatrix CCoorSys::TEME2ECF(const double &dMJD)
{
    if(dMJD < 1)
    {
        return (CMatrix::IDENTITY_MAT);
    }

    BulletinB tmpB;
    /// 判断是否读入IERS公布的文件
    if(CIRESInfo::GetInstance()->IsInit())
    {
        CIRESInfo::GetInstance()->GetBulletinB(dMJD,tmpB);
    }

    double dTemp[3][3];
    iauIr(dTemp);
    iauRz(iauGmst82(DJM0,dMJD+tmpB.dUT1_UTC/DAYSEC),dTemp);

    PoleMatrix(tmpB.dPX*DAS2R,tmpB.dPY*DAS2R,dTemp);
    return(CMatrix(dTemp[0],3,3));
}

CMatrix CCoorSys::ECF2TEME(const double &dMJD)
{
    return(CVecMat::Transp(TEME2ECF(dMJD)));
}

/// 将J2000的矢量旋转到ECF坐标系下
bool CCoorSys::J20002ECF(const double &dMJD, const CVector &vJ2000, CVector &vEcf)
{
    /// 判断两个矩阵元素是否相同
    if(!CVecMat::IsValid(vJ2000,vEcf))
    {
        return(false);
    }

    /// 进行矩阵运算
    if(!CVecMat::CalReault(J20002ECF(dMJD),vJ2000,vEcf))
    {
        return(false);
    }

    return(true);
}

/// 将ECF的矢量旋转到J2000坐标系下
bool CCoorSys::ECF2J2000(const double &dMJD, const CVector &vEcf, CVector &vJ2000)
{
    /// 判断两个矩阵元素是否相同
    if(!CVecMat::IsValid(vEcf,vJ2000))
    {
        return(false);
    }

    /// 进行矩阵运算
    if(!CVecMat::CalReault(ECF2J2000(dMJD),vEcf,vJ2000))
    {
        return(false);
    }

    return(true);
}

/// J2000到ECF的旋转矩阵
CMatrix CCoorSys::J20002ECF(const double &dMJD)
{
    if(dMJD < 1)
    {
        return (CMatrix::IDENTITY_MAT);
    }

    BulletinB tmpB;
    /// 判断是否读入IERS公布的文件
    if(CIRESInfo::GetInstance()->IsInit())
    {
        CIRESInfo::GetInstance()->GetBulletinB(dMJD,tmpB);
    }

    /// 进行时间转换
    CTimeSys timeSys(dMJD);
    double dTT = timeSys.GetTT();

    double dTemp[3][3];
    iauIr(dTemp);

    PrecMatrix(dTT,dTemp);

    double dpsi,eps,deps,dom;
    NutMatrix(dTT,dpsi,eps,deps,dom,dTemp);
    GHAMatrix(dMJD+tmpB.dUT1_UTC/DAYSEC,dpsi,deps,eps,dom,dTemp);

    PoleMatrix(tmpB.dPX*DAS2R,tmpB.dPY*DAS2R,dTemp);

    return(CMatrix(dTemp[0],3,3));
}

/// ECF到J2000的旋转矩阵
CMatrix CCoorSys::ECF2J2000(const double &dMJD)
{
    return(CVecMat::Transp(J20002ECF(dMJD)));
}

CMatrix CCoorSys::J20002MOD(const double &dMJD)
{
    if(dMJD < 1)
    {
        return (CMatrix::IDENTITY_MAT);
    }

    /// 进行时间转换
    CTimeSys timeSys(dMJD);
    double dTT = timeSys.GetTT();

    double dTemp[3][3];
    iauIr(dTemp);

    PrecMatrix(dTT,dTemp);
    return(CMatrix(dTemp,3,3));
}

CMatrix CCoorSys::MOD2J2000(const double &dMJD)
{
    return(CVecMat::Transp(J20002MOD(dMJD)));
}

CMatrix CCoorSys::J20002TOD(const double &dMJD)
{
    if(dMJD < 1)
    {
        return (CMatrix::IDENTITY_MAT);
    }

    /// 进行时间转换
    CTimeSys timeSys(dMJD);
    double dTT = timeSys.GetTT();

    double dTemp[3][3];
    iauIr(dTemp);

    PrecMatrix(dTT,dTemp);
    double dpsi,eps,deps,dom;
    NutMatrix(dTT,dpsi,eps,deps,dom,dTemp);
    return(CMatrix(dTemp,3,3));
}

CMatrix CCoorSys::TOD2J2000(const double &dMJD)
{
    return(CVecMat::Transp(J20002TOD(dMJD)));
}

CMatrix CCoorSys::J20002TEME(const double &dMJD)
{
    if(dMJD < 1)
    {
        return (CMatrix::IDENTITY_MAT);
    }

    /// 进行时间转换
    CTimeSys timeSys(dMJD);
    double dTT = timeSys.GetTT();

    double dTemp[3][3];
    iauIr(dTemp);

    PrecMatrix(dTT,dTemp);
    double dpsi,eps,deps,dom;
    NutMatrix(dTT,dpsi,eps,deps,dom,dTemp);
    iauRz(EquationOfEquinoxes(dpsi,deps,eps,dom),dTemp);

    return(CMatrix(dTemp,3,3));
}

/// 将VVLH的矢量旋转到ECI坐标系下
bool CCoorSys::VVLH2J2000(const CVector &vKep, const CVector &vVVLH, CVector &vECI)
{
    /// 判断两个矩阵元素是否相同
    if(!CVecMat::IsValid(vVVLH,vECI))
    {
        return(false);
    }

    /// 进行矩阵运算
    if(!CVecMat::CalReault(VVLH2J2000(vKep),vVVLH,vECI))
    {
        return(false);
    }

    return(true);
}

/// 将ECI的矢量旋转到VVLH坐标系下
bool CCoorSys::J20002VVLH(const CVector &vKep, const CVector &vECI, CVector &vVVLH)
{
    /// 判断两个矩阵元素是否相同
    if(!CVecMat::IsValid(vECI,vVVLH))
    {
        return(false);
    }

    /// 进行矩阵运算
    if(!CVecMat::CalReault(J20002VVLH(vKep),vECI,vVVLH))
    {
        return(false);
    }

    return(true);
}

/// 轨道坐标系到ECI的旋转矩阵
CMatrix CCoorSys::VVLH2J2000(const CVector &vKep)
{
    return(CVecMat::Transp(J20002VVLH(vKep)));
}

/// ECI到轨道坐标系的旋转矩阵
CMatrix CCoorSys::J20002VVLH(const CVector &vKep)
{
    /// 判断数据是否有效
    if(vKep.Size() < 6)
    {
        return(CMatrix::IDENTITY_MAT);
    }

    double dTemp[3][3];
    iauIr(dTemp);
    iauRz(vKep(3),dTemp);
    iauRx(vKep(2),dTemp);
    iauRz(Modulo(CKepler::TrueAnom(vKep(1),vKep(5))+vKep(4),D2PI)+DPI/2.0,dTemp);
    iauRz(-DPI/2.0,dTemp);
    return(CMatrix(dTemp,3,3));
}
