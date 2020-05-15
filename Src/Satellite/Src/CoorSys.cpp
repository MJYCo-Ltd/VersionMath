#include "CoorSys.h"
#include "VecMat.h"
#include "TimeSys.h"
#include "IRESInfo.h"
#include "sofa.h"
#include "MathCommonAlgorithm.h"
#include "CommonAlgorithm.h"
#include "Kepler.h"

using namespace Math;
using namespace Aerospace;
using namespace Satellite;

CCoorSys::CCoorSys()
{
}

/// 将ECI的矢量旋转到ECF坐标系下
bool CCoorSys::ECI2ECF(const double &dMJD, const CVector &vEci, CVector &vEcf)
{
    /// 判断两个矩阵元素是否相同
    if(!IsValid(vEci,vEcf))
    {
        return(false);
    }

    /// 计算ECI到ECF转换矩阵
    CMatrix tmpMat = ECI2ECF(dMJD);

    /// 进行矩阵运算
    if(!CalReault(tmpMat,vEci,vEcf))
    {
        return(false);
    }

    return(true);
}

/// 将ECF的矢量旋转到ECI坐标系下
bool CCoorSys::ECF2ECI(const double &dMJD, const CVector &vEcf, CVector &vEci)
{
    /// 判断两个矩阵元素是否相同
    if(!IsValid(vEcf,vEci))
    {
        return(false);
    }

    /// 计算ECI到ECF转换矩阵
    CMatrix tmpMat = ECF2ECI(dMJD);

    /// 进行矩阵运算
    if(!CalReault(tmpMat,vEcf,vEci))
    {
        return(false);
    }

    return(true);
}

/// ECI到ECF的旋转矩阵
CMatrix CCoorSys::ECI2ECF(const double &dMJD)
{
    CMatrix tmpMat;
    if(dMJD <= 0)
    {
        return (tmpMat);
    }

    CIRESInfo* pIERS = CIRESInfo::GetInstance();
    BulletinB tmpB;

    /// 判断是否读入IERS公布的文件
    if(pIERS && pIERS->IsInit())
    {
        tmpB = pIERS->GetBulletinB(dMJD);
    }
    else
    {
        return(tmpMat);
    }

    /// 进行时间转换
    CTimeSys timeSys(dMJD);

    static double dTemp[3][3];

    /// 进行运算
    iauC2t06a(DJM0,timeSys.GetTT(),DJM0,timeSys.GetUT1(),tmpB.dPX*DAS2R,tmpB.dPY*DAS2R,dTemp);
    return(CMatrix(dTemp[0],3,3));
}

/// ECF到ECI的旋转矩阵
CMatrix CCoorSys::ECF2ECI(const double &dMJD)
{
    return (CVecMat::Transp(ECI2ECF(dMJD)));
}

/// 将TEME的矢量旋转到ECI坐标系下
bool CCoorSys::TEME2ECI(const double &dMJD, const CVector &vTeme, CVector &vEci)
{
    /// 判断两个矩阵元素是否相同
    if(!IsValid(vTeme,vEci))
    {
        return(false);
    }

    CMatrix tmpMat = TEME2ECI(dMJD);

    /// 进行矩阵运算
    if(!CalReault(tmpMat,vTeme,vEci))
    {
        return(false);
    }

    return(true);
}

/// 将ECI的矢量旋转到TEME坐标系下
bool CCoorSys::ECI2TEME(const double &dMJD, const CVector &vEci, CVector &vTeme)
{
    /// 判断两个矩阵元素是否相同
    if(!IsValid(vEci,vTeme))
    {
        return(false);
    }

    CMatrix tmpMat = ECI2TEME(dMJD);

    /// 进行矩阵运算
    if(!CalReault(tmpMat,vEci,vTeme))
    {
        return(false);
    }

    return(true);
}

/// TEME到ECI的旋转矩阵
CMatrix CCoorSys::TEME2ECI(const double &dMJD)
{
    /// 正交矩阵的转置与矩阵的逆相同
    return(CVecMat::Transp(ECI2TEME(dMJD)));
}

/// ECI到TEME的旋转矩阵
CMatrix CCoorSys::ECI2TEME(const double &dMJD)
{
    CMatrix tmpMat;
    if(dMJD < 0.)
    {
        return(tmpMat);
    }

    /// 进行时间转换
    CTimeSys timeSys(dMJD);

    /// 获得TT时间
    double dTT = timeSys.GetTT();
    /// 黄道倾斜度（弧度）
    double dObecm = iauObl06(DJM0, dTT);

    double dLoec;
    double dObec;
    iauNut06a(DJM0, dTT, &dLoec, &dObec);

    double dSG = dLoec * cos(dObec + dObecm);

    CMatrix temedMat = CVecMat::R_z(dSG);

    /// 岁差章动矩阵
    double Pn_Matrix[3][3];
    iauPnm06a(DJM0, dTT, Pn_Matrix);
    CMatrix pnMatix(Pn_Matrix,3,3);
    /// 岁差章动、TEME矩阵相乘
    return(pnMatix * temedMat);
}

bool CCoorSys::TEME2ECF(const double &dMJD, const CVector &vTeme, CVector &vEcf)
{
    /// 判断两个矩阵元素是否相同
    if(!IsValid(vTeme,vEcf))
    {
        return(false);
    }

    CMatrix tmpMat = TEME2ECF(dMJD);

    /// 进行矩阵运算
    if(!CalReault(tmpMat,vTeme,vEcf))
    {
        return(false);
    }

    return(true);
}

bool CCoorSys::ECF2TEME(const double &dMJD, const CVector &vEcf, CVector &vTeme)
{
    /// 判断两个矩阵元素是否相同
    if(!IsValid(vEcf,vTeme))
    {
        return(false);
    }

    CMatrix tmpMat = ECF2TEME(dMJD);

    /// 进行矩阵运算
    if(!CalReault(tmpMat,vEcf,vTeme))
    {
        return(false);
    }

    return(true);
}

/// TEME2ECF 的旋转矩阵
CMatrix CCoorSys::TEME2ECF(const double &dMJD)
{
    CMatrix tmpMat;
    if(dMJD <= 0)
    {
        return (tmpMat);
    }

    /// 时间转换
    CTimeSys timeSys(dMJD);
    double dUt = timeSys.GetUT1();

    CIRESInfo* pIERS = CIRESInfo::GetInstance();
    BulletinB tmpB;

    /// 判断是否读入IERS公布的文件
    if(pIERS && pIERS->IsInit())
    {
        tmpB = pIERS->GetBulletinB(dMJD);
    }
    else
    {
        return(tmpMat);
    }

    return(PoleMatrix(tmpB.dPX*DAS2R,tmpB.dPY*DAS2R)*GHAMatrix(dUt));
}

CMatrix CCoorSys::ECF2TEME(const double &dMJD)
{
    return(CVecMat::Transp(TEME2ECF(dMJD)));
}

/// 将J2000的矢量旋转到ECF坐标系下
bool CCoorSys::J20002ECF(const double &dMJD, const CVector &vJ2000, CVector &vEcf)
{
    /// 判断两个矩阵元素是否相同
    if(!IsValid(vJ2000,vEcf))
    {
        return(false);
    }

    CMatrix tmpMat = J20002ECF(dMJD);

    /// 进行矩阵运算
    if(!CalReault(tmpMat,vJ2000,vEcf))
    {
        return(false);
    }

    return(true);
}

/// 将ECF的矢量旋转到J2000坐标系下
bool CCoorSys::ECF2J2000(const double &dMJD, const CVector &vEcf, CVector &vJ2000)
{
    /// 判断两个矩阵元素是否相同
    if(!IsValid(vEcf,vJ2000))
    {
        return(false);
    }

    CMatrix tmpMat = ECF2J2000(dMJD);

    /// 进行矩阵运算
    if(!CalReault(tmpMat,vEcf,vJ2000))
    {
        return(false);
    }

    return(true);
}

/// J2000到ECF的旋转矩阵
CMatrix CCoorSys::J20002ECF(const double &dMJD)
{
    CMatrix tmpMat;
    if(dMJD <= 0)
    {
        return (tmpMat);
    }

    CIRESInfo* pIERS = CIRESInfo::GetInstance();
    BulletinB tmpB;

    /// 判断是否读入IERS公布的文件
    if(pIERS && pIERS->IsInit())
    {
        tmpB = pIERS->GetBulletinB(dMJD);
    }
    else
    {
        return(tmpMat);
    }

    /// 进行时间转换
    CTimeSys timeSys(dMJD);
    double dTT = timeSys.GetTT();
    double dUt = timeSys.GetUT1();

    return(PoleMatrix(tmpB.dPX*DAS2R,tmpB.dPY*DAS2R)*GHAMatrix(dUt)*NutMatrix(dTT)*PrecMatrix(DJM00,dTT));
}

/// ECF到J2000的旋转矩阵
CMatrix CCoorSys::ECF2J2000(const double &dMJD)
{
    return(CVecMat::Transp(J20002ECF(dMJD)));
}

/// 将VVLH的矢量旋转到ECI坐标系下
bool CCoorSys::VVLH2ECI(const CVector &vKep, const CVector &vVVLH, CVector &vECI)
{
    /// 判断两个矩阵元素是否相同
    if(!IsValid(vVVLH,vECI))
    {
        return(false);
    }

    CMatrix tmpMat = VVLH2ECI(vKep);

    /// 进行矩阵运算
    if(!CalReault(tmpMat,vVVLH,vECI))
    {
        return(false);
    }

    return(true);
}

/// 将ECI的矢量旋转到VVLH坐标系下
bool CCoorSys::ECI2VVLH(const CVector &vKep, const CVector &vECI, CVector &vVVLH)
{
    /// 判断两个矩阵元素是否相同
    if(!IsValid(vECI,vVVLH))
    {
        return(false);
    }

    CMatrix tmpMat = ECI2VVLH(vKep);

    /// 进行矩阵运算
    if(!CalReault(tmpMat,vECI,vVVLH))
    {
        return(false);
    }

    return(true);
}

/// 轨道坐标系到ECI的旋转矩阵
CMatrix CCoorSys::VVLH2ECI(const CVector &vKep)
{
    return(CVecMat::Transp(ECI2VVLH(vKep)));
}

/// ECI到轨道坐标系的旋转矩阵
CMatrix CCoorSys::ECI2VVLH(const CVector &vKep)
{
    CMatrix tmp;

    /// 判断数据是否有效
    if(vKep.Size() < 6)
    {
        return(tmp);
    }

    double dTmp = Modulo(CKepler::TrueAnom(vKep(1),vKep(5))+vKep(4),D2PI)+DPI/2.0;

    return(CVecMat::R_x(-DPI/2.0)*CVecMat::R_z(dTmp)*
           CVecMat::R_x(vKep(2))*CVecMat::R_z(vKep(3)));
}

void CCoorSys::LBH2XYZ(double dL, double dB, double dH, double &dX, double &dY, double &dZ)
{
    static double dTemp[3];

    iauGd2gc(WGS84,dL,dB,dH,dTemp);

    dX = dTemp[0];
    dY = dTemp[1];
    dZ = dTemp[2];
}

void CCoorSys::XYZ2LBH(double dX, double dY, double dZ, double &dL, double &dB, double &dH)
{
    static double dTemp[3];
    dTemp[0] = dX;
    dTemp[1] = dY;
    dTemp[2] = dZ;

    iauGc2gd(WGS84,dTemp,&dL,&dB,&dH);
}

bool CCoorSys::LBH2XYZ(const CVector &vGeo, CVector &v3D)
{

    int nGeo = vGeo.Size();

    /// 判断两个矩阵元素是否相同
    if(!IsValid(vGeo,v3D))
    {
        return(false);
    }

    for(int i=0; i<nGeo; i+=3)
    {
        LBH2XYZ(vGeo(i),vGeo(i+1),vGeo(i+2),v3D(i),v3D(i+1),v3D(i+2));
    }

    return(true);
}

bool CCoorSys::XYZ2LBH(const CVector &v3D, CVector &vGeo)
{
    int n3D = v3D.Size();

    /// 判断两个矩阵元素是否相同
    if(!IsValid(v3D,vGeo))
    {
        return(false);
    }

    for(int i=0; i<n3D; i+=3)
    {
        XYZ2LBH(v3D(i),v3D(i+1),v3D(i+2),vGeo(i),vGeo(i+1),vGeo(i+2));
    }
    return(true);
}

bool CCoorSys::GLOBAL2LOCAL(double dL, double dB, const CVector &vGlobal, CVector &vLocal)
{
    /// 判断两个向量元素是否相同
    if(!IsValid(vGlobal,vLocal))
    {
        return(false);
    }

    CMatrix tmpMat = GLOBAL2LOCAL(dL,dB);

    /// 进行矩阵运算
    if(!CalReault(tmpMat,vGlobal,vLocal))
    {
        return(false);
    }

    return(true);
}

bool CCoorSys::LOCAL2GLOBAL(double dL, double dB, const CVector &vLocal, CVector &vGlobal)
{
    /// 判断两个矩阵元素是否相同
    if(!IsValid(vLocal,vGlobal))
    {
        return(false);
    }

    CMatrix tmpMat = LOCAL2GLOBAL(dL,dB);
    /// 查看矩阵是否有效
    if(!tmpMat)
    {
        return(false);
    }

    /// 进行矩阵运算
    if(!CalReault(tmpMat,vLocal,vGlobal))
    {
        return(false);
    }

    return(true);
}

CMatrix CCoorSys::GLOBAL2LOCAL(double dL, double dB)
{
    double  Aux;

    // Transformation to Zenith-East-North System
    CMatrix M = CVecMat::R_y(-dB)*CVecMat::R_z(dL);

    // Cyclic shift of rows 0,1,2 to 1,2,0 to obtain East-North-Zenith system
    for (int j=0; j<3; ++j)
    {
        Aux=M(0,j);
        M(0,j)=M(1,j);
        M(1,j)=M(2,j);
        M(2,j)= Aux;
    }

    return(M);
}

CMatrix CCoorSys::LOCAL2GLOBAL(double dL, double dB)
{
    return (CVecMat::Transp(GLOBAL2LOCAL(dL,dB)));
}
