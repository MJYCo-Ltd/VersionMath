#include <cmath>
#include "GisMath_Common.h"

/// 经纬度和xyz的转换
bool GisMath::LBH2XYZ(double dL, double dB, double dH, double &dX, double &dY, double &dZ)
{
    double dTemp[3];

    if(0 == iauGd2gc(EM_TYPE,dL,dB,dH,dTemp))
    {
        dX = dTemp[0];
        dY = dTemp[1];
        dZ = dTemp[2];
        return(true);
    }
    else
    {
        return(false);
    }
}

/// 经纬度和xyz的转换
bool GisMath::XYZ2LBH(double dX, double dY, double dZ, double &dL, double &dB, double &dH)
{
    double dTemp[3];
    dTemp[0] = dX;
    dTemp[1] = dY;
    dTemp[2] = dZ;

    return(0 == iauGc2gd(EM_TYPE,dTemp,&dL,&dB,&dH));
}

/// 经纬度和xyz的转换
bool GisMath::LBH2XYZ(const CVector &vGeo, CVector &v3D)
{

    int nGeo = vGeo.Size();

    /// 判断两个矩阵元素是否相同
    if(!CVecMat::IsValid(vGeo,v3D))
    {
        return(false);
    }

    for(int i=0; i<nGeo; i+=3)
    {
        if(!LBH2XYZ(vGeo(i),vGeo(i+1),vGeo(i+2),v3D(i),v3D(i+1),v3D(i+2)))
        {
            return(false);
        }
    }

    return(true);
}

/// 经纬度和xyz的转换
bool GisMath::XYZ2LBH(const CVector &v3D, CVector &vGeo)
{
    int n3D = v3D.Size();

    /// 判断两个矩阵元素是否相同
    if(!CVecMat::IsValid(v3D,vGeo))
    {
        return(false);
    }

    for(int i=0; i<n3D; i+=3)
    {
        if(!XYZ2LBH(v3D(i),v3D(i+1),v3D(i+2),vGeo(i),vGeo(i+1),vGeo(i+2)))
        {
            return(false);
        }
    }
    return(true);
}

/// 贝塞尔大地主题解算 正解
int  GisMath::CalBaiser(double dLon, double dLat, double dAzim, double dDist, double &dLon2, double &dLat2)
{
    CheckInit();

    /// 弧度转成度
    dLon *= DR2D;
    dLat *= DR2D;
    dAzim *= DR2D;

    geod_direct(&PJ_GEOD,dLat,dLon,dAzim,dDist,&dLat2,&dLon2,0);


    /// 度转成弧度
    dLat2 *= DD2R;
    dLon2 *= DD2R;

    return(0);
}

/// 贝塞尔大地主题解算 反解
int GisMath::CalBaiserF(double dLon1, double dLat1, double dLon2, double dLat2, double &dAzim1, double &dAzim2, double &dDist)
{
    CheckInit();

    dLon1 *= DR2D;
    dLat1 *= DR2D;
    dLon2 *= DR2D;
    dLat2 *= DR2D;

    geod_inverse(&PJ_GEOD,dLat1,dLon1,dLat2,dLon2,&dDist,&dAzim1,&dAzim2);

    dAzim1 = iauAnp(dAzim1*DD2R);
    dAzim2 = iauAnp(dAzim2*DD2R);

    return(0);
}

/// 计算两个地理坐标的弧长
double GisMath::CalArcLen(double dLon1, double dLat1, double dLon2, double dLat2, GisMath::LENGTH_TYPE lenTYpe)
{
    double dAzim1,dAzim2,dDist;
    CalBaiserF(dLon1,dLat1,dLon2,dLat2,dAzim1,dAzim2,dDist);

    switch (lenTYpe)
    {
    default:
        return(dDist);
    case KM:
        return(dDist*M2KM);
    case NMI:
        return(dDist*M2NMI);
    }
}

/// 计算多边形围成的面积
double GisMath::CalclutaGeoArea(const CVector &vGeoIn, GisMath::AREA_TYPE type)
{
    CheckInit();
    /// 不足以构成多边形
    int n_xy = vGeoIn.Size()/2;
    if(n_xy < 3)
    {
        return(-1);
    }

    /// 构建临时数据
    double dArea,dLength,*dLon=new double[n_xy](),*dLat=new double[n_xy]();
    for(short i=0; i<n_xy; ++i)
    {
        dLon[i] = vGeoIn(i*2);
        dLat[i] = vGeoIn(i*2+1);
    }

    /// 计算多边形面积
    geod_polygonarea(&PJ_GEOD,dLat,dLon,n_xy,&dArea,&dLength);

    /// 释放空间
    delete[]dLon;
    delete[]dLat;
    switch (type)
    {
    default:
        return(fabs(dArea));
    case KMsqu:
        return(fabs(dArea)*Msqu2KMsqu);
    }
}
