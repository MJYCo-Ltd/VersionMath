#include <cmath>
#include "GisMath_Common.h"

/// 全局坐标到局部坐标的转换
bool GisMath::GLOBAL2LOCAL(double dL, double dB, const Math::CVector &vGlobal, Math::CVector &vLocal)
{
    /// 判断两个向量元素是否相同
    if(!Math::CVecMat::IsValid(vGlobal,vLocal))
    {
        return(false);
    }

    Math::CMatrix tmpMat = GLOBAL2LOCAL(dL,dB);

    /// 进行矩阵运算
    if(!Math::CVecMat::CalReault(tmpMat,vGlobal,vLocal))
    {
        return(false);
    }

    return(true);
}

/// 局部坐标转到全局坐标
bool GisMath::LOCAL2GLOBAL(double dL, double dB, const Math::CVector &vLocal, Math::CVector &vGlobal)
{
    /// 判断两个矩阵元素是否相同
    if(!Math::CVecMat::IsValid(vLocal,vGlobal))
    {
        return(false);
    }

    Math::CMatrix tmpMat = LOCAL2GLOBAL(dL,dB);
    /// 查看矩阵是否有效
    if(!tmpMat)
    {
        return(false);
    }

    /// 进行矩阵运算
    if(!Math::CVecMat::CalReault(tmpMat,vLocal,vGlobal))
    {
        return(false);
    }

    return(true);
}

/// 将全局坐标转成局部坐标
Math::CMatrix GisMath::GLOBAL2LOCAL(double dL, double dB)
{
    double  Aux;

    // Transformation to Zenith-East-North System
    Math::CMatrix M = Math::CVecMat::R_y(-dB)*Math::CVecMat::R_z(dL);

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

/// 局部坐标转成全局坐标
Math::CMatrix GisMath::LOCAL2GLOBAL(double dL, double dB)
{
    return (Math::CVecMat::Transp(GLOBAL2LOCAL(dL,dB)));
}

/// 根据距离 站心 的俯仰、方位角、距离 解算在世界坐标系下的位置
void Local2Globel(double dLon, double dLat, double dAzim, double dElev, double dDist,
                  double &dX, double& dY, double& dZ)
{
    /// 用于存放临时数据
    Math::CVector vTemp(3);

    /// 计算在站心坐标系下的位置
    /// 因球坐标以逆时针为正方向，而方位角以顺时针为正方向
    /// 故此处对方位角取反；同时球坐标以正东为起始位置，此处
    /// 需要增加90度
    vTemp = Math::CVecMat::VecPolar(DPI*0.5-dAzim,dElev,dDist);

    /// 在世界坐标下 所求位置 相对于给定位置的位置
    vTemp = vTemp * GisMath::GLOBAL2LOCAL(dLon,dLat);

    dX = vTemp(0);
    dY = vTemp(1);
    dZ = vTemp(2);
}

/// 根据世界坐标系下的位置解算 在站心坐标系下 的俯仰、方位角、距离
void Globel2Local(double dLon, double dLat, double dX, double dY, double dZ,
                  double &dAzim, double &dElev, double &dDist)
{
    /// 计算 世界坐标系 到 局部天东北坐标系 旋转矩阵
    Math::CMatrix M = GisMath::GLOBAL2LOCAL(dLon,dLat);

    /// 用于存放临时数据
    Math::CVector vTemp(3);
    vTemp(0) = dX;
    vTemp(1) = dY;
    vTemp(2) = dZ;

    /// 将世界坐标转换到站心坐标
    vTemp = M * vTemp;

    /// 计算在站心坐标系的俯仰，方位
    Math::CVecMat::AzEl(vTemp,dAzim,dElev);
    dDist = vTemp.Length();
}

/// 根据给定的地理坐标，俯仰、方位、距离，计算所求点的地理坐标
int GisMath::GeoCalEndGeo(double dLon, double dLat, double dHeight, double dAzim, double dElev, double dDist, double &dLon2, double &dLat2, double &dHeight2)
{
    /// 用于存放世界坐标
    double dTemp[3];

    CheckInit();

    /// 根据地理坐标计算世界坐标
    iauGd2gc(EM_TYPE,dLon,dLat,dHeight,dTemp);

    /// 计算方位、俯仰、距离对应的数据
    Local2Globel(dLon,dLat,dAzim,dElev,dDist,dLon2,dLat2,dHeight2);

    /// 计算所求点的世界坐标
    dTemp[0] += dLon2;
    dTemp[1] += dLat2;
    dTemp[2] += dHeight2;

    /// 将世界坐标转换成 地理坐标
    iauGc2gd(EM_TYPE,dTemp,&dLon2,&dLat2,&dHeight2);

    return(0);
}

/// 根据给定经纬度，俯仰、方位、翻滚、偏移地址，计算所求点的地理坐标
int GisMath::GeoCalEndGeo(double dLon, double dLat, double dHeight,
                          double dAzim, double dElev, double dRoll,
                          double dX, double dY, double dZ,
                          double &dLon2, double &dLat2, double &dHeight2)
{
    /// 用于存放世界坐标
    double dTemp[3];

    CheckInit();

    /// 根据地理坐标计算世界坐标
    iauGd2gc(EM_TYPE,dLon,dLat,dHeight,dTemp);


    /// 构建偏移向量
    Math::CVector offset(dX,dY,dZ);

    /// 方位旋转
    Math::CMatrix mAzim = Math::CVecMat::R_z(dAzim);

    /// 俯仰旋转
    Math::CMatrix mElev = Math::CVecMat::R_x(-dElev);

    /// 横滚旋转
    Math::CMatrix mRoll = Math::CVecMat::R_y(-dRoll);

    Math::CMatrix tmpMatrix = mAzim*mElev*mRoll;

    offset = tmpMatrix*offset;

    /// 求 东北天 到 地心 的 旋转矩阵
    Math::CMatrix M = LOCAL2GLOBAL(dLon,dLat);

    /// 在世界坐标下 所求位置 相对于给定位置的位置
    offset = M * offset;

    dTemp[0] += offset(0);
    dTemp[1] += offset(1);
    dTemp[2] += offset(2);

    /// 将世界坐标转换成 地理坐标
    iauGc2gd(EM_TYPE,dTemp,&dLon2,&dLat2,&dHeight2);

    return(0);
}

/// 根据给定的世界坐标，俯仰、方位、距离，计算所求点的世界坐标
int GisMath::XYZCalEndXYZ(double dX, double dY, double dZ, double dAzim, double dElev, double dDist, double &dX2, double &dY2, double &dZ2)
{
    /// 用于存放世界坐标
    double dTemp[3],dLon,dLat,dHeight;

    CheckInit();

    /// 赋值
    dTemp[0] = dX;
    dTemp[1] = dY;
    dTemp[2] = dZ;

    /// 计算经纬高
    iauGc2gd(EM_TYPE,dTemp,&dLon,&dLat,&dHeight);

    /// 计算方位、俯仰、距离对应的数据
    Local2Globel(dLon,dLat,dAzim,dElev,dDist,dX,dY,dZ);

    /// 计算所求点的世界坐标
    dX2 = dTemp[0] + dX;
    dY2 = dTemp[1] + dY;
    dZ2 = dTemp[2] + dZ;

    return(0);
}

/// 根据两个地理坐标计算第二个相对于第一个的方位角、俯仰角
int GisMath::CalAzElGeo(double dLon1, double dLat1, double dHeight1, double dLon2, double dLat2, double dHeight2, double &dAzim, double &dElev, double &dDist)
{
    double dTemp1[3],dTemp2[3],dx,dy,dz;

    CheckInit();

    /// 根据地理坐标计算世界坐标
    iauGd2gc(EM_TYPE,dLon1,dLat1,dHeight1,dTemp1);
    iauGd2gc(EM_TYPE,dLon2,dLat2,dHeight2,dTemp2);

    /// 计算在 世界坐标系下 的相对坐标
    dx = dTemp2[0] - dTemp1[0];
    dy = dTemp2[1] - dTemp1[1];
    dz = dTemp2[2] - dTemp1[2];

    /// 计算 在站心坐标系下的 俯仰、方位
    Globel2Local(dLon1,dLat1,dx,dy,dz,dAzim,dElev,dDist);

    return(0);
}

int GisMath::CalAzElGeo(double dLon1, double dLat1, double dHeight1,
                        double dLon2, double dLat2, double dHeight2,
                        double dAzim, double dElev, double dRoll,
                        double &dX, double &dY, double &dZ)
{
    double dTemp1[3],dTemp2[3],dx,dy,dz;

    CheckInit();

    /// 根据地理坐标计算世界坐标
    iauGd2gc(EM_TYPE,dLon1,dLat1,dHeight1,dTemp1);
    iauGd2gc(EM_TYPE,dLon2,dLat2,dHeight2,dTemp2);

    /// 计算在 世界坐标系下 的相对坐标
    dx = dTemp2[0] - dTemp1[0];
    dy = dTemp2[1] - dTemp1[1];
    dz = dTemp2[2] - dTemp1[2];

    Math::CVector vTemp(dx,dy,dz);

    /// 计算 世界坐标系 到 局部天东北坐标系 旋转矩阵
    Math::CMatrix M = LOCAL2GLOBAL(dLon1,dLat1);

    /// 将世界坐标转换到站心坐标
    vTemp = M * vTemp;

    /// 方位旋转
    Math::CMatrix mAzim = Math::CVecMat::R_z(dAzim);

    /// 俯仰旋转
    Math::CMatrix mElev = Math::CVecMat::R_x(-dElev);

    /// 横滚旋转
    Math::CMatrix mRoll = Math::CVecMat::R_y(-dRoll);


    vTemp = vTemp * mAzim * mElev * mRoll;

    dX = vTemp(0);
    dY = vTemp(1);
    dZ = vTemp(2);

    return(0);
}

/// 根据两个世界坐标计算第二个相对于第一个的方位角、俯仰角
int GisMath::CalAzElXYZ(double dX1, double dY1, double dZ1, double dX2, double dY2, double dZ2, double &dAzim, double &dElev, double &dDist)
{
    double dTemp[3],dLon,dLat,dHeight,dx,dy;

    CheckInit();

    /// 赋值
    dTemp[0] = dX1;
    dTemp[1] = dY1;
    dTemp[2] = dZ1;

    /// 根据地理坐标计算世界坐标
    iauGc2gd(EM_TYPE,dTemp,&dLon,&dLat,&dHeight);

    /// 计算在 世界坐标系下 的相对坐标
    dx = dX2 - dX1;
    dy = dY2 - dY1;
    dHeight = dZ2 - dZ1;

    /// 计算 在站心坐标系下的 俯仰、方位
    Globel2Local(dLon,dLat,dx,dy,dHeight,dAzim,dElev,dDist);

    return(0);
}

bool GisMath::WGS842GJC02(double& dLon,double& dLat)
{
    double x=dLon - 105.0;
    double y=dLat - 35.0;

    double dTempLat = -100.0 + 2.0 * x + 3.0 * y + 0.2 * y * y + 0.1 * x * y + 0.2 * sqrt(abs(x));
    dTempLat += (20.0 * sin(6.0 * x * DPI) + 20.0 * sin(2.0 * x * DPI)) * 2.0 / 3.0;
    dTempLat += (20.0 * sin(y * DPI) + 40.0 * sin(y / 3.0 * DPI)) * 2.0 / 3.0;
    dTempLat += (160.0 * sin(y / 12.0 * DPI) + 320 * sin(y * DPI / 30.0)) * 2.0 / 3.0;

    double dTempLon = 300.0 + x + 2.0 * y + 0.1 * x * x + 0.1 * x * y + 0.1 * sqrt(abs(x));
    dTempLon += (20.0 * sin(6.0 * x * DPI) + 20.0 * sin(2.0 * x * DPI)) * 2.0 / 3.0;
    dTempLon += (20.0 * sin(x * DPI) + 40.0 * sin(x / 3.0 * DPI)) * 2.0 / 3.0;
    dTempLon += (150.0 * sin(x / 12.0 * DPI) + 300.0 * sin(x / 30.0 * DPI)) * 2.0 / 3.0;

    double radLat = dLat *DD2R;
    double magic = sin(radLat);
    magic = 1 - 0.00669342162296594323 * magic * magic;
    double sqrtMagic = sqrt(magic);
    dTempLat = (dTempLat * 180.0) / ((6378245.0 * (1 - 0.00669342162296594323)) / (magic * sqrtMagic) * DPI);
    dTempLon = (dTempLon * 180.0) / (6378245.0 / sqrtMagic * cos(radLat) * DPI);
    dLat += dTempLat;
    dLon += dTempLon;
    return(true);
}
