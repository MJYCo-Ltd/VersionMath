#include "GisMath_Common.h"

/// 全局坐标到局部坐标的转换
bool GisMath::GLOBAL2LOCAL(double dL, double dB, const CVector &vGlobal, CVector &vLocal)
{
    /// 判断两个向量元素是否相同
    if(!CVecMat::IsValid(vGlobal,vLocal))
    {
        return(false);
    }

    CMatrix tmpMat = GLOBAL2LOCAL(dL,dB);

    /// 进行矩阵运算
    if(!CVecMat::CalReault(tmpMat,vGlobal,vLocal))
    {
        return(false);
    }

    return(true);
}

/// 局部坐标转到全局坐标
bool GisMath::LOCAL2GLOBAL(double dL, double dB, const CVector &vLocal, CVector &vGlobal)
{
    /// 判断两个矩阵元素是否相同
    if(!CVecMat::IsValid(vLocal,vGlobal))
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
    if(!CVecMat::CalReault(tmpMat,vLocal,vGlobal))
    {
        return(false);
    }

    return(true);
}

/// 将全局坐标转成局部坐标
CMatrix GisMath::GLOBAL2LOCAL(double dL, double dB)
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

/// 局部坐标转成全局坐标
CMatrix GisMath::LOCAL2GLOBAL(double dL, double dB)
{
    return (CVecMat::Transp(GLOBAL2LOCAL(dL,dB)));
}

/// 根据距离 站心 的俯仰、方位角、距离 解算在世界坐标系下的位置
void Local2Globel(double dLon, double dLat, double dAzim, double dElev, double dDist,
                  double &dX, double& dY, double& dZ)
{
    /// 用于存放临时数据
    static CVector vTemp(3);

    /// 计算在站心坐标系下的位置
    /// 因球坐标以逆时针为正方向，而方位角以顺时针为正方向
    /// 故此处对方位角取反；同时球坐标以正东为起始位置，此处
    /// 需要增加90度
    vTemp = CVecMat::VecPolar(DPI*0.5-dAzim,dElev,dDist);

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
    CMatrix M = GisMath::GLOBAL2LOCAL(dLon,dLat);

    /// 用于存放临时数据
    static CVector vTemp(3);
    vTemp(0) = dX;
    vTemp(1) = dY;
    vTemp(2) = dZ;

    /// 将世界坐标转换到站心坐标
    vTemp = M * vTemp;

    /// 计算在站心坐标系的俯仰，方位
    CVecMat::AzEl(vTemp,dAzim,dElev);
    dDist = vTemp.Length();
}

/// 根据给定的地理坐标，俯仰、方位、距离，计算所求点的地理坐标
int GisMath::GeoCalEndGeo(double dLon, double dLat, double dHeight, double dAzim, double dElev, double dDist, double &dLon2, double &dLat2, double &dHeight2)
{
    /// 用于存放世界坐标
    static double dTemp[3];

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
    static double dTemp[3];

    CheckInit();

    /// 根据地理坐标计算世界坐标
    iauGd2gc(EM_TYPE,dLon,dLat,dHeight,dTemp);


    /// 构建偏移向量
    CVector offset(dX,dY,dZ);

    /// 方位旋转
    CMatrix mAzim = CVecMat::R_z(dAzim);

    /// 俯仰旋转
    CMatrix mElev = CVecMat::R_x(-dElev);

    /// 横滚旋转
    CMatrix mRoll = CVecMat::R_y(-dRoll);

    CMatrix tmpMatrix = mAzim*mElev*mRoll;

    offset = tmpMatrix*offset;

    /// 求 东北天 到 地心 的 旋转矩阵
    CMatrix M = LOCAL2GLOBAL(dLon,dLat);

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
    static double dTemp[3],dLon,dLat,dHeight;

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
    static double dTemp1[3],dTemp2[3],dx,dy,dz;

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
    static double dTemp1[3],dTemp2[3],dx,dy,dz;

    CheckInit();

    /// 根据地理坐标计算世界坐标
    iauGd2gc(EM_TYPE,dLon1,dLat1,dHeight1,dTemp1);
    iauGd2gc(EM_TYPE,dLon2,dLat2,dHeight2,dTemp2);

    /// 计算在 世界坐标系下 的相对坐标
    dx = dTemp2[0] - dTemp1[0];
    dy = dTemp2[1] - dTemp1[1];
    dz = dTemp2[2] - dTemp1[2];

    CVector vTemp(dx,dy,dz);

    /// 计算 世界坐标系 到 局部天东北坐标系 旋转矩阵
    CMatrix M = LOCAL2GLOBAL(dLon1,dLat1);

    /// 将世界坐标转换到站心坐标
    vTemp = M * vTemp;

    /// 方位旋转
    CMatrix mAzim = CVecMat::R_z(dAzim);

    /// 俯仰旋转
    CMatrix mElev = CVecMat::R_x(-dElev);

    /// 横滚旋转
    CMatrix mRoll = CVecMat::R_y(-dRoll);


    vTemp = vTemp * mAzim * mElev * mRoll;

    dX = vTemp(0);
    dY = vTemp(1);
    dZ = vTemp(2);

    return(0);
}

/// 根据两个世界坐标计算第二个相对于第一个的方位角、俯仰角
int GisMath::CalAzElXYZ(double dX1, double dY1, double dZ1, double dX2, double dY2, double dZ2, double &dAzim, double &dElev, double &dDist)
{
    static double dTemp[3],dLon,dLat,dHeight,dx,dy;

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
