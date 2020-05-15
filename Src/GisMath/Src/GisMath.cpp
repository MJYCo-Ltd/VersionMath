#include "GisMath.h"
#include "VecMat.h"
#include "../Inc/geodesic.h"

static geod_geodesic PJ_WGS84;
static GisMath::ELLIPSOID  EM_TYPE;
static bool          S_INIT(false);

/*****************************************************
 * 东北天坐标系 笛卡尔坐标系，以物体的重心为坐标原点，物体的正东方向为X轴正方向，
 *            物体重心与地心连线的反方向(天)为Z轴正方向，
 *            根据右手定则，确定Y轴指向(指向正北)。
 *            该坐标系又称站心坐标系，因该坐标常用于地面站观测物体，因此原点又称站心。
 *****************************************************/
#include "sofa.h"

/// 检查是否初始化
inline void CheckInit()
{
    if(!S_INIT)
    {
        GisMath::InitGis(GisMath::WGS_84);
    }
}

/// 初始化地理信息系统
void GisMath::InitGis(ELLIPSOID typeEllipsoid)
{
    static double dA,dF;

    switch (typeEllipsoid)
    {
    case BJ_54:
        dA = 6378245;
        dF = 1/298.3;
        break;
    case WGS_72:
        iauEform(WGS72,&dA,&dF);
        EM_TYPE = typeEllipsoid;
        break;
    case GRS_80:
        iauEform(GRS80,&dA,&dF);
        EM_TYPE = typeEllipsoid;
        break;
    case CGCS_2000:
        dA = 6378137;
        dF = 1.0 / 298.257222101;
        break;
    default:
        iauEform(WGS84,&dA,&dF);
        EM_TYPE = WGS_84;
        break;
    }

    /// 初始化系数
    geod_init(&PJ_WGS84,dA,dF);
    S_INIT = true;
}

void GisMath::LBH2XYZ(double dL, double dB, double dH, double &dX, double &dY, double &dZ)
{
    static double dTemp[3];

    iauGd2gc(WGS84,dL,dB,dH,dTemp);

    dX = dTemp[0];
    dY = dTemp[1];
    dZ = dTemp[2];
}

void GisMath::XYZ2LBH(double dX, double dY, double dZ, double &dL, double &dB, double &dH)
{
    static double dTemp[3];
    dTemp[0] = dX;
    dTemp[1] = dY;
    dTemp[2] = dZ;

    iauGc2gd(WGS84,dTemp,&dL,&dB,&dH);
}

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
        LBH2XYZ(vGeo(i),vGeo(i+1),vGeo(i+2),v3D(i),v3D(i+1),v3D(i+2));
    }

    return(true);
}

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
        XYZ2LBH(v3D(i),v3D(i+1),v3D(i+2),vGeo(i),vGeo(i+1),vGeo(i+2));
    }
    return(true);
}

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

CMatrix GisMath::LOCAL2GLOBAL(double dL, double dB)
{
    return (CVecMat::Transp(GLOBAL2LOCAL(dL,dB)));
}

/**
 * @brief 根据距离 站心 的俯仰、方位角、距离 解算在世界坐标系下的位置
 * @param dLon   经度            [rad]
 * @param dLat   纬度            [rad]
 * @param dAzim  相对于站的方位角  [rad]
 * @param dElev  相对于站的俯仰角  [rad]
 * @param dDist  相对于站的距离    [rad]
 * @param dX     在世界坐标系下X位置 [m]
 * @param dY     在世界坐标系下Y位置 [m]
 * @param dZ     在世界坐标系下Z位置 [m]
 */
inline void Local2Globel(double dLon, double dLat, double dAzim, double dElev, double dDist, double &dX, double& dY, double& dZ)
{
    /// 用于存放临时数据
    static CVector vTemp(3);

    /// 计算在站心坐标系下的位置
    /// 因球坐标以逆时针为正方向，而方位角以顺时针为正方向
    /// 故此处对方位角取反；同时球坐标以正东为起始位置，此处
    /// 需要增加90度
    vTemp = CVecMat::VecPolar(DPI*0.5-dAzim,dElev,dDist);

    static double  dAux;

    /// 计算 世界坐标系 到 局部天东北坐标系 旋转矩阵
    CMatrix M = CVecMat::R_y(-dLat)*CVecMat::R_z(dLon);

    /// 移动矩阵列 形成 世界坐标系 到 站心坐标系 旋转矩阵
    for (int j=0; j<3; ++j)
    {
        dAux=M(0,j);
        M(0,j)=M(1,j);
        M(1,j)=M(2,j);
        M(2,j)= dAux;
    }

    /// 求 东北天 到 地心 的 旋转矩阵
    M = CVecMat::Transp(M);

    /// 在世界坐标下 所求位置 相对于给定位置的位置
    vTemp = M * vTemp;

    dX = vTemp(0);
    dY = vTemp(1);
    dZ = vTemp(2);
}

/**
 * @brief 根据世界坐标系下的位置解算 在站心坐标系下 的俯仰、方位角、距离
 * @param dLon   经度            [rad]
 * @param dLat   纬度            [rad]
 * @param dX     在世界坐标系下X位置 [m]
 * @param dY     在世界坐标系下Y位置 [m]
 * @param dZ     在世界坐标系下Z位置 [m]
 * @param dAzim  相对于站的方位角  [rad]
 * @param dElev  相对于站的俯仰角  [rad]
 * @param dDist  相对于站的距离    [rad]
 */
inline void Globel2Local(double dLon, double dLat, double dX, double dY, double dZ, double &dAzim, double &dElev, double &dDist)
{
    /// 临时变量
    static double  dAux;

    /// 计算 世界坐标系 到 局部天东北坐标系 旋转矩阵
    CMatrix M = CVecMat::R_y(-dLat)*CVecMat::R_z(dLon);

    /// 移动矩阵列 形成 世界坐标系 到 站心坐标系 旋转矩阵
    for (int j=0; j<3; ++j)
    {
        dAux=M(0,j);
        M(0,j)=M(1,j);
        M(1,j)=M(2,j);
        M(2,j)= dAux;
    }

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

    static double  dAux;

    /// 计算 世界坐标系 到 局部天东北坐标系 旋转矩阵
    CMatrix M = CVecMat::R_y(-dLat)*CVecMat::R_z(dLon);

    /// 移动矩阵列 形成 世界坐标系 到 站心坐标系 旋转矩阵
    for (int j=0; j<3; ++j)
    {
        dAux=M(0,j);
        M(0,j)=M(1,j);
        M(1,j)=M(2,j);
        M(2,j)= dAux;
    }

    /// 求 东北天 到 地心 的 旋转矩阵
    M = CVecMat::Transp(M);

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

    /// 临时变量
    static double  dAux;

    /// 计算 世界坐标系 到 局部天东北坐标系 旋转矩阵
    CMatrix M = CVecMat::R_y(-dLat1)*CVecMat::R_z(dLon1);

    /// 移动矩阵列 形成 世界坐标系 到 站心坐标系 旋转矩阵
    for (int j=0; j<3; ++j)
    {
        dAux=M(0,j);
        M(0,j)=M(1,j);
        M(1,j)=M(2,j);
        M(2,j)= dAux;
    }

    /// 将世界坐标转换到站心坐标
    vTemp = M * vTemp;

    /// 方位旋转
    CMatrix mAzim = CVecMat::R_z(dAzim);

    /// 俯仰旋转
    CMatrix mElev = CVecMat::R_x(-dElev);

    /// 横滚旋转
    CMatrix mRoll = CVecMat::R_y(-dRoll);

    CMatrix tmpMatrix = mAzim*mElev*mRoll;

    /// 旋转矩阵的逆为矩阵的转置
    CMatrix endMatrix = CVecMat::Transp(tmpMatrix);

    vTemp = endMatrix*vTemp;

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

/// 计算两个地理坐标的弧长
double GisMath::CalArcLen(double dLon1, double dLat1, double dLon2, double dLat2)
{
    static double dAzim1,dAzim2,dDist;
    CalBaiserF(dLon1,dLat1,dLon2,dLat2,dAzim1,dAzim2,dDist);
    return(dDist);
}

/// 贝塞尔大地主题解算 正解
int  GisMath::CalBaiser(double dLon, double dLat, double dAzim, double dDist, double &dLon2, double &dLat2)
{
    CheckInit();

    /// 弧度转成度
    dLon *= DR2D;
    dLat *= DR2D;
    dAzim *= DR2D;

    geod_direct(&PJ_WGS84,dLat,dLon,dAzim,dDist,&dLat2,&dLon2,0);

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

    geod_inverse(&PJ_WGS84,dLat1,dLon1,dLat2,dLon2,&dDist,&dAzim1,&dAzim2);
    if(dAzim1 < 0)
    {
        dAzim1 += 360;
    }

    if(dAzim2 < 0)
    {
        dAzim2 += 360;
    }

    dAzim1 *= DD2R;
    dAzim2 *= DD2R;

    return(0);
}
