#ifndef MY_GISMATH_H
#define MY_GISMATH_H

#include "GisMath_Global.h"
#include "VecMat.h"

/***
 * rad 表示 弧度
 * deg 表示 度
 * m   表示 米
 ***/

namespace GisMath
{
    /**
     * 椭球类型
     */
    enum ELLIPSOID
    {
        WGS_84=1,
        GRS_80,
        WGS_72,
        BJ_54,
        CGCS_2000
    };

    /**
     * @brief 面积类型
     */
    enum AREA_TYPE
    {
        KMsqu,
        Msqu
    };

    /**
     * @brief 长度单位
     */
    enum LENGTH_TYPE
    {
        KM,
        M,
        NMI     /// 海里
    };

    /**
     * @brief 初始化地理信息系统
     * @attention 调用坐标转换之前进行初始化，未初始化将按照WGS_84处理
     * @param typeEllipsoid 椭球类型
     */
    GISMATHSHARED_EXPORT void InitGis(ELLIPSOID typeEllipsoid);

    /**
     * @brief 计算椭球和射线的交点
     * @param pt            射线起点 [m,m,m]
     * @param stDir         射线方向 相对于起点的一个全局坐标系下的方向向量
     * @param rInsertPos    交点     [m,m,m]
     * @return
     */
    GISMATHSHARED_EXPORT bool CalLineInterEllipsoid(const CVector& pt,const CVector& stDir,CVector& rInsertPos);

    /************************ 表述方式的转换 ************************/
    /**
     * @brief 地心笛卡尔坐标 转 地理坐标
     * @param dX   X轴坐标[m]
     * @param dY   Y轴坐标[m]
     * @param dZ   Z轴坐标[m]
     * @param dL   经度[rad]
     * @param dB   纬度[rad]
     * @param dH   高度[m]
     */
    GISMATHSHARED_EXPORT bool XYZ2LBH(double dX, double dY, double dZ, double& dL, double& dB, double& dH);
    GISMATHSHARED_EXPORT bool LBH2XYZ(double dL, double dB, double dH, double& dX, double& dY, double& dZ);

    /// @attention 用向量表示(适合数据为3的倍数) 支持批量处理
    GISMATHSHARED_EXPORT bool XYZ2LBH(const CVector& v3D, CVector& vGeo);
    GISMATHSHARED_EXPORT bool LBH2XYZ(const CVector& vGeo, CVector& v3D);

    /************************ 观察方式的转换 ************************/
    /**
     * @brief 全局坐标转换成局部东北天坐标
     * @param dL       站的经度            [rad]
     * @param dB       站的纬度            [rad]
     * @param vGlobal  全局(ECF)下笛卡尔坐标 [m,m,m]
     * @param vLocal   站心坐标系系的坐标    [m,m,m]
     */
    GISMATHSHARED_EXPORT bool GLOBAL2LOCAL(double dL, double dB,const CVector& vGlobal, CVector& vLocal);
    GISMATHSHARED_EXPORT bool LOCAL2GLOBAL(double dL, double dB,const CVector& vLocal, CVector& vGlobal);

    /**
     * @brief 全局坐标转换成局部东北天坐标转换的矩阵
     * @param dL 站的经度 [rad]
     * @param dB 站的纬度 [rad]
     * @attention vLocal = GLOBAL2LOCAL(dLon,dLat) * vGlobal;
     * @return
     */
    GISMATHSHARED_EXPORT CMatrix GLOBAL2LOCAL(double dL, double dB);
    GISMATHSHARED_EXPORT CMatrix LOCAL2GLOBAL(double dL, double dB);


    /**
     * @brief 根据给定的地理坐标，俯仰、方位、距离，计算所求点的地理坐标
     * @attention 此处的距离为 直线距离
     *
     * @param dLon      给定点经度      [rad]
     * @param dLat      给定点纬度      [rad]
     * @param dHeight   给定点高度      [m]
     * @param dAzim     方位角         [rad]
     * @param dElev     俯仰角         [rad]
     * @param dDist     距离给定点的直线距离 [m]
     * @param dLon2     所求点的经度    [rad]
     * @param dLat2     所求点的纬度    [rad]
     * @param dHeight2  所求点的高度    [m]
     * @return 返回值为 0 表示计算正确
     */
    GISMATHSHARED_EXPORT int GeoCalEndGeo(double dLon, double dLat, double dHeight
                                        , double dAzim, double dElev, double dDist
                                        , double& dLon2, double& dLat2, double &dHeight2);

    /**
     * @brief 根据给定的地理坐标，俯仰、方位、距离，计算所求点的地理坐标
     * @attention 此处的距离为 直线距离
     *
     * @param dLon      给定点经度      [rad]
     * @param dLat      给定点纬度      [rad]
     * @param dHeight   给定点高度      [m]
     * @param dAzim     方位角         [rad]
     * @param dElev     俯仰角         [rad]
     * @param dRoll      翻滚角         [rad]
     * @param dDist     距离给定点的直线距离 [m]
     * @param dLon2     所求点的经度    [rad]
     * @param dLat2     所求点的纬度    [rad]
     * @param dHeight2  所求点的高度    [m]
     * @return 返回值为 0 表示计算正确
     */
    GISMATHSHARED_EXPORT int GeoCalEndGeo(double dLon,double dLat,double dHeight
                                        ,double dAzim,double dElev,double dRoll, double dX,double dY, double dZ
                                        ,double& dLon2, double& dLat2,double& dHeight2);

    /**
     * @brief 根据给定的世界坐标，俯仰、方位、距离，计算所求点的世界坐标
     * @param dX       给定点X轴分量     [m]
     * @param dY       给定点Y轴分量     [m]
     * @param dZ       给定点Z轴分量     [m]
     * @param dAzim    方位角           [rad]
     * @param dElev    俯仰角           [rad]
     * @param dDist    距离给定点的直线距离 [m]
     * @param dX2      所求点的X轴分量    [m]
     * @param dY2      所求点的Y轴分量    [m]
     * @param dZ2      所求点的Z轴分量    [m]
     * @return 返回值为 0 表示计算正确
     */
    GISMATHSHARED_EXPORT int XYZCalEndXYZ(double dX, double dY, double dZ
                                       , double dAzim, double dElev, double dDist
                                      , double& dX2, double& dY2, double &dZ2);

    /**
     * @brief 根据两个地理坐标计算第二个相对于第一个的方位角、俯仰角、距离
     * @param dLon1        第1个坐标的经度   [rad]
     * @param dLat1         第1个坐标的纬度   [rad]
     * @param dHeight1    第1个坐标的高度   [m]
     * @param dLon2        第2个坐标的经度   [rad]
     * @param dLat2         第2个坐标的纬度   [rad]
     * @param dHeight2    第2个坐标的高度   [m]
     * @param dAzim        第2个点相对于第1个点的方位角 [rad]
     * @param dElev         第2个点相对于第1个点的俯仰角 [rad]
     * @param dDist         第2个点相对于第1个点的直线距离 [m]
     * @return
     */
    GISMATHSHARED_EXPORT int CalAzElGeo(double dLon1, double dLat1, double dHeight1
                                  , double dLon2, double dLat2, double dHeight2,
                                    double &dAzim, double &dElev,double &dDist);

    GISMATHSHARED_EXPORT int CalAzElGeo(double dLon1,double dLat1,double dHeight1
                                  ,double dLon2,double dLat2,double dHeight2,
                                   double dAzim, double dElev, double dRoll,
                                    double &dX, double& dY, double &dZ);

    /**
     * @brief 根据两个世界坐标计算第二个相对于第一个的方位角、俯仰角
     * @param dX1      第1个点的X轴分量   [m]
     * @param dY1      第1个点的Y轴分量   [m]
     * @param dZ1      第1个点的Z轴分量   [m]
     * @param dX1      第2个点的X轴分量   [m]
     * @param dY1      第2个点的Y轴分量   [m]
     * @param dZ1      第2个点的Z轴分量   [m]
     * @param dAzim  第2个点相对于第1个点的方位角 [rad]
     * @param dElev    第2个点相对于第1个点的俯仰角 [rad]
     * @param dDist    第2个点相对于第1个点的直线距离 [m]
     * @return
     */
    GISMATHSHARED_EXPORT int CalAzElXYZ(double dX1, double dY1, double dZ1
                                  , double dX2, double dY2, double dZ2,
                                    double &dAzim, double &dElev, double& dDist);

    /**
     * @brief 计算两个地理坐标的球面弧长
     * @param dLon1   第1个点的经度 [rad][-π,π)
     * @param dLat1    第1个点的纬度 [rad][-π/2,π/2]
     * @param dLon2   第2个点的经度 [rad][-π,π)
     * @param dLat2    第2个点的纬度 [rad][-π/2,π/2]
     * @return 球面弧长  [m]
     */
    GISMATHSHARED_EXPORT double CalArcLen(double dLon1, double dLat1,
                                          double dLon2, double dLat2,LENGTH_TYPE lenTYpe);

    /**
     * @brief 贝塞尔大地主题解算 正解
     * @attention 根据给定的经度、纬度，及相对于该经度纬度的
     *                   方位、距离，求另一个点的经纬度，
     *                  距离：椭球的大地线长度
     * @param dLon      给定点经度      [rad][-π,π)
     * @param dLat      给定点纬度      [rad][-π/2,π/2]
     * @param dAzim     方位角             [rad][0,2π)
     * @param dDist     距离给定点的距离 [m]
     * @param dLon2     所求的经度值    [rad][-π,π)
     * @param dLat2     所求的纬度值   [rad][-π/2,π/2]
     * @return 0 表示没有错误
     */
    GISMATHSHARED_EXPORT int CalBaiser(double dLon, double dLat, double dAzim
                               , double dDist, double &dLon2, double &dLat2);

    /**
     * @brief 贝塞尔大地主题解算 反解
     * @attention 方位角均为真方位角
     *                   考虑的是椭球
     * @param dLon1   第1个点的经度 [rad][-π,π)
     * @param dLat1    第1个点的纬度 [rad][-π/2,π/2]
     * @param dLon2   第2个点的经度 [rad][-π,π)
     * @param dLat2    第2个点的纬度 [rad][-π/2,π/2]
     * @param dAzim1 第2个点相对于第1个点的方位角 [rad][0,2π)
     * @param dAzim2 第1个点相对于第2个点的方位角 [rad][0,2π)
     * @param dDist    两点之间的大地线长度 [m]
     * @return 0 表示没有错误
     */
    GISMATHSHARED_EXPORT int CalBaiserF(double dLon1, double dLat1, double dLon2, double dLat2
                                 , double& dAzim1, double &dAzim2, double &dDist);

    /**
     * @brief 根据经纬度计算面积
     * @param vGeoIn  经纬度类型(lon,lat,lon,lat……)[deg,deg,deg,deg,……]
     * @param type    输出面积类型
     * @return
     */
    GISMATHSHARED_EXPORT double CalclutaGeoArea(const CVector& vGeoIn,AREA_TYPE type);
}

#endif // GISMATH_H
