#ifndef GISMATH_COMMON_H
#define GISMATH_COMMON_H
#include "geodesic.h"
#include <VersionMathCommon.h>
#include <GisMath/GisMath.h>

#ifdef __cplusplus
extern "C" {
#endif
double iauAnp(double a);
int iauEform(int n, double *a, double *f);
int iauGd2gc(int n, double elong, double phi, double height,
             double xyz[3] );
int iauGc2gd(int n, double xyz[3],
             double *elong, double *phi, double *height);
#ifdef __cplusplus
}
#endif

extern geod_geodesic       PJ_GEOD;
extern GisMath::ELLIPSOID  EM_TYPE;
extern bool                S_INIT;
extern double              EARTH_DA;
extern double              EARTH_DB;
extern double              EARTH_DF;
const double        M2KM(1e-3);    /// 米到千米的转换
const double        M2NMI(1/1852.);/// 米到海里的转换
const double        Msqu2KMsqu(1e-6);  /// 平方米到平方千米的转换

/**
 * @brief 判断是否初始化
 */
void CheckInit();

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
void Local2Globel(double dLon, double dLat, double dAzim, double dElev, double dDist,
                  double &dX, double& dY, double& dZ);

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
void Globel2Local(double dLon, double dLat, double dX, double dY, double dZ,
                  double &dAzim, double &dElev, double &dDist);

#endif // GISMATH_COMMON_H
