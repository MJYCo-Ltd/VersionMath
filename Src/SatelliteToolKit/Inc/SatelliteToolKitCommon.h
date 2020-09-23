#ifndef SATELLITETOOLKITCOMMON_H
#define SATELLITETOOLKITCOMMON_H
#include <vector>
#include <VecMat.h>
/**
 * @brief 时间和角度
 */
struct TimeElev
{
    double dElev;
    double dMJDTime;
};

struct STK_Point
{
    double dX;
    double dY;
};

using namespace std;
using namespace Math;

#include "SatelliteToolKit_TypeDef.h"
/**
 * @brief 判断时间是否有效
 * @param stStartTime 开始有效
 * @param stEndTime   结束有效
 * @param dMJDStart   开始的约简儒略日
 * @param dMJDEnd     结束的约简儒略日
 * @return
 */
bool JudgeTimeValid(const BJTime &stStartTime, const BJTime &stEndTime,
                    double& dMJDStart, double& dMJDEnd);

/**
 * @brief 判断数据的趋势
 * @return
 */
int  JudgeDataTrend(const vector<TimeElev>& vElev, float dMin);

/**
 * @brief 判断点是否与载荷相交
 * @param satECFPV
 * @param station3DPos
 * @param satPRY
 * @param eRotate
 * @return
 */
bool JudgeIsInsert(const PV& satECFPV, const Pos& station3DPos, const Pos& satPRY,
                   RotateType eRotate, CVector& vInsert, double &dAngle);

///**
// * @brief 判断两点是否被地球遮挡
// * @param rPos1
// * @param rPos2
// * @return
// */
//bool  InsertEarth(const Pos& rPos1,const Pos& rPos2);
//bool  InsertEarth(const CVector& vR1, const CVector& vR2);

/**
 * @brief 计算卫星的旋转矩阵
 * @param satPV       卫星位置速度
 * @param satPRY      卫星姿态
 * @param eRotate     卫星旋转顺序
 * @return
 */
CMatrix CalSatRoteMatrix(const PV& satPV, const Pos& satPRY, RotateType eRotate);

///**
// * @brief 通过卫星的位置姿态计算旋转矩阵
// * @param satPV 卫星的位置 速度
// * @return
// */
//CMatrix CalSatMatrix(const PV& satPV);
//CMatrix CalSatMatrix(const CVector& vPV);

/**
 * @brief 通过旋转模型计算旋转矩阵
 * @param satPRY
 * @param eROtate
 * @return
 */
CMatrix CalRotateMatrix(const Pos& satPRY, RotateType eRotate);

/**
 * @brief 计算太阳高度角
 * @param dMJD               约简儒略日
 * @param stGruondStation    地面站经纬度信息
 * @return
 */
double CalSolarAltitude(double dMJD, const Pos& gruondStationGeo);

/**
 * @brief 计算二维向量的叉乘
 * @param rP1
 * @param rP2
 * @param rP
 * @return
 */
double CalCross(const STK_Point& rP1, const STK_Point& rP2, const STK_Point& rP);

/**
 * @brief 计算射线和椭球相交
 * @param pt           射线的起始点
 * @param stDir        射线的方向
 * @param rInsertPos   射线与椭球的交点
 * @return 相交返回true不相交返回false
 */
bool CalLineInterEllipsoid(const Pos& pt,const Pos& stDir,Pos& rInsertPos);

/**
 * @brief 是否相交
 * @param satPos    卫星位置
 * @param rDir      旋转后载荷指向
 * @param dMaxAngle 载荷最大的角度
 * @return
 */
bool CanIntersection(const Pos& satPos, const CVector &rDir, const double& dMaxAngle);

/**
 * @brief 判断角度是否在椭圆内
 * @param dHAngle
 * @param dVAngle
 * @param dAngle
 * @return
 */
bool EllipseVisible(const double& dVisibleHAngle, const double& dVisibleVAngle, const double& dAngle, const CVector& vInsert);

/**
 * @brief 判断角度是否在矩形内
 * @param dVisibleHAngle
 * @param dVisibleVAngle
 * @param dAngle
 * @param vInsert
 * @return
 */
bool RectangleVisible(const double& dVisibleHAngle, const double& dVisibleVAngle, const double& dAngle, const CVector& vInsert);
#endif // SATELLITETOOLKITCOMMON_H
