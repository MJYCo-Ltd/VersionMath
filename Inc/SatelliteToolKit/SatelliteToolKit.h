#ifndef SATELLITETOOLKIT_H
#define SATELLITETOOLKIT_H

#include<string>
#include<vector>
#include "SatelliteToolKit_global.h"
#include "SatelliteToolKit_TypeDef.h"

using namespace  std;

namespace Math
{
    class CMatrix;
}

/**
 * @brief 初始化卫星工具包的环境
 * @param sErrorInfo 错误信息
 * @return
 */
bool SATELLITETOOLKIT_EXPORT InitSatelliteToolKit(string& sErrorInfo);

/**
 * @brief 释放卫星工具包占用的资源
 */
void SATELLITETOOLKIT_EXPORT CloseSatelliteToolKit();

/**
 * @brief 根据开普勒六根数计算卫星轨道
 * @param stStartTime 开始计算时间
 * @param stEndTime   结束时间
 * @param nStep       计算步长 [秒]
 * @param stKepler    开普勒六根数和根数的历元时间
 * @param stSatPos    计算出来的卫星轨道数据
 * @return
 */
bool SATELLITETOOLKIT_EXPORT TwoBody(const BJTime& stStartTime,
                                     const BJTime& stEndTime,
                                     unsigned int nStep,
                                     const Kepler&  stKepler,
                                     SatellitePos& stSatPos);

/**
 * @brief 根据两行星历计算卫星轨道
 * @param stStartTime 开始计算时间
 * @param stEndTime   结束时间
 * @param nStep       计算步长 [秒]
 * @param sLine1       第一行星历
 * @param sLine2       第二行星历
 * @param stSatPos    计算出来的卫星轨道数据
 * @return
 */
bool SATELLITETOOLKIT_EXPORT SGP4(const BJTime& stStartTime,
                                  const BJTime& stEndTime,
                                  unsigned int nStep,
                                  const string&  sLine1,
                                  const string&  sLine2,
                                  SatellitePos& stSatPos);

/**
 * @brief 计算满足要求的太阳高度角的时段
 * @param stStartTime         开始时间
 * @param stEndTime           结束时间
 * @param fIncidentMinAngle   最小入射角度[弧度][0~fIncidentMaxAngle)
 * @param fIncidentMaxAngle   最大入射角度[弧度](fIncidentMinAngle,π/2]
 * @param stGeo               地面目标位置
 *                            经度、纬度、高度 [弧度、弧度、米]
 * @return
 */
vector<Period> SATELLITETOOLKIT_EXPORT SolarAltitude(const BJTime& stStartTime,
                                                     const BJTime& stEndTime,
                                                     float   fIncidentMinAngle,
                                                     float   fIncidentMaxAngle,
                                                     const Pos&   stGeo);

/**
 * @brief 计算太阳高度角
 * @param stTime  北京时间
 * @param stGeo   地面目标位置 经度、纬度、高度 [弧度、弧度、米]
 * @return
 */
double SATELLITETOOLKIT_EXPORT CalSolarAltitude(const BJTime& stTime,const Pos& stGeo);

/**
 * @brief 计算卫星相对地面可见弧段
 * @param stStartTime  开始时间     [北京时间]
 * @param stEndTime    结束时间     [北京时间]
 * @param stSatElement  卫星根数数据
 * @param stPos        地面站位置 经度、纬度、高度 [rad、rad、m]
 * @param emRotate     载荷旋转顺序
 * @param stPRY        载荷旋转角度  [rad,rad,rad][绕X轴旋转,绕Y轴旋转,绕Z轴旋转]
 * @param fHAngle      水平视角     [rad] X轴方向
 * @param fVangle      垂直视角     [rad] Y轴方向
 * @param emType       形状类型
 * @return
 */
vector<Period> SATELLITETOOLKIT_EXPORT VisiblePeriod(const BJTime& stStartTime,
                                                     const BJTime& stEndTime,
                                                     const Satellite_Element& stSatElement,
                                                     const Pos&   stPos,
                                                     RotateType emRotate,
                                                     const Pos& satPRY,
                                                     double  dHAngle,
                                                     double  dVAngle,
                                                     ShapeType emType);

/**
 * @brief 判断卫星和地面站是否可见
 * @param satECFPos          卫星位置    wgs84 坐标[m,m,m]
 * @param station3D          地面站位置  wgs84 坐标[m,m,m]
 * @param fVisibleAngle   最小可见角度  [弧度][0~π/2]
 * @return
 */
SATELLITETOOLKIT_EXPORT bool IsVisible(const Pos& satECFPos, const Pos& station3D, double dVisibleAngle);

/**
 * @brief 判断传感器是否可视
 * @param satPV           卫星的位置速度[ECF m,m,m m/s,m/s,m/s]
 * @param stationPos      地面位置     [m,m,m]
 * @param satPRY          卫星的姿态   [弧度,弧度,弧度][绕X轴旋转,绕Y轴旋转,绕Z轴旋转]
 * @param fVisibleHAngle  水平方向上的角度   [弧度] X轴方向
 * @param fVisibleVAngle  垂直方向上的角度   [弧度] Y轴方向
 * @param eRotate         旋转顺序
 * @return
 */
SATELLITETOOLKIT_EXPORT bool EllipseVisible(const PV& satECFPV, const Pos& station3DPos, const Pos& satPRY,
                                       double dVisibleHAngle, double dVisibleVAngle, RotateType eRotate);
SATELLITETOOLKIT_EXPORT bool RectangleVisible(const PV& satECFPV, const Pos& station3DPos, const Pos& satPRY,
                                              double dVisibleHAngle, double dVisibleVAngle, RotateType eRotate);

/**
 * @brief Target相对于Eye可见时段
 * @param stSatPosEye    观测星
 * @param stSatPosTarget 被观测星
 * @param fVisibleAngle  观测星的半角[弧度](0~π/2)
 * @return
 */
vector<Period> SATELLITETOOLKIT_EXPORT VisiblePeriod(const SatellitePos& stSatPosEye,
                                                     const SatellitePos& stSatPosTarget,
                                                     float fVisibleAngle);

/**
 * @brief 计算某位置与地面相交区
 * @param satPV     卫星的位置速度[ECF m,m,m m/s,m/s,m/s]
 * @param satPRY    卫星的姿态   [rad,rad,rad][绕X轴旋转,绕Y轴旋转,绕Z轴旋转]
 * @param vAngle    水平、垂直方向的半角 [弧度(0~π/2)，弧度(0~π/2)]
 * @return
 */
vector<Pos> SATELLITETOOLKIT_EXPORT Intersect(const PV& satPV,
                                              const Math::CMatrix &rotateMatrix,
                                              const vector<Pos>& vAngle);

/**
 * @brief 计算某位置圆锥与地面相交区
 * @param stPos     地球外面某一点ECF位置 x,y,z[米，米，米]
 * @param fPitch    俯仰角[弧度](-π~π]
 * @param fRoll     翻滚角[弧度](-π~π]
 * @param fRaw      方位角[弧度][0~2π)
 * @param fAngle    圆锥半角[弧度](0~π/2)
 * @return
 */
vector<Pos> SATELLITETOOLKIT_EXPORT IntersectCircle(const PV& satPV,
                                                    const Pos& satPRY,
                                                    RotateType eRotate,
                                                    double dHAngle,
                                                    double dVAngle);

/**
 * @brief 卫星圆形载荷过某一点的时间段
 * @param stSatPos     卫星位置数据
 * @param stPos        地面位置 经度、纬度、高度 [弧度、弧度、米]
 * @param fAngle       圆锥半角[弧度](0~π/2)
 * @return
 */
vector<Period> SATELLITETOOLKIT_EXPORT IntersectCircle(const SatellitePos& stSatPos,
                                                       const Pos& stPos,
                                                       float fAngle);

/**
 * @brief 计算某位置方锥与地面相交区
 * @param stPos     地球外面某一点ECF位置 x,y,z[米，米，米]
 * @param fHAngle   方锥水平半角[弧度](0~π/2)
 * @param fVAngle   方锥垂直半角[弧度](0~π/2)
 * @return
 */
vector<Pos> SATELLITETOOLKIT_EXPORT IntersectRectangle(const PV& satPV,
                                                       const Pos& satPRY,
                                                       RotateType eRotate,
                                                       double dHAngle,
                                                       double dVAngle);

/**
 * @brief 根据地面位置计算卫星的姿态
 * @param satPV       卫星的位置速度[ECF m,m,m m/s,m/s,m/s]
 * @param rPos        地面位置     [m,m,m]
 * @param eRoteType   旋转顺序
 * @param satPRY      所求的姿态值
 * @return            可以计算姿态值为true，不可以计算为false
 */
bool SATELLITETOOLKIT_EXPORT CalPRY(const PV& satPV,
                                    const Pos& rPos,
                                    RotateType eRoteType,
                                    Pos& satPRY);

/**
 * @brief 卫星方形载荷过某一点的时间段
 * @param stSatPos     卫星位置数据
 * @param stPos        地面位置 经度、纬度、高度 [弧度、弧度、米]
 * @param fHAngle      方锥水平半角[弧度](0~π/2)
 * @param fVAngle      方锥垂直半角[弧度](0~π/2)
 * @return
 */
vector<Period> SATELLITETOOLKIT_EXPORT IntersectRectangle(const SatellitePos& stSatPos,
                                                       const Pos& stPos,
                                                       float fHAngle,
                                                       float fVAngle);

vector<Satellite_Element> SATELLITETOOLKIT_EXPORT CreateConstellatory(Satellite_Element satTemplet,
                                              int nPlanes,
                                              int nNumSats);

#endif
