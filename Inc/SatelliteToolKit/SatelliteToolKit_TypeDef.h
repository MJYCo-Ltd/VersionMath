#ifndef SATELLITETOOLKIT_TYPEDEF_H
#define SATELLITETOOLKIT_TYPEDEF_H

/******************************************************
 *  本库 所有的年月日时分秒 是 北京时间      【东8区】
 *      所有的double     是  UTC约简儒略日 【0时区】
 ******************************************************/

#include<vector>
using namespace std;

/**
 * 绕X轴旋转,又称滚动旋转Roll,简写R (-π~π]
 * 绕Y轴旋转,又称俯仰旋转Picth,简写P (-π~π]
 * 绕Z轴旋转,又称方位旋转Yaw,简写Y [0~2π)
 * 每一行的两个旋转顺序类型是等同的
 */
enum RotateType
{
    Rota_PRY,Rota_YXZ,
    Rota_PYR,Rota_YZX,
    Rota_RYP,Rota_XZY,
    Rota_RPY,Rota_XYZ,
    Rota_YPR,Rota_ZYX,
    Rota_YRP,Rota_ZXY
};

/**
 * @brief 卫星根数类型
 */
enum SatElementType
{
    SAT_TLE,     /// TLE
    SAT_TWOBODY, /// 六根数
    SAT_PV       /// 位置速度
};

/**
 * @brief 形状类型
 */
enum ShapeType
{
    eRectangle,  /// 矩形
    eEllipse     /// 椭圆
};

struct BJTime   /// 北京时间
{	
	BJTime()
	{
		nYear = 0;
		uMonth = uDay = uHour = uMinute = uSecond = 0;
		uMSecond = 0;
	}
	BJTime(short _nYear, unsigned char _uMonth, unsigned char _uDay, unsigned char _uHour, unsigned char _uMinute, unsigned char _uSecond, unsigned short _uMSecond)
	{
		nYear = _nYear; uMonth = _uMonth; uDay = _uDay; uHour = _uHour; uMinute = _uMinute; uSecond = _uSecond, uMSecond = _uMSecond;
	}
    short         nYear;
    unsigned char uMonth;
    unsigned char uDay;
    unsigned char uHour;
    unsigned char uMinute;
    unsigned char uSecond;
    unsigned short uMSecond;  /// 毫秒
};

struct Pos /// 位置结构体
{
    double dX;
    double dY;
    double dZ;
};

struct PV  /// 位置速度结构体
{
    Pos stP; /// 位置
    Pos stV; /// 速度
};


struct Kepler   /// 开普勒六根数
{
    BJTime stEpoch;  /// 历元时间
    double dA;       ///  长半轴  [米]
    double dE;       /// 离心率
    double dI;        /// 轨道倾角 [弧度]
    double dRAAN;   /// 升交点赤经 [弧度]
    double dW;       /// 近地点幅角 [弧度]
    double dMA;      /// 平近点角 [弧度]
};

struct SatellitePV /// 卫星某一时刻的状态
{
    BJTime stEpoch;  /// 历元时间
    PV     stPV;     /// ECI下卫星的位置和速度[m,m,m m/s,m/s,m/s]
};

struct SatelliteTle/// 卫星的TLE数据
{
    string sLine1;
    string sLine2;
};


struct Satellite_Element /// 卫星根数
{
    SatElementType elemType; /// 根数类型
    SatellitePV    stSatPV;  /// 卫星瞬时状态
    Kepler         stKepler; /// 开普勒六根数
    SatelliteTle   stTLE;    /// tle数据
};


struct SatellitePos /// 卫星位置结构体
{
    BJTime stStart;
    BJTime stEnd;
    unsigned int uStepMs;	/// 计算步长 [毫秒]
    vector<double> vTimes;	/// 约简儒略日
    vector<PV> vJ2000;		/// J2000下位置
    vector<PV> vECF;		/// 地固系下位置坐标
    vector<Pos> vLLA;		/// 经纬高 [deg deg m]
};

struct Period /// 时间段结构体
{
    BJTime stStart;       ///开始时间
    BJTime stEnd;         ///结束时间
    double dDurationTime; ///持续时间
    vector<Pos> vLLa;      /// 该段时间对应的地面区域范围
};

#endif
