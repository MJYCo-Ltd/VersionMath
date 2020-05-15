#ifndef YTY_JPL_H
#define YTY_JPL_H

/*****************************************
  作用：本文件负责读取DExxx 文件，并解析文件内容，
       进行坐标运算
  注意：此类是非线程安全的，在多线程中使用需要对变量
       进行加锁保护，否则运算结果不能保证
  结构：采用单例设计模式进行设计
数据获取：地址：ftp://ssd.jpl.nasa.gov/；如果地址不可达
        则搜索"JPL星历表"关键字
   注意：
        DE118 采用的是 B1950 参考系
        DE2XX 采用的是 J2000 参考系
        DE4XX 采用的是 ICRF  参考系
 *****************************************/

#include <string>
using namespace std;
#include "SAT_global.h"

namespace Math{
class CVector;
}

namespace Aerospace{
using namespace Math;
enum PLANET_TYPE
{
    Mercury 				= 1, /// 水星
    Venus 					= 2, /// 金星
    Earth 					= 3, /// 地球
    Mars 					= 4, /// 火星
    Jupiter 				= 5, /// 木星
    Saturn 					= 6, /// 土星
    Uranus 					= 7, /// 天王星
    Neptune 				= 8, /// 海王星
    Pluto 					= 9, /// 冥王星
    Moon 					= 10,/// 月球
    Sun 					= 11,/// 太阳
    SolarSystemBarycenter 	= 12,/// 太阳系质心
    EarthMoonBarycenter 	= 13,/// 地月质心
    Nutation 				= 14,/// 章动数据
    Libration 				= 15 /// 月球平动
};

class ALGORITHM_EXPORT CJPLEphemeris
{
public:

    /**
     * @brief 获取单例指针
     * @return 返回单例指针
     */
    static CJPLEphemeris* GetInstance();

    /**
     * @brief 释放单例占用的空间
     */
    static void Realse();

    /**
     * @brief 初始化
     * @param sFileName DExxx文件名称
     * @return
     */
    bool Init(const string& sFileName);

    /**
     * @brief 查看星历表是否初始化
     * @return
     */
    inline bool IsInit()const {return (m_bInit);}

    /**
     * @brief 获取太阳在以地球为原点的坐标系下的位置
     * @param dMjdTT [约简儒略日TT时间]
     * @return xyz [m][ECI]
     */
    CVector GetSunPos(const double& dMjdTT) const;

    /**
     * @brief 获取月球在以地球为原点的坐标系下的位置
     * @param dMjdTT [约简儒略日TT时间]
     * @return xyz [m][ECI]
     */
    CVector GetMoonPos(const double& dMjdTT) const;

    /**
     * @brief 获取任一行星以其他行星为原点的坐标系下的位置
     * @param dMjdTT    [MJD TT]
     * @param planet1 所要计算的行星
     * @param centerPlanet 中心行星
     * @return (x,y,z) [m][ECI]
     */
    CVector GetPos(const double& dMjdTT, const PLANET_TYPE &planet1
                   , const PLANET_TYPE& centerPlanet)const;

    /**
     * @brief 获取地球到太阳的平均距离
     * @return [m]
     */
    double GetAu(){return(m_dAu);}

private:
    /// 构造函数
    CJPLEphemeris();
    /// 析构函数
    ~CJPLEphemeris();

private:
    void*  m_pJplDate; /// 存放jpl星历表数据
    bool   m_bInit;    /// 保存初始化信息
    double m_dAu;      /// 太阳到地球的距离
    static CJPLEphemeris* m_pInstance;/// 单例指针
};
}
#endif // YTY_JPL_H
