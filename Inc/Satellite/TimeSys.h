#ifndef YTY_TIMESYS_H
#define YTY_TIMESYS_H

/*****************************************
  作用：完成各种时间的相互转换
  注意：此类是非线程安全的，在多线程中使用需要对变量
       进行加锁保护，否则运算结果不能保证
名词解释：
        TCB(Barycentric Coordinate Time) 质心坐标时
        TDB(Barycentric Dynamical  Time) 质心力学时
        TCG(Geocentric Coordinate  Time) 地心坐标时
        TAI(International Atomic   Time) 国际原子时(GPS也是原子时钟的)
        TT(Terrestrial             Time) 地球时
        UTC(Coordinated Universal  Time) 世界协调时(我们生活中使用的时间)

        ///不再使用的时间系统
        ET(Ephemeris               Time) 历书时(1984年以前用作太阳系引力理论自变量的时间尺度，
        其单位和零点有协议的定义，已为TT 和 TDB所代替)
        TDT(Terrestrial Dynamical  Time) 地球力学时(1979年IAU决议定义的视地心历表使用的时
        间，1991年为地球时（TT）所代替)
 *****************************************/
#include "SAT_global.h"
#include "Date.h"
namespace Aerospace{
class ALGORITHM_EXPORT CTimeSys
{
public:
    CTimeSys();
    CTimeSys(double dMjdUTC);
    CTimeSys(const CDate& rDate);

    /**
     * @brief 获取约简儒略日格式的 TT 时间
     */
    double GetTT() const;

    /**
     * @brief 获取约简儒略日格式的 TAI 时间
     */
    double GetTAI() const;

    /**
     * @brief 获取约简儒略日格式的 UT1 时间
     */
    double GetUT1() const;

    /**
     * @brief 获取约简儒略日格式的 TDB 时间
     */
    double GetTDB() const;

    /**
     * @brief 获取约简儒略日格式的 TCB 时间
     */
    double GetTCB() const;

    /**
     * @brief 获取约简儒略日格式的 TCG 时间
     */
    double GetTCG() const;

    /**
     * @brief 获取约简儒略日格式的 GPS 时间
     */
    double GetGPS() const;

private:
    CDate m_dateUtc;
};
}
#endif // YTY_TIMESYS_H
