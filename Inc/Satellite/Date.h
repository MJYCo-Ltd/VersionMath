#ifndef YTY_DATE_H
#define YTY_DATE_H

/*****************************************
  作用：进行时间的运算，采用约简儒略日(UTC时间)，
       不建议通过CDate保存北京时间，因为轨道
       运算采用的都是UTC时间。
  注意：此类是非线程安全的，在多线程中使用需要对变量
       进行加锁保护，否则运算结果不能保证
非法数据：简约儒略日小于0的数据认为非法
 *****************************************/

#include <iostream>
#include "SAT_global.h"

namespace Aerospace{
enum DATE_TYPE
{
    UTC,BJ
};

class ALGORITHM_EXPORT CDate
{
public:
    /// 默认构造函数
    CDate(DATE_TYPE emType=UTC);
    CDate(int nYear, int nMonth, int nDay
          ,int nHour, int nMinute, double dSecond
          ,DATE_TYPE emType=BJ);
    CDate(double dMjd,DATE_TYPE emType=UTC);
    /// 复制构造
    CDate(const CDate& rDate);

    /// 设置时区
    void   SetDateType(const DATE_TYPE& emType);
    inline DATE_TYPE GetDateType()const{return (m_emType);}

    /// 获取日历信息
    bool GetDate(int& nYear, int& nMonth, int& nDay
                 ,int& nHour, int& nMinute, double& dSecond) const;

    /// 设置获取秒小数点后的有效位数
    /// @param nDecimal 范围 [1~9]
    void SetSecondDecimal(int nDecimal);
    inline int  GetSecondDecimal() const {return (m_nDecimal);}

    /// 获取UTC时间的约简儒略日
    /// 获取北京时间的约简儒略日，没有意义
    /// 因此不提供获取北京时间的约简儒略日
    inline double GetMJD() const{return(m_dMjd);}

    /// 设置/获取儒略日时间
    /// @Param emType 儒略日表示的类型
    void   SetJD(double dJD,DATE_TYPE emType=UTC);
    double GetJD()const;

    /// 重载操作符
    /// @param dTime 用天数表示的时间
    void operator +=(double dTime);
    void operator -=(double dTime);
    bool operator ==(const CDate& rDate);
    bool operator !=(const CDate& rDate);
    operator bool()const{return (0 < m_dMjd);}

    /// 重载赋值操作符
    CDate& operator=(const double dMjd);
    CDate& operator=(const CDate& rDate);
private:
    double    m_dMjd;    /// 约简儒略日
    DATE_TYPE m_emType;  /// 日期的时区
    int       m_nDecimal;/// 秒的位数
};
}
using namespace Aerospace;
ALGORITHM_EXPORT std::ostream& operator<< (std::ostream& os, const CDate& D);
#endif // YTY_DATE_H
