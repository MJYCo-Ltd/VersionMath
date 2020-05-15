#include <iostream>
#include "IRESInfo.h"
#include "Date.h"
#include "TimeSys.h"
#include "sofa.h"
#include "CommonAlgorithm.h"

/*****************************************
  各时统的转换关系
  TT  = TAI + 32.184s
  GPS = TAI - 19s
  TAI = UTC + Ns     (N 为闰秒数，从IERS公报C获取)
  UTC = UT1 + ns     (n 大于0,小于1,从IERS公报B获取)
  TCG = TT + LG×(JD-2443144.5)×86400s
 (
  LG = 6.9692903×10e-10

                《卫星轨道模型方法和应用》
  (德)门斯布吕克(Oliver Montenbruck),(德)吉尔(Eberhard Gill)著
                    王家松        祝开建        胡小工      译
  国防工业出版社,2012.4第1版第1次印刷
  156页
 )
  TCB = TCG + LC×(JD-2443144.5)×86400s + P
 (
  LC = 1.4808268457×10e-8
  P ≈ 0.0016568×sin(35999°.37×T+357°.5)
    + 0.0000224×sin(32964°.5×T+246°)
    + 0.0000138×sin(71998°.7×T+355°)
    + 0.0000048×sin(3034°.9×T+25°)
    + 0.0000047×sin(34777°.3×T+230°)
  T = (JD-2451545.0)/36525
                《卫星轨道模型方法和应用》
  (德)门斯布吕克(Oliver Montenbruck),(德)吉尔(Eberhard Gill)著
                    王家松        祝开建        胡小工      译
  国防工业出版社,2012.4第1版第1次印刷
  157页
 )
  TDB = TCB - LB×(JD-2443144.5)×86400s
 (
  LB = LC + LG = 1.5505197487×10e-8
                《卫星轨道模型方法和应用》
  (德)门斯布吕克(Oliver Montenbruck),(德)吉尔(Eberhard Gill)著
                    王家松        祝开建        胡小工      译
  国防工业出版社,2012.4第1版第1次印刷
  158页
 )
 *****************************************/
using namespace Aerospace;
static const double TAI_GPS(19./DAYSEC);
static const double ELC(1.480868457e-8);

CTimeSys::CTimeSys():m_dateUtc(0.0)
{
}

CTimeSys::CTimeSys(double dMjdUTC):m_dateUtc(dMjdUTC)
{
}

CTimeSys::CTimeSys(const CDate &rDate)
{
    m_dateUtc = rDate;

    /// 强制设置时间为UTC时间
    if(UTC != rDate.GetDateType())
    {
        m_dateUtc.SetDateType(UTC);
    }
}

double CTimeSys::GetTT() const
{
    double dTai = GetTAI();
    if(0 == dTai)
    {
        return 0.0;
    }
    return (dTai+TTMTAI/DAYSEC);
}

double CTimeSys::GetTAI() const
{
    if(!m_dateUtc)
    {
        return(0);
    }

    int nYear, nMonth, nDay, nHour, nMinute;
    double dF,dSecond,dTAI_UTC;

    /// 获取年月日
    if(!m_dateUtc.GetDate(nYear,nMonth,nDay,nHour,nMinute,dSecond))
    {
        return(0);
    }

    iauTf2d('+',nHour,nMinute,dSecond,&dF);

    /// 判断输入是否合法
    if(iauDat(nYear,nMonth,nDay,dF,&dTAI_UTC)<0)
    {
        return (0);
    }

    /// 返回TAI时间
    return (m_dateUtc.GetMJD()+dTAI_UTC/DAYSEC);
}

double CTimeSys::GetUT1() const
{
    if(!m_dateUtc)
    {
        return(0);
    }

    double dUtc = m_dateUtc.GetMJD();

    /// 获取 读取IERS 文件的类的实例
    CIRESInfo* pIERS = CIRESInfo::GetInstance();

    /// 判断是否初始化
    if(!pIERS || !pIERS->IsInit())
    {
        std::cerr<<" ERROR IERS file Read faild!"<<std::endl;
        return(0);
    }

    /// 获取数据
    BulletinB rB = pIERS->GetBulletinB(m_dateUtc.GetMJD());

    return(dUtc+rB.dUT1_UTC/DAYSEC);
}

double CTimeSys::GetTDB() const
{
    double dTcb = GetTCB();
    if(0 == dTcb)
    {
        return(0);
    }

    double dJd1,dJd2;
    /// 将TCB转换成TDB
    iauTcbtdb(DJM0,dTcb,&dJd1,&dJd2);

    return(dJd1+dJd2-DJM0);
}

double CTimeSys::GetTCB() const
{
    double dTcg = GetTCG();
    if(0 == dTcg)
    {
        return(0);
    }

    static const double t77t = DJM77 + TTMTAI/DAYSEC;

    double dT = (dTcg - DJM00)/DJC;
    double dTcb = dTcg +(dTcg-t77t)*ELC+CalP(dT)/DAYSEC;

    return(dTcb);
}

double CTimeSys::GetTCG() const
{
    double dJd1,dJd2,dTT;
    dTT = GetTT();
    if(0 == dTT)
    {
        return(0);
    }

    iauTttcg(DJM0,dTT,&dJd1,&dJd2);
    /// 返回TCG
    return(dJd1+dJd2-DJM0);
}

double CTimeSys::GetGPS() const
{
    double dTai = GetTAI();
    if(0 == dTai)
    {
        return(0);
    }

    return (dTai- TAI_GPS);
}
