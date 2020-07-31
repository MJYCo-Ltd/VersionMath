#include <iomanip>
using namespace std;
#include "Date.h"
#include "sofa.h"

using namespace Aerospace;
static const double BJ_UTC(1./3.);
static const char SCUTC[] = "UTC";

/// 此处将m_dMjd放置在前面，防止警告
CDate::CDate(DATE_TYPE emType):m_dMjd(0),m_emType(emType)
{
    m_nDecimal = 3;
}

CDate::CDate(int nYear, int nMonth, int nDay
                   , int nHour, int nMinute, double dSecond
                   , DATE_TYPE emType):m_emType(emType)
{
    double d1,d2;
    int nStatus;
    nStatus = iauDtf2d(SCUTC,nYear,nMonth,nDay,nHour,nMinute,dSecond,&d1,&d2);
    if(nStatus>=0)
    {
        m_dMjd = d1 + d2 - DJM0;
    }
    else
    {
        m_dMjd = 0;
    }

    /// 北京时间较UTC时间早8个小时
    if(BJ == m_emType && 0 != m_dMjd)
    {
        m_dMjd -= BJ_UTC;
    }

    m_nDecimal = 3;
}

CDate::CDate(double dMjd,DATE_TYPE emType):m_dMjd(dMjd),m_emType(emType)
{
    m_nDecimal = 3;
}

CDate::CDate(const CDate &rDate)
{
    this->m_dMjd = rDate.m_dMjd;
    this->m_emType = rDate.m_emType;
    SetSecondDecimal(rDate.m_nDecimal);
}

void CDate::SetDateType(const DATE_TYPE& emType)
{
    if(emType == m_emType)
    {
        return;
    }
    else
    {
        m_emType = emType;
    }
}

/// 获取日历信息
bool CDate::GetDate(int &nYear, int &nMonth, int &nDay
                      , int &nHour, int &nMinute, double &dSecond) const
{
    /// 查看数据是否有效
    if(!this->operator bool())
    {
        return(false);
    }

    double dMjd = m_dMjd;
    if(BJ == m_emType)
    {
        dMjd += BJ_UTC;
    }

    int hms[4];
    iauD2dtf(SCUTC,m_nDecimal,DJM0,dMjd,&nYear,&nMonth,&nDay,hms);

    double d = pow(.1,m_nDecimal);
    /// 设置时分秒
    nHour = hms[0];
    nMinute = hms[1];
    dSecond = hms[2] + double(hms[3])*d;

    return(true);
}

void CDate::SetSecondDecimal(int nDecimal)
{
    if(nDecimal < 1 || nDecimal > 9)
    {
        return;
    }
    else
    {
        m_nDecimal = nDecimal;
    }
}

void CDate::SetJD(double dJD, DATE_TYPE emType)
{
    m_dMjd = dJD - DJM0;
    if(BJ == emType)
    {
        m_dMjd -= BJ_UTC;
    }
}

double CDate::GetJD() const
{
    return(m_dMjd+DJM0);
}

void CDate::operator +=(double dTime)
{
    m_dMjd += dTime;
}

void CDate::operator -=(double dTime)
{
    m_dMjd -= dTime;
}

bool CDate::operator ==(const CDate& rDate)
{
    /// 判断地址是否相同
    if(this == &rDate)
    {
        return(true);
    }

    /// 判断约简儒略日是否相同
    if(fabs(m_dMjd - rDate.m_dMjd) > 1e-16)
    {
        return(false);
    }

    /// 判断时间类型是否相同
    if(m_emType != rDate.m_emType)
    {
        return(false);
    }

    return(true);
}

bool CDate::operator !=(const CDate& rDate)
{
    return !(this->operator ==(rDate));
}

CDate& CDate::operator =(const double dMjd)
{
    m_dMjd = dMjd;
    return(*this);
}

CDate& CDate::operator =(const CDate& rDate)
{
    if(this == &rDate)
    {
        return(*this);
    }
    else
    {
        m_dMjd = rDate.m_dMjd;
        m_emType = rDate.m_emType;
        SetSecondDecimal(rDate.m_nDecimal);
    }

    return(*this);
}

std::ostream& operator << (std::ostream& os, const CDate& rDate)
{
  // Variables
  int    Year, Month, Day;
  int    H, M;
  double S;
  if(!rDate)
  {
      os << "Invalid Date";
      return os;
  }
  else
  {
      rDate.GetDate(Year,Month,Day,H,M,S);
  }

  DATE_TYPE emType = rDate.GetDateType();
  switch(emType)
  {
  case UTC:
      os << "UTC ";
      break;
  case BJ:
      os << "BJ  ";
      break;
  default:
      os << "UNKNOWN";
  }

  os << setfill('0')
     << setw(4) << Year  << "/"
     << setw(2) << Month << "/"
     << setw(2) << Day   << "  ";
  os << setw(2) << H << ":"
     << setw(2) << M << ":"
     << fixed << setprecision(rDate.GetSecondDecimal())
     << setw(6) << S
     << setfill(' ');

  return os;
}
