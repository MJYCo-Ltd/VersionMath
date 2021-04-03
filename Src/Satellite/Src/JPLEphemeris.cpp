#include <string>
using namespace std;
#include <Math/VecMat.h>
#include <Satellite/JPLEphemeris.h>
#include "jpl_eph/jpleph.h"
#include "jpl_eph/jpl_int.h"
#include "sofam.h"
using namespace Aerospace;
using namespace Math;

CJPLEphemeris* CJPLEphemeris::m_pInstance = 0;

CJPLEphemeris* CJPLEphemeris::GetInstance()
{
    if(0 == m_pInstance)
    {
        m_pInstance = new CJPLEphemeris;
    }

    return(m_pInstance);
}

void CJPLEphemeris::Realse()
{
    delete(m_pInstance);
    m_pInstance=nullptr;
}

CJPLEphemeris::CJPLEphemeris():m_pJplDate(0),m_bInit(false),m_dAu(0)
{
}

CJPLEphemeris::~CJPLEphemeris()
{
    /// 关闭数据
    if(0 != m_pJplDate)
    {
        jpl_close_ephemeris(m_pJplDate);
    }

    m_pJplDate = 0;
    m_dAu = 0;
}

bool CJPLEphemeris::Init(const string &sFileName)
{
    /// 如果数据不为空，则关闭星历表
    if(0 != m_pJplDate)
    {
        jpl_close_ephemeris(m_pJplDate);
    }

    /// 将数据置空
    m_pJplDate = 0;

    /// 初始化星历表
    m_pJplDate = jpl_init_ephemeris(sFileName.c_str(),0,0);

    if(0 == m_pJplDate)
    {
        m_bInit = false;
    }
    else
    {
        m_dAu = jpl_get_double(m_pJplDate,JPL_EPHEM_AU_IN_KM) * 1000.;
        m_bInit = true;
    }

    return(m_bInit);
}

CVector CJPLEphemeris::GetSunPos(const double &dMjdTT) const
{
    return(GetPos(dMjdTT,Sun,Earth));
}

CVector CJPLEphemeris::GetMoonPos(const double &dMjdTT) const
{
    return(GetPos(dMjdTT,Moon,Earth));
}

CVector CJPLEphemeris::GetPos(const double &dMjdTT, const PLANET_TYPE& planet1
                                                          , const PLANET_TYPE &centerPlanet) const
{
    CVector tmpVet;

    if(!m_bInit)
    {
        return tmpVet;
    }

    tmpVet.Resize(3);
    if(planet1 == centerPlanet)
    {
        return(tmpVet);
    }

    double tmpRRD[6];
    /// 最后一个 如果为 0 则只计算位置； 非0 则计算速度
    int nError = jpl_pleph(m_pJplDate,dMjdTT+DJM0,planet1,centerPlanet,tmpRRD,0);
    if(0 != nError)
    {
        return (tmpVet);
    }

    /// 转换成标准单位
    tmpVet(0) = tmpRRD[0]*m_dAu;
    tmpVet(1) = tmpRRD[1]*m_dAu;
    tmpVet(2) = tmpRRD[2]*m_dAu;
    return (tmpVet);
}
