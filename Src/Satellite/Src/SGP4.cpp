#include <iomanip>
#include <Math/VecMat.h>
#include <Satellite/Date.h>
#include <Satellite/SGP4.h>

#include "SGP4/sgp4io.h"
#include "sofam.h"
using namespace Satellite;
using namespace Math;

CSGP4::CSGP4(const string &strLine1, const string &strLine2):
    m_bValid(false),m_dEpochMJD(-1.)
{
    m_pTle = new elsetrec;
    SetTLE(strLine1,strLine2);
}

CSGP4::CSGP4(const char strLine1[], const char strLine2[])
{
    m_pTle = new elsetrec;
    if(nullptr == strLine1 || nullptr == strLine2)
    {
        memset(m_strTLE[0],0,TLELENGTH);
        memset(m_strTLE[1],0,TLELENGTH);
        m_dEpochMJD = -1;
        m_Oribit.Resize(0);

        m_bValid = false;
        return;
    }
    SetTLE(strLine1,strLine2);
}

CSGP4::~CSGP4()
{
    delete m_pTle;
    m_pTle = nullptr;
}

/// 设置TLE数据
///
bool CSGP4::SetTLE(const string &strLine1, const string &strLine2)
{
    /// 判断数据是否非法
    if(strLine1.size()<FIRSTVALUELENGTH || strLine2.size()<SECONDVALUELENGTH)
    {
        memset(m_strTLE[0],0,TLELENGTH);
        memset(m_strTLE[1],0,TLELENGTH);
        m_dEpochMJD = -1;
        m_Oribit.Resize(0);

        m_bValid = false;
        return(m_bValid);
    }
    else
    {
        strncpy(m_strTLE[0],strLine1.c_str(),TLELENGTH);
        strncpy(m_strTLE[1],strLine2.c_str(),TLELENGTH);

        /// 解析tle数据信息
        twoline2rv(m_strTLE[0],m_strTLE[1],wgs72,*m_pTle);
        m_bValid = m_pTle != nullptr;

        if(!m_bValid)
        {
            return(m_bValid);
        }
    }

    if(m_Oribit.Size() < 6)
    {
        m_Oribit.Resize(6);
    }

    double mu, radiusearthkm, tumin, xke, j2, j3, j4, j3oj2;

    /// 获取wgs84的相关系数
    getgravconst( wgs72, tumin, mu, radiusearthkm, xke, j2, j3, j4, j3oj2 );

    /// 计算相关系数
    m_Oribit(0) = m_pTle->a*radiusearthkm*1e3;
    m_Oribit(1) = m_pTle->ecco;
    m_Oribit(2) = m_pTle->inclo;
    m_Oribit(3) = m_pTle->nodeo;
    m_Oribit(4) = m_pTle->argpo;
    m_Oribit(5) = m_pTle->mo;

    /// 设置约简儒略日时间
    m_dEpochMJD = m_pTle->jdsatepoch - DJM0;

    return(m_bValid);
}

/// 计算位置、速度
CVector CSGP4::CalPV(const double &dMJD)
{
    CVector tempVector;

    if(!m_bValid)
    {
        return(tempVector);
    }

    double dMinutes = (dMJD + DJM0 - m_pTle->jdsatepoch)*1440.;

    double dR[3],dV[3];
    if(!sgp4(wgs72,*m_pTle,dMinutes,dR,dV))
    {
        return(tempVector);
    }

    tempVector.Resize(6);
    const static double dKM = 1.0e3;
    tempVector(0) = dR[0]*dKM;
    tempVector(1) = dR[1]*dKM;
    tempVector(2) = dR[2]*dKM;
    tempVector(3) = dV[0]*dKM;
    tempVector(4) = dV[1]*dKM;
    tempVector(5) = dV[2]*dKM;

    return(tempVector);
}

/// 通过位置速度计算经典根数(开普勒六根数)
CVector CSGP4::ClassicalElements(const CVector &vRV)
{
    CVector tmpVector;
    if(vRV.Size() < 6)
    {
        return tmpVector;
    }
    double rv[6];

    const double dRad = 1e-3;
    for(int i=0; i<6; ++i)
    {
        rv[i] = vRV(i)*dRad;
    }

    double mu, radiusearthkm, tumin, xke, j2, j3, j4, j3oj2;

    getgravconst( wgs72, tumin, mu, radiusearthkm, xke, j2, j3, j4, j3oj2 );

    tmpVector.Resize(6);

    double p,nu,d1,d2,d3;
    rv2coe(rv,rv+3,mu,p,tmpVector(0),tmpVector(1),tmpVector(2),tmpVector(3),
           tmpVector(4),nu,tmpVector(5),d1,d2,d3);

    tmpVector(0) *= 1e3;
    return (tmpVector);
}
