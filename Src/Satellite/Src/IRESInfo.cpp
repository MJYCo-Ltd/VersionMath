// YTY_IRESInfo.cpp: implementation of the CVVP_IRESInfo class.
//
//////////////////////////////////////////////////////////////////////
#include <fstream>
#include <cstring>
using namespace std;
#include <Satellite/IRESInfo.h>
#include <Math/Intpol.h>
#include "sofa.h"
#include "sofam.h"

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////
using namespace Aerospace;
using namespace Numerical;

CIRESInfo* CIRESInfo::GetInstance()
{
    static CIRESInfo s_IresInfo;

    return(&s_IresInfo);
}

CIRESInfo::CIRESInfo()
{
}

CIRESInfo::~CIRESInfo()
{
}

bool CIRESInfo::Init(const string& strFileName)
{
    /// 判断文件是否为空
    if(strFileName.empty())
    {
        return(false);
    }
    else
    {
        ifstream fileIn;
        fileIn.open(strFileName.c_str());

        /// 查看文件是否存在
        if(!fileIn.is_open())
        {
            return(false);
        }

        char strBuffer[188];

        /// 临时变量用于保存读取文件
        double dMjd,dPx,dPy,dUt1_Utc;
        int nReadSize=0;

        /// 清空数据
        ClearData();
        m_bInit = false;

        while(!fileIn.eof())
        {
            fileIn.getline(strBuffer,188);

            if(strlen(strBuffer) < 187)
            {
                continue;
            }

            /// 根据 readme.finals2000A 文件 读取文件
            nReadSize = sscanf(strBuffer+7,"%lf",&dMjd);
            nReadSize = sscanf(strBuffer+134,"%10lf%10lf%11lf"
                               ,&dPx,&dPy,&dUt1_Utc);
            if(3 == nReadSize)
            {
                m_vMJD.push_back(dMjd);
                m_vPx.push_back(dPx);
                m_vPy.push_back(dPy);
                m_vUt1_Utc.push_back(dUt1_Utc);
                ++m_unCount;
            }
            else
            {
                nReadSize = 0;
                nReadSize += sscanf(strBuffer+18,"%9lf",&dPx);
                nReadSize += sscanf(strBuffer+37,"%9lf",&dPy);
                nReadSize += sscanf(strBuffer+58,"%9lf",&dUt1_Utc);
                if(3 == nReadSize)
                {
                    m_vMJD.push_back(dMjd);
                    m_vPx.push_back(dPx);
                    m_vPy.push_back(dPy);
                    m_vUt1_Utc.push_back(dUt1_Utc);
                    ++m_unCount;
                }
            }
        }

        if(0 != m_unCount)
        {
            m_bInit = true;
        }
        /// 关闭文件
        fileIn.close();
        return(m_bInit);
    }
}

bool CIRESInfo::GetBulletinB(double dMJD,BulletinB& bulletinB) const
{
    /// 如果数据为空 则返回 空值
    if(0 == m_unCount || dMJD < m_vMJD[0] || dMJD > m_vMJD[m_unCount-1])
    {
        return false;
    }


    if (m_unCount == 2)
    {
        bulletinB.dPX = Intpol::ItLagrange(2,&m_vMJD[0],
                &m_vPx[0],dMJD);
        bulletinB.dPY = Intpol::ItLagrange(2,&m_vMJD[0],
                &m_vPy[0],dMJD);
        bulletinB.dUT1_UTC = Intpol::ItLagrange(2,&m_vMJD[0],
                &m_vUt1_Utc[0],dMJD);
    }
    else if(m_unCount >= 3)
    {
        int nPos = (int)(dMJD-m_vMJD[0]);

        /// 防止向上访问越界
        if (nPos > int(m_unCount - 3))
        {
            nPos=int(m_unCount-3);
        }

        /// 防止向下访问越界
        if (nPos < 0)
        {
            nPos = 0;
        }

        bulletinB.dPX = Intpol::ItNewton(3,&m_vMJD[nPos],
                                        &m_vPx[nPos],dMJD);
        bulletinB.dPY = Intpol::ItNewton(3,&m_vMJD[nPos],
                                        &m_vPy[nPos],dMJD);
        bulletinB.dUT1_UTC = Intpol::ItNewton(3,&m_vMJD[nPos],
                                              &m_vUt1_Utc[nPos],dMJD);
    }

    return(true);
}


double CIRESInfo::GetUT2_UT1(double dMJD) const
{
    double dTemp = iauEpb(DJM0,dMJD);
    double dUt2_Ut1 = 0.022*sin(2*DPI*dTemp) - 0.012*cos(2*DPI*dTemp) - 0.006*sin(4*DPI*dTemp) + 0.007*cos(4*DPI*dTemp);
    return(dUt2_Ut1);
}

void CIRESInfo::ClearData()
{
    m_vMJD.clear();
    m_vPx.clear();
    m_vPy.clear();
    m_vUt1_Utc.clear();
    m_unCount = 0;
}
