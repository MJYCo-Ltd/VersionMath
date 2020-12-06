// YTY_IRESInfo.cpp: implementation of the CVVP_IRESInfo class.
//
//////////////////////////////////////////////////////////////////////
#include <fstream>
#include <cstring>
using namespace std;
#include <Satellite/IRESInfo.h>
#include "sofa.h"

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////
using namespace Aerospace;
/// 拉格朗日插值
template<typename T, int N>
T ItLagrange(double t, const double *x, const T *y)
{
    int i, j;
    T s;
    memset(&s,0,sizeof(s));

    for (i = 0; i < N; ++i)
    {
        if (fabs(t - x[i]) < 1e-8)
            return y[i];
        double X = 1;
        for (j = 0; j < N; ++j)
        {
            if (i != j)
                X = X * (t - x[j]) / (x[i] - x[j]);
        }
        s += X * y[i];
    }
    return s;
}

/// 牛顿插值
template<typename T, int N>
T ItNewton(double t, const double *x, const T *y)
{
    int i, j;
    T S[N + 1], p;

    //复制数据到S
    for (i = 0; i < N; ++i)
        S[i] = y[i];

    for (j = 0; j < N - 1; ++j)
        for (i = N - 1; i > j; --i)
            S[i] = (S[i] - S[i - 1]) / (x[i] - x[i - j - 1]);
        p = S[N - 1];
        for (i = N - 2; i >= 0; --i)
            p = p * (t - x[i]) + S[i];
        return p;
}

CIRESInfo* CIRESInfo::m_pThis = nullptr;
CIRESInfo* CIRESInfo::GetInstance()
{
    if(nullptr== m_pThis)
	{
        m_pThis = new CIRESInfo;
		return m_pThis;
	}
	else
	{
		return m_pThis;
	}
}

void CIRESInfo::Realse()
{
    delete m_pThis;
    m_pThis = nullptr;
}

CIRESInfo::CIRESInfo()
    :m_unCount(0)
    ,m_bInit(false)
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

BulletinB CIRESInfo::GetBulletinB(double dMJD) const
{
	BulletinB tempBulletinB;
    /// 如果数据为空 则返回 空值
	if(m_vMJD.empty())
	{
		memset(&tempBulletinB,0,sizeof(tempBulletinB));
		return tempBulletinB;
	}

    /// 如果小于第一个变量,则返回第一个值
	if(dMJD < m_vMJD[0])
	{
        tempBulletinB.dPX      = m_vPx[0];
        tempBulletinB.dPY      = m_vPy[0];
		tempBulletinB.dUT1_UTC = m_vUt1_Utc[0];
	}
    /// 如果时间大于最后一个变量,则返回最后一个值
	else if(dMJD > m_vMJD[m_unCount-1])
	{
        tempBulletinB.dPX      = m_vPx[m_unCount-1];
        tempBulletinB.dPY      = m_vPy[m_unCount-1];
		tempBulletinB.dUT1_UTC = m_vUt1_Utc[m_unCount-1];
	}
	else
	{
		if (m_unCount == 2)
		{
            tempBulletinB.dPX = ItLagrange<double, 2>(dMJD, &m_vMJD[0],
                                                      &m_vPx[0]);
            tempBulletinB.dPY = ItLagrange<double, 2>(dMJD, &m_vMJD[0],
                                                      &m_vPy[0]);
            tempBulletinB.dUT1_UTC = ItLagrange<double, 2>(dMJD, &m_vMJD[0],
                                                           &m_vUt1_Utc[0]);
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

			tempBulletinB.dPX = ItNewton<double, 3>(dMJD, &m_vMJD[nPos],
				&m_vPx[nPos]);
			tempBulletinB.dPY = ItNewton<double, 3>(dMJD, &m_vMJD[nPos],
				&m_vPy[nPos]);
			tempBulletinB.dUT1_UTC = ItNewton<double, 3>(dMJD, &m_vMJD[nPos],
				&m_vUt1_Utc[nPos]);
		}
        else
        {
            tempBulletinB.dPX = m_vPx[0];
            tempBulletinB.dPY = m_vPy[0];
            tempBulletinB.dUT1_UTC = m_vUt1_Utc[0];
        }
	}

	return tempBulletinB;
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
