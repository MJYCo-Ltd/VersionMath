#include <limits>
#include <cstring>
#include <iostream>
#include <cmath>
#include "Vector.h"
using namespace std;
using namespace Math;

CVector::CVector():m_nDim(0),m_pdV(0)
{
}

/// 构造一个全为零的向量
CVector::CVector (int nSize)
    : m_nDim(nSize)
{
    /// 开辟空间
    m_pdV = new double[nSize]();
}

/// 复制构造
CVector::CVector (const CVector& rV)
    : m_nDim(rV.m_nDim)
{
    m_pdV = new double [rV.m_nDim];

    /// 复制信息
    memcpy(m_pdV,rV.m_pdV,sizeof(*m_pdV)*rV.m_nDim);
}

/// 通过数组构造
CVector::CVector (const double* pdArry, int nDim)
    : m_nDim(nDim)
{
    m_pdV = new double [m_nDim];

    /// 复制信息
    memcpy(m_pdV,pdArry,sizeof(*m_pdV)*nDim);
}

/// 3维向量
CVector::CVector (double dX, double dY, double dZ)
    : m_nDim(3)
{
    m_pdV = new double[m_nDim];

    m_pdV[0]=dX;
    m_pdV[1]=dY;
    m_pdV[2]=dZ;
}

/// 6维向量
CVector::CVector (double dx, double dy, double dz
                          ,double dX, double dY, double dZ)
    : m_nDim(6)
{
    m_pdV = new double [m_nDim];

    m_pdV[0]=dx;
    m_pdV[1]=dy;
    m_pdV[2]=dz;

    m_pdV[3]=dX;
    m_pdV[4]=dY;
    m_pdV[5]=dZ;
}

CVector::~CVector()
{
    if(0 != m_pdV)
    {
        delete [] m_pdV;
    }
}

/// Set方法
void CVector::Set(double xT, double yT, double zT)
{
    if(m_nDim<3)
    {
        Resize(3);
    }

    m_pdV[0] = xT;
    m_pdV[1] = yT;
    m_pdV[2] = zT;
}

/// 重设大小
CVector& CVector::Resize(int nSize)
{
    /// 如果相等则直接返回
    if (m_nDim==nSize) return (*this);

    int    i       /// 游标
          ,i_max;  /// 需要复制的值的个数

    /// 开辟新的空间
    double *v_new = new double[nSize];

    i_max = ((nSize<m_nDim)? nSize : m_nDim);

    /// 复制旧值
    for (i=0;i<i_max;++i)
    {
        v_new[i]=m_pdV[i];
    }

    /// 填充 0
    for (i=i_max;i<nSize;++i)
    {
        v_new[i]=0.0;
    }

    /// 释放旧空间
    delete [] m_pdV;

    /// 内存指向新的地址
    m_pdV = v_new;
    m_nDim = nSize;

    return (*this);
}

/// 计算向量长度
double CVector::Length() const
{
    double target = 0.0;
    for(int i = 0; i < m_nDim; ++i)
    {
        target += m_pdV[i] * m_pdV[i];
    }

    return sqrt(target);
}

/// 将向量归一化
void CVector::Normalize()
{
    double r = Length();

    this->operator /=(r);
}

/// 访问部分数据
CVector CVector::slice (int nFirst, int nLast) const
{
    /// 设置
    int i
       ,nMax= nLast < m_nDim ? nLast : m_nDim;

    CVector Aux(nLast-nFirst+1);

    /// 填充
    for (i=nFirst; i<=nMax; ++i)
    {
        Aux.m_pdV[i-nFirst]=m_pdV[i];
    }

    /// 填充 防止越界访问
    for(i=nMax+1; i<=nLast; ++i)
    {
        Aux.m_pdV[i-nFirst] = 0;
    }

    return Aux;
}


/// 获取向量的开平方向量
CVector CVector::Sqrt() const
{
    CVector Aux(m_nDim);
    for (int i=0; i<m_nDim; ++i)
    {
        Aux.m_pdV[i]=sqrt(m_pdV[i]);
    }
    return Aux;
}


/// 赋值
CVector& CVector::operator=(const double value)
{
    for (int i=0; i<m_nDim; ++i)
    {
        m_pdV[i]=value;
    }

    return (*this);
}

/// 赋值
CVector& CVector::operator=(const CVector& V)
{
    if (this == &V)
    {
        return (*this);
    }

    /// 检查向量是否为空
    if (IsEmpty())
    {
        m_nDim = V.m_nDim;
        m_pdV = new double [V.m_nDim];
    };

    /// 检查维度是否相同
    if (m_nDim != V.m_nDim)
    {
        cerr << "ERROR: Incompatible sizes in Vector operator=(Vector)" << endl;
        return (*this);
    };

    /// 复制数据
    memcpy(m_pdV,V.m_pdV,sizeof(*m_pdV)*m_nDim);

    return (*this);
}


// Vector addition/subtraction with assignment

void CVector::operator += (const CVector& V)
{
    /// 如果维数不相等则不进行操作
    if (m_nDim!=V.m_nDim)
    {
        cerr << "ERROR: Incompatible shape in Vector operator+=(Vector)" << endl;
        return;
    };

    for (int i=0; i<m_nDim; ++i)
    {
        m_pdV[i]+=V.m_pdV[i];
    }
}

void CVector::operator -= (const CVector& V)
{
    if (m_nDim!=V.m_nDim)
    {
        cerr << "ERROR: Incompatible shape in Vector operator-=(Vector)" << endl;
        return;
    };

    for (int i=0; i<m_nDim; ++i)
    {
        m_pdV[i]-=V.m_pdV[i];
    }
}

bool CVector::operator ==(const CVector& V)const
{
    /// 判断地址是否相同
    if(&V == this)
    {
        return (true);
    }
    else
    {
        /// 判断维数是否相同
        if(m_nDim != V.m_nDim)
        {
            return(false);
        }
        else
        {
            /// 判断各值是否相同
            for(int i=0; i<m_nDim; ++i)
            {
                if(fabs(m_pdV[i]-V.m_pdV[i]) > numeric_limits<double>::min())
                {
                    return(false);
                }
            }

            return(true);
        }
    }
}

bool CVector::operator != (const CVector& V)const
{
    return !(this->operator ==(V));
}


CVector& CVector::operator /=(const double value)
{
    double dInv = 1.0;

    /// 防止除零
    if(0 != value)
    {
        dInv = 1.0 / value;
    }

    /// 调用乘法运算进行计算
    this->operator *=(dInv);

    return (*this);
}

CVector& CVector::operator *=(const double value)
{
    for(int i = 0; i < m_nDim; ++i)
    {
        m_pdV[i] *= value;
    }
    return (*this);
}
