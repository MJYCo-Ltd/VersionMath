#include <limits>
#include <cstring>
#include <iostream>
#include <cmath>
#include <VersionMathCommon.h>
#include <Math/Vector.h>
#include <Math/MemPool.h>
using namespace std;
using namespace Math;

const CVector CVector::X_AXIS(1,0,0,true);
const CVector CVector::Y_AXIS(0,1,0,true);
const CVector CVector::Z_AXIS(0,0,1,true);
const CVector CVector::NULL_VECTOR;

CVector::CVector()
{
}

/// 构造一个全为零的向量
CVector::CVector (unsigned int nSize)
{
    m_unDim = nSize;
    /// 开辟空间
    m_pdV = CMemPool::GetInstance()->Create<double>(m_unDim);
}

/// 复制构造
CVector::CVector (const CVector& rV)
{
    if(rV.m_unDim < 1) return;

    m_unDim = rV.m_unDim;
    m_pdV = CMemPool::GetInstance()->Create<double>(m_unDim);

    /// 复制信息
    memcpy(m_pdV,rV.m_pdV,sizeof(*m_pdV)*m_unDim);
}

/// 通过数组构造
CVector::CVector (const double* pdArry, unsigned int nDim)
{
    if(nDim < 1 || nullptr == pdArry) return;

    m_unDim = nDim;
    m_pdV = CMemPool::GetInstance()->Create<double>(m_unDim);

    /// 复制信息
    memcpy(m_pdV,pdArry,sizeof(*m_pdV)*nDim);
}

/// 3维向量
CVector::CVector (double dX, double dY, double dZ, bool bIsNormal)
    : m_unDim(3),m_bIsNormal(bIsNormal)
{
    m_pdV = CMemPool::GetInstance()->Create<double>(m_unDim);

    m_pdV[0]=dX;
    m_pdV[1]=dY;
    m_pdV[2]=dZ;
}

/// 6维向量
CVector::CVector (double dx, double dy, double dz
                          ,double dX, double dY, double dZ)
    : m_unDim(6)
{
    m_pdV = CMemPool::GetInstance()->Create<double>(m_unDim);

    m_pdV[0]=dx;
    m_pdV[1]=dy;
    m_pdV[2]=dz;

    m_pdV[3]=dX;
    m_pdV[4]=dY;
    m_pdV[5]=dZ;
}

CVector::~CVector()
{
    if(nullptr != m_pdV)
    {
        CMemPool::GetInstance()->Remove(m_pdV);
        m_pdV = nullptr;
    }

    m_unDim = 0;
}

/// Set方法
void CVector::Set(double xT, double yT, double zT)
{
    if(m_unDim<3)
    {
        Resize(3);
    }

    m_pdV[0] = xT;
    m_pdV[1] = yT;
    m_pdV[2] = zT;
}

/// 重设大小
CVector& CVector::Resize(unsigned int unSize)
{
    /// 如果相等则直接返回
    if (m_unDim==unSize) return (*this);

    if(unSize < 1)
    {
        m_unDim = 0;
        CMemPool::GetInstance()->Remove(m_pdV);
        return(*this);
    }

    /// 开辟新的空间
    double *pNew = CMemPool::GetInstance()->Create<double>(unSize);;

    unsigned int unMax = ((unSize<m_unDim)? unSize : m_unDim);

    if(unMax > 0)
    {
        memcpy(pNew,m_pdV,unMax*sizeof(*pNew));
    }

    /// 释放旧空间
    CMemPool::GetInstance()->Remove(m_pdV);

    /// 内存指向新的地址
    m_pdV = pNew;
    m_unDim = unSize;

    return (*this);
}

/// 计算向量长度
double CVector::Length() const
{
    double target = 0.0;
    for(unsigned int i = 0; i < m_unDim; ++i)
    {
        target += m_pdV[i] * m_pdV[i];
    }

    return sqrt(target);
}

/// 将向量归一化
void CVector::Normalize()
{
    static const double epsilon(0.0000001);
    double r = Length();
    if(r==1)
    {
        m_bIsNormal = true;
        return;
    }

    if(r > epsilon)
    {

        this->operator /=(r);
        m_bIsNormal = true;
    }
}

/// 访问部分数据
CVector CVector::slice (unsigned int nFirst,unsigned int nLast) const
{
    /// 设置
    unsigned int i
       ,nMax= nLast < m_unDim ? nLast : m_unDim;

    CVector Aux(nLast-nFirst+1);

    /// 填充
    for (i=nFirst; i<=nMax; ++i)
    {
        Aux.m_pdV[i-nFirst]=m_pdV[i];
    }

    return Aux;
}


/// 获取向量的开平方向量
CVector CVector::Sqrt() const
{
    CVector Aux(m_unDim);
    for (unsigned int i=0; i<m_unDim; ++i)
    {
        Aux.m_pdV[i]=sqrt(m_pdV[i]);
    }
    return Aux;
}


/// 赋值
CVector& CVector::operator=(const double value)
{
    for (unsigned int i=0; i<m_unDim; ++i)
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

    Resize(V.m_unDim);

    /// 复制数据
    memcpy(m_pdV,V.m_pdV,sizeof(*m_pdV)*m_unDim);

    return (*this);
}


// Vector addition/subtraction with assignment

void CVector::operator += (const CVector& V)
{
    /// 如果维数不相等则不进行操作
    if (m_unDim!=V.m_unDim)
    {
        cerr << "ERROR: Incompatible shape in Vector operator+=(Vector)" << endl;
        return;
    };

    for (unsigned int i=0; i<m_unDim; ++i)
    {
        m_pdV[i]+=V.m_pdV[i];
    }
}

void CVector::operator -= (const CVector& V)
{
    if (m_unDim!=V.m_unDim)
    {
        cerr << "ERROR: Incompatible shape in Vector operator-=(Vector)" << endl;
        return;
    };

    for (unsigned int i=0; i<m_unDim; ++i)
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
        if(m_unDim != V.m_unDim)
        {
            return(false);
        }
        else
        {
            /// 判断各值是否相同
            for(unsigned int i=0; i<m_unDim; ++i)
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
    for(unsigned int i = 0; i < m_unDim; ++i)
    {
        m_pdV[i] *= value;
    }
    return (*this);
}
