#include <iostream>
#include <cstring>
#include <limits>
#include <cmath>
#include <Math/Matrix.h>
#include <Math/Vector.h>
#include <Math/MemPool.h>
using namespace std;
using namespace Math;
static double sBuffer[9]={1,0,0,0,1,0,0,0,1};
const CMatrix CMatrix::IDENTITY_MAT(sBuffer,3,3);
const CMatrix CMatrix::NULL_MATRIX;

/// 默认构造函数
CMatrix::CMatrix()
{
    InitAttr();
}

CMatrix::CMatrix (unsigned int dim1, unsigned int dim2)
    : m_unRow(dim1), m_unCol(dim2)
{
    /// 防止索引值的出现
    if(m_unRow < 1 || m_unCol < 1)
    {
        cerr << "ERROR: Incompatible shape in Matrix Construction" << endl;
        InitAttr();
        return;
    }

    /// 开辟空间
    m_ppdM = CMemPool::GetInstance()->Create<double*>(m_unRow);
    for (unsigned int i=0; i<m_unRow; ++i)
    {
        m_ppdM[i] = CMemPool::GetInstance()->Create<double>(m_unCol);
    }
}

/// 复制构造函数
CMatrix::CMatrix (const CMatrix& rM)
{
    /// 如果要复制的对象为空矩阵，则直接将矩阵置空
    if(rM.IsEmpty())
    {
        InitAttr();
        return;
    }

    unsigned int i,j;
    m_unRow = rM.m_unRow;
    m_unCol = rM.m_unCol;
    /// 开辟空间
    m_ppdM = CMemPool::GetInstance()->Create<double*>(m_unRow);
    for (i=0; i<m_unRow; ++i)
    {
        m_ppdM[i] = CMemPool::GetInstance()->Create<double>(m_unCol);
    }

    j = sizeof(**m_ppdM)*m_unCol;
    /// 初始化
    for (i=0; i<m_unRow; ++i)
    {
        /// 一行一行的复制
        memcpy(m_ppdM[i],rM.m_ppdM[i],j);
    }
}

/// 通过数组构造
CMatrix::CMatrix (const double* p, unsigned int dim1, unsigned int dim2)
{
    /// 防止数组指向非法地址，或者索引值错误
    if(0 == p || dim1 < 1 || dim2 < 1)
    {
        InitAttr();
        return;
    }

    unsigned int i,j;
    m_unRow = dim1;
    m_unCol = dim2;
    /// 开辟空间
    m_ppdM = CMemPool::GetInstance()->Create<double*>(m_unRow);
    for (i=0; i<m_unRow; ++i)
    {
        m_ppdM[i] = CMemPool::GetInstance()->Create<double>(m_unCol);
    }
    j = sizeof(**m_ppdM)*m_unCol;
    /// 初始化
    for (i=0; i<m_unRow; ++i)
    {
        memcpy(m_ppdM[i],p+i*m_unCol,j);
    }

}

CMatrix::CMatrix(const double p[3][3], unsigned int dim1, unsigned int dim2)
{
    /// 防止数组指向非法地址，或者索引值错误
    if(0 == p || dim1 < 1 || dim2 < 1)
    {
        InitAttr();
        return;
    }

    unsigned int i,j;
    m_unRow = dim1;
    m_unCol = dim2;
    /// 开辟空间
    m_ppdM = CMemPool::GetInstance()->Create<double*>(m_unRow);
    for (i=0; i<m_unRow; ++i)
    {
        m_ppdM[i] = CMemPool::GetInstance()->Create<double>(m_unCol);
    }
    j = sizeof(**m_ppdM)*m_unCol;
    /// 初始化
    for (i=0; i<m_unRow; ++i)
    {
        memcpy(m_ppdM[i],p[i],j);
    }
}

CMatrix::~CMatrix()
{
    Empty();
}

/// 重置矩阵大小
CMatrix& CMatrix::Resize(unsigned int nRow, unsigned int nCol)
{
    /// 如果行列相等则直接返回
    if (m_unRow==nRow && m_unCol==nCol)
    {
        return (*this);
    }

    /// 判断范围是否合法
    if(nRow < 1 || nCol < 1)
    {
        cerr << "ERROR: Incompatible shape in Matrix Resize(int,int)" << endl;
        /// 清空数据
        this->Empty();
        return (*this);
    }

    unsigned int    i,j,i_max,j_max;
    /// 开辟一个新的空间
    double **M_new = CMemPool::GetInstance()->Create<double*>(nRow);
    for (i=0; i<nRow; ++i)
    {
        M_new[i] = CMemPool::GetInstance()->Create<double>(nCol);
    }

    /// 复制已有的数据，如果大于原来的数据则置零
    i_max = ((nRow<m_unRow)? nRow : m_unRow);
    j_max = ((nCol<m_unCol)? nCol : m_unCol);

    /// 复制原来的数据
    j = j_max * sizeof(**M_new);
    for (i=0;i<i_max; ++i)
    {
        /// 一行一行的复制
        memcpy(M_new[i],m_ppdM[i],j);
    }

    /// 填充下方区域
    j = nCol * sizeof(**M_new);
    for (i=i_max;i<nRow;++i)
    {
        /// 一行一行的置零
        memset(M_new[i],0,j);
    }

    /// 填充右侧区域
    for (i=0;i<i_max;++i)
    {
        for (j=j_max;j<nCol;++j)
        {
            M_new[i][j]=0.0;
        }
    }

    /// 释放原来的空间
    for (i=0; i<m_unRow; ++i)
    {
        CMemPool::GetInstance()->Remove(m_ppdM[i]);
    }
    CMemPool::GetInstance()->Remove(m_ppdM);

    /// 将指针指向新开辟的空间
    m_ppdM = M_new;
    m_unRow = nRow;
    m_unCol = nCol;
    return (*this);
}

bool CMatrix::IsId() const
{
    /// 判断是否为空 是否为方阵
    if(IsEmpty() || !IsSquare())
    {
        return(false);
    }

    /// 判断对角线上的值是否都为1.0
    for(unsigned int i=0; i<m_unRow; ++i)
    {
        for(unsigned int j=0; j<m_unCol; ++j)
        {
            if(i == j)
            {
                /// 对角线上的值 为 1.0
                if(1.0 != m_ppdM[i][j])
                {
                    return(false);
                }
            }
            else
            {
                /// 其他位置的值 为 0.0
                if(0.0 == m_ppdM[i][j])
                {
                    return(false);
                }
            }
        }
    }

    return(true);
}

bool CMatrix::Translate(const CVector &vIn, CVector &vOut) const
{
    if(vIn.Size() > 2 && 3 == m_unCol && 3 == m_unRow)
    {
        unsigned int nCount = vIn.Size() / 3;
        double dvX,dvY,dvZ;

        /// 如果输出大小不足以存储则扩展空间
        if(vOut.Size() < nCount*3)
        {
            vOut.Resize(nCount*3);
        }

        for(unsigned int i=0,j=0;i<nCount;++i,j+=3)
        {
            dvX = vIn(j);
            dvY = vIn(j+1);
            dvZ = vIn(j+2);
            vOut(j)   = m_ppdM[0][0]*dvX + m_ppdM[0][1]*dvY + m_ppdM[0][2]*dvZ;
            vOut(j+1) = m_ppdM[1][0]*dvX + m_ppdM[1][1]*dvY + m_ppdM[1][2]*dvZ;
            vOut(j+2) = m_ppdM[2][0]*dvX + m_ppdM[2][1]*dvY + m_ppdM[2][2]*dvZ;
        }

        return(true);
    }
    else
    {
        return(false);
    }
}

/// 赋值操作
CMatrix& CMatrix::operator=(const double value)
{
    for (unsigned int i=0; i<m_unRow; ++i)
    {
        for (unsigned int j=0; j<m_unCol; ++j)
        {
            m_ppdM[i][j]=value;
        }
    }
    return (*this);
}

CMatrix& CMatrix::operator=(const CMatrix& rM)
{
    if (this == &rM)
    {
        return (*this);
    }

    unsigned int i,j;
    /// 如果矩阵为空则开辟空间
    if (IsEmpty())
    {
        m_unRow = rM.m_unRow;
        m_unCol = rM.m_unCol;
        m_ppdM = CMemPool::GetInstance()->Create<double*>(m_unRow);
        for (i=0; i<m_unRow; ++i)
        {
            m_ppdM[i] = CMemPool::GetInstance()->Create<double>(m_unCol);
        }
    }

    /// 如果矩阵的行列不同则直接返回
    if ( (m_unRow!=rM.m_unRow) || (m_unCol!=rM.m_unCol) )
    {
        cerr << "ERROR: Incompatible shapes in Matrix operator=(Matrix)" << endl;
        return(*this);
    }

    /// 计算每一行的空间大小
    j = m_unCol * sizeof(**m_ppdM);
    for (i=0; i<m_unRow; ++i)
    {
        /// 一行一行的复制
        memcpy(m_ppdM[i],rM.m_ppdM[i],j);
    }
    return (*this);
}


/// 返回指定列数据
CVector CMatrix::GetCol(unsigned int j) const
{
    CVector Res(m_unRow);
    for (unsigned int i=0; i<m_unRow; ++i)
    {
        Res(i)=m_ppdM[i][j];
    }
    return Res;
}

/// 返回指定行数据
CVector CMatrix::GetRow(unsigned int i) const
{
    CVector Res(m_unCol);
    for (unsigned int j=0; j<m_unCol; ++j)
    {
        Res(j)=m_ppdM[i][j];
    }
    return Res;
}

CVector CMatrix::Diag() const
{
    if (IsEmpty() || !IsSquare())
    {
        cerr << "ERROR: Invalid shape in Matrix.Diag()" << endl;
        return(CVector::NULL_VECTOR);
    }

    CVector Vec(m_unRow);
    /// 获取对角线上的值
    for (unsigned int i=0; i<m_unRow; ++i)
    {
        Vec(i) = m_ppdM[i][i];
    }
    return (Vec);
}

double CMatrix::Trace() const
{
    return this->Trace(0,m_unRow-1);
}

double CMatrix::Trace(unsigned int low, unsigned int upp) const
{
    double tmp = 0.0;

    /// 不是方阵则退出
    if (!IsSquare())
    {
        cerr << "ERROR: Invalid shape in Matrix.Trace()" << endl;
        return(tmp);
    }

    /// 判断参数是否合法
    if (m_unRow<=upp)
    {
        cerr << "ERROR: Invalid arguments in Matrix.Trace()" << endl;
        return(tmp);
    }

    for (unsigned int i=low; i<=upp; ++i)
    {
        tmp += m_ppdM[i][i];
    }
    return (tmp);
}

/// 获取矩阵的部分
CMatrix CMatrix::slice(unsigned int first_row, unsigned int last_row, unsigned int first_col, unsigned int last_col)
{
    /// 判断边界值是否满足
    if (last_row<first_row || m_unRow-1<last_row
     || last_col<first_col || m_unCol-1<last_col)
    {
        cerr << "ERROR: Invalid arguments in Matrix.slice()" << endl;
        return(NULL_MATRIX);
    }

    /// 重置矩阵大小
    CMatrix Aux(last_row-first_row+1,last_col-first_col+1);

    /// 获取值
    for (unsigned int i=0;i<=last_row-first_row;++i)
    {
        for(unsigned int j=0;j<=last_col-first_col;++j)
        {
            Aux(i,j) = m_ppdM[i+first_row][j+first_col];
        }
    }

    return (Aux);
}


bool CMatrix::SetCol(unsigned int j, const CVector& Col)
{
    /// 判断向量是否可以放入矩阵
    if (Col.Size()!=m_unRow)
    {
        cerr << "ERROR: Incompatible shapes in Matrix.SetCol()" << endl;
        return(false);
    }

    /// 判断列是否超出范围
    if (m_unCol<=j)
    {
        cerr << "ERROR: Column index out of range in Matrix.SetCol()" << endl;
        return(false);
    }

    for (unsigned int i=0; i<m_unRow; ++i)
    {
        m_ppdM[i][j]=Col(i);
    }

    return(true);
}

bool CMatrix::SetRow(unsigned int i, const CVector& Row)
{
    /// 判断向量是否可以放入矩阵
    if (Row.Size()!=m_unCol)
    {
        cerr << "ERROR: Incompatible shapes in Matrix.SetRow()" << endl;
        return(false);
    }

    /// 判断行是否超出范围
    if (m_unRow<=i)
    {
        cerr << "ERROR: Row index out of range in Matrix.SetRow()" << endl;
        return(false);
    }

    for (unsigned int j=0; j<m_unCol; ++j)
    {
        m_ppdM[i][j]=Row(j);
    }

    return(true);
}



/// 矩阵的加法减法
bool CMatrix::operator +=(const CMatrix& rM)
{
    /// 检查矩阵的行列是否相同
    if ( (m_unRow!=rM.m_unRow) || (m_unCol!=rM.m_unCol) )
    {
        cerr << "ERROR: Incompatible shape in Matrix operator+=(Matrix)" << endl;
        return(false);
    };

    /// 进行加法运算
    for (unsigned int i=0; i<m_unRow; ++i)
    {
        for (unsigned int j=0; j<m_unCol; ++j)
        {
            this->m_ppdM[i][j]+=rM.m_ppdM[i][j];
        }
    }

    return(true);
}

bool CMatrix::operator -=(const CMatrix& rM)
{
    /// 检查矩阵的行列是否相同
    if ( (m_unRow!=rM.m_unRow) || (m_unCol!=rM.m_unCol) )
    {
        cerr << "ERROR: Incompatible shape in Matrix operator-=(Matrix)" << endl;
        return(false);
    };

    /// 进行各值取反
    for (unsigned int i=0; i<m_unRow; ++i)
    {
        for (unsigned int j=0; j<m_unCol; ++j)
        {
            this->m_ppdM[i][j]-=rM.m_ppdM[i][j];
        }
    }

    return(true);
}

bool CMatrix::operator ==(const CMatrix& rM) const
{
    /// 如果地址相等则认为是同一个对象
    if(this == &rM)
    {
        return(true);
    }

    /// 首先判断行列是否相等
    if(m_unRow != rM.m_unRow || m_unCol != rM.m_unCol)
    {
        return(false);
    }

    /// 判断各值是否相同
    for(unsigned int i=0; i<m_unRow; ++i)
    {
        for(unsigned int j=0; j<m_unCol; ++j)
        {
            if(fabs(m_ppdM[i][j]-rM.m_ppdM[i][j]) > numeric_limits<double>::min())
            {
                return(false);
            }
        }
    }

    return (true);
}

bool CMatrix::operator !=(const CMatrix& rM)const
{
    /// 直接调用 == 操作符，减少错误的发生
    return(!(this->operator ==(rM)));
}

CMatrix& CMatrix::operator /= (const double value)
{
    /// 防止除零
    if(0 != value)
    {
        double dInv = 1.0/value;
        /// 直接调用 *= 操作符，减少错误的发生
        return (this->operator *=(dInv));
    }

    return (*this);
}

CMatrix& CMatrix::operator *=(const double value)
{
    /// 将所有的元素进行 缩放
    for(unsigned int i=0; i<m_unRow; ++i)
    {
        for(unsigned int j=0; j<m_unCol; ++j)
        {
            m_ppdM[i][j] *= value;
        }
    }

    return (*this);
}

void CMatrix::Empty()
{
    /// 如果缓存为空，则直接返回
    if(0 == m_ppdM)
    {
        InitAttr();
        return;
    }

    /// 删除空间
    for(unsigned int i=0; i<m_unRow; ++i)
    {
        CMemPool::GetInstance()->Remove(m_ppdM[i]);
    }

    CMemPool::GetInstance()->Remove(m_ppdM);
    InitAttr();
}

void CMatrix::InitAttr()
{
    m_unRow = m_unCol = 0;
    m_ppdM = 0;
}
