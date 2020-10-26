#include "Matrix.h"
#include "Vector.h"
#include <iostream>
#include <cstring>
#include <limits>
#include <cmath>
using namespace std;
using namespace Math;
/// 默认构造函数
CMatrix::CMatrix ()
{
    InitAttr();
}

CMatrix::CMatrix (int dim1, int dim2)
    : m_nRow(dim1), m_nCol(dim2)
{
    /// 防止索引值的出现
    if(m_nRow < 1 || m_nCol < 1)
    {
        cerr << "ERROR: Incompatible shape in Matrix Construction" << endl;
        InitAttr();
        return;
    }

    int i,j;
    /// 开辟空间
    m_ppdM = new double*[m_nRow];
    for (i=0; i<m_nRow; ++i)
    {
        m_ppdM[i] = new double[m_nCol];
    }

    j = sizeof(**m_ppdM)*m_nCol;
    /// 初始化
    for (i=0; i<m_nRow; ++i)
    {
        memset(m_ppdM[i],0,j);
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

    int i,j;
    m_nRow = rM.m_nRow;
    m_nCol = rM.m_nCol;
    /// 开辟空间
    m_ppdM = new double*[m_nRow];
    for (i=0; i<m_nRow; ++i)
    {
        m_ppdM[i] = new double[m_nCol];
    }

    j = sizeof(**m_ppdM)*m_nCol;
    /// 初始化
    for (i=0; i<m_nRow; ++i)
    {
        /// 一行一行的复制
        memcpy(m_ppdM[i],rM.m_ppdM[i],j);
    }
}

/// 通过数组构造
CMatrix::CMatrix (const double* p, int dim1, int dim2)
{
    /// 防止数组指向非法地址，或者索引值错误
    if(0 == p || dim1 < 1 || dim2 < 1)
    {
        InitAttr();
        return;
    }

    int i,j;
    m_nRow = dim1;
    m_nCol = dim2;
    /// 开辟空间
    m_ppdM = new double*[m_nRow];
    for (i=0; i<m_nRow; ++i)
    {
        m_ppdM[i] = new double[m_nCol];
    }
    j = sizeof(**m_ppdM)*m_nCol;
    /// 初始化
    for (i=0; i<m_nRow; ++i)
    {
        memcpy(m_ppdM[i],p+i*m_nCol,j);
    }

}

CMatrix::CMatrix(const double p[3][3], int dim1, int dim2)
{
    /// 防止数组指向非法地址，或者索引值错误
    if(0 == p || dim1 < 1 || dim2 < 1)
    {
        InitAttr();
        return;
    }

    int i,j;
    m_nRow = dim1;
    m_nCol = dim2;
    /// 开辟空间
    m_ppdM = new double*[m_nRow];
    for (i=0; i<m_nRow; ++i)
    {
        m_ppdM[i] = new double[m_nCol];
    }
    j = sizeof(**m_ppdM)*m_nCol;
    /// 初始化
    for (i=0; i<m_nRow; ++i)
    {
        memcpy(m_ppdM[i],p[i],j);
    }
}

CMatrix::~CMatrix()
{
    Empty();
}

/// 重置矩阵大小
CMatrix& CMatrix::Resize(int nRow, int nCol)
{
    /// 如果行列相等则直接返回
    if (m_nRow==nRow && m_nCol==nCol)
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

    int    i,j,i_max,j_max;
    /// 开辟一个新的空间
    double **M_new = new double*[nRow];
    for (i=0; i<nRow; ++i)
    {
        M_new[i] = new double[nCol];
    }

    /// 复制已有的数据，如果大于原来的数据则置零
    i_max = ((nRow<m_nRow)? nRow : m_nRow);
    j_max = ((nCol<m_nCol)? nCol : m_nCol);

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
    for (i=0; i<m_nRow; ++i)
    {
        delete[] m_ppdM[i];
    }
    delete [] m_ppdM;

    /// 将指针指向新开辟的空间
    m_ppdM = M_new;
    m_nRow = nRow;
    m_nCol = nCol;
    return (*this);
}

bool CMatrix::IsId() const
{
    /// 判断是否为空 是否为方阵
    if(IsEmpty() || IsSquare())
    {
        return(false);
    }

    /// 判断对角线上的值是否都为1.0
    for(int i=0; i<m_nRow; ++i)
    {
        for(int j=0; j<m_nCol; ++j)
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

/// 赋值操作
CMatrix& CMatrix::operator=(const double value)
{
    for (int i=0; i<m_nRow; ++i)
    {
        for (int j=0; j<m_nCol; ++j)
        {
            m_ppdM[i][j]=value;
        }
    }
    return (*this);
}

CMatrix& CMatrix::operator=(const CMatrix& rM)
{
    int i,j;
    if (this == &rM)
    {
        return (*this);
    }

    /// 如果矩阵为空则开辟空间
    if (IsEmpty())
    {
        m_nRow = rM.m_nRow;
        m_nCol = rM.m_nCol;
        m_ppdM = new double*[m_nRow];
        for (i=0; i<m_nRow; ++i)
        {
            m_ppdM[i] = new double[m_nCol];
        }
    }

    /// 如果矩阵的行列不同则直接返回
    if ( (m_nRow!=rM.m_nRow) || (m_nCol!=rM.m_nCol) )
    {
        cerr << "ERROR: Incompatible shapes in Matrix operator=(Matrix)" << endl;
        return(*this);
    }

    /// 计算每一行的空间大小
    j = m_nCol * sizeof(**m_ppdM);
    for (i=0; i<m_nRow; ++i)
    {
        /// 一行一行的复制
        memcpy(m_ppdM[i],rM.m_ppdM[i],j);
    }
    return (*this);
}


/// 返回指定列数据
CVector CMatrix::GetCol(int j) const
{
    CVector Res(m_nRow);
    for (int i=0; i<m_nRow; ++i)
    {
        Res(i)=m_ppdM[i][j];
    }
    return Res;
}

/// 返回指定行数据
CVector CMatrix::GetRow(int i) const
{
    CVector Res(m_nCol);
    for (int j=0; j<m_nCol; ++j)
    {
        Res(j)=m_ppdM[i][j];
    }
    return Res;
}

CVector CMatrix::Diag() const
{
    CVector Vec;
    if (IsEmpty() || !IsSquare())
    {
        cerr << "ERROR: Invalid shape in Matrix.Diag()" << endl;
        return(Vec);
    }

    Vec.Resize(m_nRow);
    /// 获取对角线上的值
    for (int i=0; i<m_nRow; ++i)
    {
        Vec(i) = m_ppdM[i][i];
    }
    return (Vec);
}

double CMatrix::Trace() const
{
    return this->Trace(0,m_nRow-1);
}

double CMatrix::Trace(int low, int upp) const
{
    double tmp = 0.0;

    /// 不是方阵则退出
    if (!IsSquare())
    {
        cerr << "ERROR: Invalid shape in Matrix.Trace()" << endl;
        return(tmp);
    }

    /// 判断参数是否合法
    if (low<0 || m_nRow<=upp)
    {
        cerr << "ERROR: Invalid arguments in Matrix.Trace()" << endl;
        return(tmp);
    }

    for (int i=low; i<=upp; ++i)
    {
        tmp += m_ppdM[i][i];
    }
    return (tmp);
}

/// 获取矩阵的部分
CMatrix CMatrix::slice(int first_row, int last_row, int first_col, int last_col)
{
    CMatrix Aux;

    /// 判断边界值是否满足
    if (first_row<0 || last_row<first_row || m_nRow-1<last_row
     || first_col<0 || last_col<first_col || m_nCol-1<last_col)
    {
        cerr << "ERROR: Invalid arguments in Matrix.slice()" << endl;
        return(Aux);
    }

    /// 重置矩阵大小
    Aux.Resize(last_row-first_row+1,last_col-first_col+1);

    /// 获取值
    for (int i=0;i<=last_row-first_row;++i)
    {
        for (int j=0;j<=last_col-first_col;++j)
        {
            Aux(i,j) = m_ppdM[i+first_row][j+first_col];
        }
    }

    return (Aux);
}


bool CMatrix::SetCol(int j, const CVector& Col)
{
    /// 判断向量是否可以放入矩阵
    if (Col.Size()!=m_nRow)
    {
        cerr << "ERROR: Incompatible shapes in Matrix.SetCol()" << endl;
        return(false);
    }

    /// 判断列是否超出范围
    if (j<0 || m_nCol<=j)
    {
        cerr << "ERROR: Column index out of range in Matrix.SetCol()" << endl;
        return(false);
    }

    for (int i=0; i<m_nRow; ++i)
    {
        m_ppdM[i][j]=Col(i);
    }

    return(true);
}

bool CMatrix::SetRow(int i, const CVector& Row)
{
    /// 判断向量是否可以放入矩阵
    if (Row.Size()!=m_nCol)
    {
        cerr << "ERROR: Incompatible shapes in Matrix.SetRow()" << endl;
        return(false);
    }

    /// 判断行是否超出范围
    if (i<0 || m_nRow<=i)
    {
        cerr << "ERROR: Row index out of range in Matrix.SetRow()" << endl;
        return(false);
    }

    for (int j=0; j<m_nCol; ++j)
    {
        m_ppdM[i][j]=Row(j);
    }

    return(true);
}



/// 矩阵的加法减法
bool CMatrix::operator +=(const CMatrix& rM)
{
    /// 检查矩阵的行列是否相同
    if ( (m_nRow!=rM.m_nRow) || (m_nCol!=rM.m_nCol) )
    {
        cerr << "ERROR: Incompatible shape in Matrix operator+=(Matrix)" << endl;
        return(false);
    };

    /// 进行加法运算
    for (int i=0; i<m_nRow; ++i)
    {
        for (int j=0; j<m_nCol; ++j)
        {
            this->m_ppdM[i][j]+=rM.m_ppdM[i][j];
        }
    }

    return(true);
}

bool CMatrix::operator -=(const CMatrix& rM)
{
    /// 检查矩阵的行列是否相同
    if ( (m_nRow!=rM.m_nRow) || (m_nCol!=rM.m_nCol) )
    {
        cerr << "ERROR: Incompatible shape in Matrix operator-=(Matrix)" << endl;
        return(false);
    };

    /// 进行各值取反
    for (int i=0; i<m_nRow; ++i)
    {
        for (int j=0; j<m_nCol; ++j)
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
    if(m_nRow != rM.m_nRow || m_nCol != rM.m_nCol)
    {
        return(false);
    }

    /// 判断各值是否相同
    for(int i=0; i<m_nRow; ++i)
    {
        for(int j=0; j<m_nCol; ++j)
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
    for(int i=0; i<m_nRow; ++i)
    {
        for(int j=0; j<m_nCol; ++j)
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
    for(int i=0; i<m_nRow; ++i)
    {
        delete[] m_ppdM[i];
    }

    delete[] m_ppdM;
    InitAttr();
}

void CMatrix::InitAttr()
{
    m_nRow = m_nCol = 0;
    m_ppdM = 0;
}
