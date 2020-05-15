#include <iostream>
#include <iomanip>
#include <cmath>
using namespace std;
#include "VecMat.h"
#include "Matrix.h"
#include "Vector.h"
#include "sofam.h"
#include "MathCommonAlgorithm.h"

using namespace Math;

CVecMat::CVecMat()
{
}

/// 合并
CVector CVecMat::Stack (const CVector& a, const CVector& b)
{
    return (a&b);
}

/// 通过球极坐标构建向量
CVector CVecMat::VecPolar (double azim, double elev, double r)
{
    return CVector(r*cos(azim)*cos(elev),r*sin(azim)*cos(elev),r*sin(elev));
}

void CVecMat::AzEl(const CVector &vecLtc, double &dAzim, double &dElev)
{
    dAzim = atan2(vecLtc(0),vecLtc(1));
    dAzim = ((dAzim<0.0)? dAzim+D2PI : dAzim);
    dElev = atan ( vecLtc(2) / sqrt(vecLtc(0)*vecLtc(0)+vecLtc(1)*vecLtc(1)) );
}

/// 点乘
double CVecMat::Dot (const CVector& left, const CVector& right)
{
    int nLeftDim = left.Size();
    if (right.Size() != nLeftDim)
    {
        cerr << "ERROR: Incompatible shape in Dot(Vector,Vector)" << endl;
        return(0.0);
    }

    double Sum = 0.0;
    for (int i=0; i<nLeftDim; ++i)
    {
        Sum+=left(i)*right(i);
    }

    return(Sum);
}

/// 计算模
double CVecMat::Norm (const CVector& V)
{
    return V.Length();
}

/// 向量的叉乘
CVector CVecMat::Cross (const CVector& left, const CVector& right)
{
    CVector Result;

    if ( (left.Size()!=3) || (right.Size()!=3) )
    {
        cerr << "ERROR: Invalid dimension in Cross(Vector,Vector)" << endl;
        return (Result);
    }

    Result.Resize(3);

    Result(0) = left(1)*right(2) - left(2)*right(1);
    Result(1) = left(2)*right(0) - left(0)*right(2);
    Result(2) = left(0)*right(1) - left(1)*right(0);
    return (Result);
}
/// /////////////////////// 向量操作 ///////////////////////////////////////////////

/// /////////////////////// 矩阵操作 ///////////////////////////////////////////////

CMatrix CVecMat::Diag(const CVector& Vec)
{
    /// 获取向量的维数
    int nSize = Vec.Size();

    /// 设置向量
    CMatrix Mat(nSize,nSize);

    for (int i=0; i<nSize; ++i)
    {
        Mat(i,i) = Vec(i);
    }

    return Mat;
}

/// 构建单位矩阵
CMatrix CVecMat::Id(int Size)
{
    /// 开辟空间
    CMatrix Aux(Size,Size);

    /// 设置对角线为1.0
    for (int i=0; i<Size; ++i)
    {
        Aux(i,i) = 1.0;
    }

    return Aux;
}

/// 通过向量相乘构建矩阵
CMatrix CVecMat::Dyadic(const CVector &left, const CVector &right)
{
    return (left*right);
}

/// 绕X轴旋转
CMatrix CVecMat::R_x(double dAngle)
{
    const double C = cos(dAngle);
    const double S = sin(dAngle);
    CMatrix U(3,3);
    U(0,0) = 1.0;  U(0,1) = 0.0;  U(0,2) = 0.0;
    U(1,0) = 0.0;  U(1,1) =  +C;  U(1,2) =  +S;
    U(2,0) = 0.0;  U(2,1) =  -S;  U(2,2) =  +C;
    return U;
}

/// 绕Y轴旋转
CMatrix CVecMat::R_y(double dAngle)
{
    const double C = cos(dAngle);
    const double S = sin(dAngle);
    CMatrix U(3,3);
    U(0,0) =  +C;  U(0,1) = 0.0;  U(0,2) =  -S;
    U(1,0) = 0.0;  U(1,1) = 1.0;  U(1,2) = 0.0;
    U(2,0) =  +S;  U(2,1) = 0.0;  U(2,2) =  +C;
    return U;
}

/// 绕Z轴旋转
CMatrix CVecMat::R_z(double dAngle)
{
    const double C = cos(dAngle);
    const double S = sin(dAngle);
    CMatrix U(3,3);
    U(0,0) =  +C;  U(0,1) =  +S;  U(0,2) = 0.0;
    U(1,0) =  -S;  U(1,1) =  +C;  U(1,2) = 0.0;
    U(2,0) = 0.0;  U(2,1) = 0.0;  U(2,2) = 1.0;
    return U;
}

CMatrix CVecMat::Transp(const CMatrix &rMat)
{
    /// 以指定矩阵的行为列，列为行
    int nRow = rMat.Col();
    int nCol = rMat.Row();

    /// 复制原来的数据
    CMatrix T(nRow,nCol);
    for ( int i=0; i<nRow; ++i)
    {
        for ( int j=0; j<nCol; ++j)
        {
            T(i,j) = rMat(j,i);
        }
    }

    return (T);
}

CMatrix CVecMat::Inv(const CMatrix &rMat)
{
    const int nRow = rMat.Row();

    CMatrix Inverse;

    if (!rMat.IsSquare())
    {
        cerr << "ERROR: Invalid shape in Inv(Matrix)" << endl;
        return(Inverse);
    }

    CMatrix LU(nRow,nRow);
    CVector b(nRow), Indx(nRow);
    // LU decomposition

    LU = rMat;
    if(!LU_Decomp ( LU, Indx ))
    {
        return (Inverse);
    }

    // Solve Ax=b for  unit vectors b_1..b_n

    Inverse.Resize(nRow,nRow);
    for (int j=0; j<nRow; ++j)
    {
        b=0.0; b(j)= 1.0;                     // Set b to j-th unit vector
        LU_BackSub ( LU, Indx, b );           // Solve Ax=b
        Inverse.SetCol(j,b);                  // Copy result
    }

    return (Inverse);
}

/// /////////////////////// 矩阵操作 end////////////////////////////////////////////
/// 合并操作
CVector operator &(const CVector& a, double b)
{
    int    n=a.Size();
    CVector tmp(n+1);

    for(int i=0; i<n; ++i)
    {
        tmp(i) = a(i);
    }
    tmp(n) = b;
    return tmp;
}

CVector operator &(double a, const CVector& b)
{
    int    n=b.Size();
    CVector tmp(n+1);
    tmp(0) = a;
    for (int i=1;i<=n;++i)
    {
        tmp(i)=b(i-1);
    }
    return tmp;
}

CVector operator &(const CVector& a, const CVector& b)
{
    int i;
    int nADim = a.Size();
    int nBDim = b.Size();
    CVector c(nADim+nBDim);
    for (i=0;i<nADim;++i)
    {
        c(i)=a(i);
    }

    for (i=0;i<nBDim;++i)
    {
        c(i+nADim)=b(i);
    }

    return c;
}

/// 对向量进行缩放
CVector operator * (double value, const CVector& V)
{
    CVector Aux = V;
    Aux *= value;
    return Aux;
}

CVector operator * (const CVector& V, double value)
{
    return value*V;
}

CVector operator / (const CVector& V, double value)
{
    CVector Aux = V;

    Aux /= value;

    return Aux;
}


/// 向量取负
CVector operator - (const CVector& V)
{
    int nSize = V.Size();
    CVector Aux(nSize);

    /// 向量取负
    for (int i=0; i<nSize; ++i)
    {
        Aux(i)=-V(i);
    }

    return Aux;
}


/// 向量的加法
CVector operator + (const CVector& left, const CVector& right)
{
    CVector Aux = left;
    Aux += right;
    return (Aux);
}

/// 向量的减法
CVector operator - (const CVector& left, const CVector& right)
{
    CVector Aux = left;
    Aux -= right;
    return (Aux);
}
/// /////////////////////// 向量操作 end///////////////////////////////////////////////

/// 以行合并，行数增加
CMatrix operator &(const CMatrix& A, const CVector& Row)
{
    int    nRow=A.Row();
    int    nCol=A.Col();
    CMatrix tmp;

    if ( nCol!=Row.Size() )
    {
        cerr << "ERROR: Incompatible shape in Matrix&Vector concatenation" << endl;
        return(tmp);
    }

    tmp.Resize(nRow+1,nCol);

    for (int j=0;j<nCol;++j)
    {
        for (int i=0;i<nRow;++i)
        {
            tmp(i,j)=A(i,j);
        }
        tmp(nRow,j) = Row(j);
    }
    return (tmp);
}

CMatrix operator &(const CVector& Row, const CMatrix& A)
{
    int    nRow=A.Row();
    int    nCol=A.Col();
    CMatrix tmp;
    if ( nCol!=Row.Size() )
    {
        cerr << "ERROR: Incompatible shape in Vector&Matrix concatenation" << endl;
        return(tmp);
    }

    tmp.Resize(nRow+1,nCol);
    for (int j=0;j<nCol;++j)
    {
        tmp(0,j) = Row(j);
        for (int i=0;i<nRow;++i)
        {
            tmp(i+1,j)=A(i,j);
        }
    }
    return (tmp);
}

CMatrix operator &(const CMatrix& A, const CMatrix& B)
{
    /// 获取需要的信息
    int nACol = A.Col();
    int nARow = A.Row();
    int nBRow = B.Row();
    int nBCol = B.Col();

    CMatrix tmp;
    if ( nACol!=nBCol)
    {
        cerr << "ERROR: Incompatible shape in Matrix&Matrix concatenation" << endl;
        return (tmp);
    }

    tmp.Resize(nARow+nBRow,nACol);

    int    i;

    for (int j=0;j<nACol;++j)
    {
        for (i=0;i<nARow;++i)
        {
            tmp(i,j)=A(i,j);
        }
        for (i=0;i<nBRow;++i)
        {
            tmp(i+nARow,j)=B(i,j);
        }
    }
    return (tmp);
}

/// 以列合并，列数增加
CMatrix operator |(const CMatrix& A, const CVector& Col)
{
    int    nRow=A.Row();
    int    nCol=A.Col();
    CMatrix tmp;
    if ( nRow!=Col.Size() )
    {
        cerr << "ERROR: Incompatible shape in Matrix|Vector concatenation" << endl;
        return(tmp);
    }

    /// 重置大小
    tmp.Resize(nRow,nCol+1);
    for (int i=0;i<nRow;++i)
    {
        for (int j=0;j<nCol;++j)
        {
            tmp(i,j)=A(i,j);
        }
        tmp(i,nCol) = Col(i);
    }
    return (tmp);
}

CMatrix operator |(const CVector& Col, const CMatrix& A)
{
    int    nRow=A.Row();
    int    nCol=A.Col();
    CMatrix tmp;
    if ( nRow!=Col.Size() )
    {
        cerr << "ERROR: Incompatible shape in Vector|Matrix concatenation" << endl;
        return(tmp);
    }

    /// 重置大小
    tmp.Resize(nRow,nCol+1);
    for (int i=0;i<nRow;++i)
    {
        tmp(i,0) = Col(i);
        for (int j=0;j<nCol;++j)
        {
            tmp(i,j+1)=A(i,j);
        }
    }
    return (tmp);
}

CMatrix operator |(const CMatrix& A, const CMatrix& B)
{

    int    j;
    int nARow = A.Row();
    int nACol = A.Col();
    int nBRow = B.Row();
    int nBCol = B.Col();
    CMatrix tmp;
    if ( nARow != nBRow )
    {
        cerr << "ERROR: Incompatible shape in Matrix|Matrix concatenation" << endl;
        return(tmp);
    }

    /// 重设大小
    tmp.Resize(nARow,nACol+nBCol);
    for (int i=0;i<nARow;++i)
    {
        for (j=0;j<nACol;++j)
        {
            tmp(i,j)=A(i,j);
        }

        for (j=0;j<nBCol;++j)
        {
            tmp(i,j+nACol)=B(i,j);
        }
    }
    return (tmp);
}

/// 输出向量
ostream& operator << (ostream& os, const CVector& rVec)
{
    int w = os.width();
    for (int i=0; i<rVec.Size(); ++i)
    {
        os << setw(w) << rVec(i);
    }
    os << endl;
    return os;
}
/// 输出矩阵
ostream& operator << (ostream& os, const CMatrix& rMat)
{
    int w = os.width();
    for (int i=0; i<rMat.Row(); ++i)
    {
        for (int j=0; j<rMat.Col(); ++j)
        {
            os << setw(w) << rMat(i,j);
        }
        os << endl;
    }
    return os;
}
/// 缩放矩阵
CMatrix operator * (double value, const CMatrix& Mat)
{
    CMatrix Aux = Mat;
    Aux *= value;

    return Aux;
}

CMatrix operator * (const CMatrix& Mat, double value)
{
    return value*Mat;
}

/// 缩放矩阵
CMatrix operator / (const CMatrix& Mat, double value)
{
    CMatrix Aux = Mat;
    Aux /= value;

    return Aux;
}


/// 取负
CMatrix operator - (const CMatrix& Mat)
{
    int nRow = Mat.Row();
    int nCol = Mat.Col();

    CMatrix Aux(nRow,nCol);
    for (int i=0; i<nRow; ++i)
    {
        for (int j=0; j<nCol; ++j)
        {
            Aux(i,j)=-Mat(i,j);
        }
    }

    return Aux;
}


/// 矩阵的加减
CMatrix operator + (const CMatrix& left, const CMatrix& right)
{
    CMatrix Aux = left;
    Aux += right;
    return Aux;
}

CMatrix operator - (const CMatrix& left, const CMatrix& right)
{
    CMatrix Aux = left;
    Aux -= right;
    return Aux;
}


/// 矩阵的乘法
CMatrix operator * (const CMatrix& left, const CMatrix& right)
{
    int nLRow = left.Row();
    int nLCol = left.Col();
    int nRRow = right.Row();
    int nRCol = right.Col();

    CMatrix Aux;

    if (nLCol!=nRRow)
    {
        cerr << "ERROR: Incompatible shape in *(Matrix,Matrix)" << endl;
        return Aux;
    }
    Aux.Resize(nLRow,nRCol);

    double Sum;
    for (int i=0; i<nLRow; ++i)
    {
        for (int j=0; j<nRCol; ++j)
        {
            Sum = 0.0;
            for (int k=0; k<nLCol; ++k)
            {
                Sum += left(i,k) * right(k,j);
            }
            Aux(i,j) = Sum;
        }
    }
    return Aux;
}


/// 向量与矩阵的乘法
CVector operator * (const CMatrix& Mat, const CVector& Vec)
{
    /// 获取矩阵的 行数 列数
    int nRow = Mat.Row();
    int nCol = Mat.Col();

    CVector Aux;
    if (Vec.Size() != nCol)
    {
        cerr << "ERROR: Incompatible shape in *(Matrix,Vector)" << endl;
        return (Aux);
    }

    /// 设置成矩阵的行
    Aux.Resize(nRow);
    double Sum;
    for (int i=0; i<nRow; ++i)
    {
        Sum = 0.0;
        for (int j=0; j<nCol; ++j)
        {
            Sum += Mat(i,j) * Vec(j);
        }
        Aux(i) = Sum;
    }
    return Aux;
}

CVector operator * (const CVector& Vec, const CMatrix& Mat)
{
    /// 获取矩阵的 行数 列数
    int nRow = Mat.Row();
    int nCol = Mat.Col();

    CVector Aux;

    /// 如果失败返回空的向量
    if (Vec.Size() != nRow)
    {
        cerr << "ERROR: Incompatible shape in *(Vector,Matrix)" << endl;
        return(Aux);
    }

    /// 设置成矩阵的列
    Aux.Resize(nCol);
    double Sum;
    for (int j=0; j<nCol; ++j)
    {
        Sum = 0.0;
        for (int i=0; i<nRow; ++i)
        {
            Sum += Vec(i) * Mat(i,j);
        }
        Aux(j) = Sum;
    }
    return Aux;
}

CMatrix operator *(const CVector& left, const CVector& right)
{
    /// 用左向量的维数作为行数
    int nRow = left.Size();

    /// 用第二个向量的维数作为列数
    int nCol = right.Size();

    /// 构建矩阵
    CMatrix Mat(nRow,nCol);
    for (int i=0;i<nRow;++i)
    {
        for (int j=0;j<nCol;++j)
        {
            Mat(i,j) = left(i)*right(j);
        }
    }
    return Mat;
}
