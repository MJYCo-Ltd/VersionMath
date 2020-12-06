#ifndef YTY_VECMAT_H
#define YTY_VECMAT_H
/*****************************************
  作用：进行矩阵的运算、向量的运算、矩阵与向量的运算
  注意：所有的操作不成功会返回一个空的向量或者矩阵，
       使用时需要注意，如果调用函数时不能保证输入
       数值的正确性，需要检测对象是否可为空，或者
       直接判断对象是否合理。例如：
       CVector tmpV;
       CMatrix tmpMat;
       …
       if(tmpMat) 或者 if(tmpMat.IsEmpty())
       {
         代码;
       }
       …
       if(tmpV) 或者 if(tmpV.IsEmpty())
       {
         代码;
       }
 *****************************************/

#include <iostream>
#include <Math/Vector.h>
#include <Math/Matrix.h>

namespace Math{
class MATH_EXPORT CVecMat
{
public:
    //////////////////////////////向量操作///////////////////////////////////////

    /**
     * @brief 根据输入检查输出是否合理
     * @param vIn  输入项
     * @param vOut 输出项
     * @return 返回检查结果
     * @attention 1、如果 vOut是一个空的向量，会自动开辟空间
     *            以匹配 vIn 的输入项
     *            2、如果 vOut 非空，则判断 vIn.Size() == vOut.Size()
     */
    static bool IsValid(const CVector& vIn, CVector& vOut);

    /**
     * @brief 批量转换 vIn 到 vOut
     * @param rMatrix 转换矩阵
     * @param vIn     输入项
     * @param vOut    输出项
     * @return 返回处理结果
     * @attention 没有对 vOut的数据大小进行检查
     *            必须先调用 IsValid 进行数据校验
     */
    static bool CalReault(const CMatrix& rMatrix, const CVector& vIn, CVector& vOut);

    /**
     * @brief 将两个向量合并成一个向量
     * @param 向量 a [1,2,3,4]
     * @param 向量 b [5]
     * @return 向量 a+b [1,2,3,4,5]
     */
    static CVector Stack(const CVector& a, const CVector& b);

    /// 通过球极坐标构建向量
    /// azim 方位角(azimuth) [rad]
    /// elev 仰角(elevation) [rad]
    static CVector VecPolar (double azim, double elev, double r=1.0);

    /**
     * @brief 根据站心坐标计算 俯仰角 方位角
     * @param vecLtc 站心坐标 [m]
     * @param dAzim  方位角  [rad]
     * @param dElev  仰角    [rad]
     * @attention 该方法与VecPolar非逆过程，用VecPolar构建向量时
     *            方位角须取反 再加 PI/2. 即(VecPolar(pi/2-dAzim,dElev,rRange))
     */
    static void AzEl(const CVector& vecLtc, double& dAzim, double& dElev);

    /// 向量点乘
    static double Dot (const CVector& left, const CVector& right);
    /// 向量的模
    static double Norm (const CVector& V);
    /// 向量的叉乘
    static CVector Cross (const CVector& left, const CVector& right);
    //////////////////////////////向量操作 end///////////////////////////////////////

    //////////////////////////////矩阵操作 //////////////////////////////////////////
    /// 通过向量构建一个对角矩阵
    static CMatrix Diag(const CVector& Vec);

    /// 构建一个单位矩阵
    /// @param Size 单位矩阵的维数
    static CMatrix Id(int Size);

    /// 通过两个向量相乘构建一个矩阵
    static CMatrix Dyadic (const CVector& left, const CVector& right);

    /// 构建旋转矩阵
    /// @param Angle 需要旋转的角度值 [rad]
    static CMatrix R_x(double dAngle);
    static CMatrix R_y(double dAngle);
    static CMatrix R_z(double dAngle);

    /// 矩阵的转置
    static CMatrix Transp(const CMatrix& rMat);

    /// 矩阵的逆
    static CMatrix Inv(const CMatrix& rMat);

    //////////////////////////////矩阵操作 end///////////////////////////////////////

private:
    CVecMat();
};
}

using namespace Math;
//////////////////////////////向量操作符///////////////////////////////////////
MATH_EXPORT CVector operator &(const CVector& a, double b);
MATH_EXPORT CVector operator &(double a, const CVector& b);
MATH_EXPORT CVector operator &(const CVector& a, const CVector& b);
MATH_EXPORT CVector operator * (double value, const CVector& V);
MATH_EXPORT CVector operator * (const CVector& V, double value);
MATH_EXPORT CVector operator / (const CVector& V, double value);
MATH_EXPORT CVector operator + (const CVector& left, const CVector& right);
MATH_EXPORT CVector operator - (const CVector& left, const CVector& right);
/// 将向量取负
MATH_EXPORT CVector operator - (const CVector& V);
//////////////////////////////向量操作符///////////////////////////////////////
/// 输出
MATH_EXPORT std::ostream& operator << (std::ostream& os, const CVector& rVec);
MATH_EXPORT std::ostream& operator << (std::ostream& os, const CMatrix& rMat);
//////////////////////////////矩阵操作 //////////////////////////////////////////
/// 以行合并，行数增加
MATH_EXPORT CMatrix operator &(const CMatrix& A, const CMatrix& B);

/// 以列合并，列数增加
MATH_EXPORT CMatrix operator |(const CMatrix& A, const CMatrix& B);

/// 缩放矩阵
MATH_EXPORT CMatrix operator * (double value, const CMatrix& Mat);
MATH_EXPORT CMatrix operator * (const CMatrix& Mat, double value);
MATH_EXPORT CMatrix operator / (const CMatrix& Mat, double value);

/// 矩阵取反
MATH_EXPORT CMatrix operator - (const CMatrix& Mat);

/// 矩阵的加减
MATH_EXPORT CMatrix operator + (const CMatrix& left, const CMatrix& right);
MATH_EXPORT CMatrix operator - (const CMatrix& left, const CMatrix& right);

/// 矩阵相乘
MATH_EXPORT CMatrix operator * (const CMatrix& left, const CMatrix& right);
//////////////////////////////矩阵操作 //////////////////////////////////////////
/// 以行合并，行数增加
MATH_EXPORT CMatrix operator &(const CMatrix& A, const CVector& Row);
MATH_EXPORT CMatrix operator &(const CVector& Row, const CMatrix& A);

/// 以列合并，列数增加
MATH_EXPORT CMatrix operator |(const CMatrix& A, const CVector& Col);
MATH_EXPORT CMatrix operator |(const CVector& Col, const CMatrix& A);

/// 矩阵乘以一个列向量
MATH_EXPORT CVector operator * (const CMatrix& Mat, const CVector& Vec);

/// 一个行向量乘以矩阵
MATH_EXPORT CVector operator * (const CVector& Vec, const CMatrix& Mat);

/// 通过两个向量相乘构建矩阵
MATH_EXPORT CMatrix operator * (const CVector& left, const CVector& right);

#endif // YTY_VECMAT_H
