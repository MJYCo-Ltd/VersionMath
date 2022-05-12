#ifndef YTY_COMMONALGORITHM_MATH_H
#define YTY_COMMONALGORITHM_MATH_H

/*****************************************
  作用：各个需要的基础算法，通用算法
  注意：此类是非线程安全的，在多线程中使用需要对变量
       进行加锁保护，否则运算结果不能保证
       对外没有接口，让外部使用，只能内部使用
 *****************************************/

#include <Math/VecMat.h>

namespace Math{
////////////基础数学////////////

/**
 * @brief 将角度值限定到[0~2PI)中
 * @param angle_in_radians 角度值 [rad]
 * @return 在[0~2PI) 中的角度值 [rad]
 */
double MATH_EXPORT zero_to_two_pi( double angle_in_radians);
/**
 * @brief 保证a的正负与b相同
 * @param a
 * @param b
 * @return
 * @attention b>=0  |a|
 *             b<0 -|a|
 */
double MATH_EXPORT sign(double a, double b);

/**
 * @brief y=x-[x]
 */
double MATH_EXPORT Frac (double x);

/**
 * @brief y%x
 */
double MATH_EXPORT Modulo(double x, double y);

/**
 * @brief arccos
 * @param arg
 * @return 角度值 [rad]
 * @attention 相较与acos进行了边界值的判断
 */
double MATH_EXPORT acose(double arg);

/**
 * @brief arcsin
 * @param arg
 * @return 角度值 [rad]
 * @attention 相较与asin进行了边界值的判断
 */
double MATH_EXPORT asine(double arg);
/**
 * @brief LU_BackSub
 * @param A
 * @param Indx
 * @param b
 */
void MATH_EXPORT LU_BackSub( CMatrix& A, CVector& Indx, CVector& b );

/**
 * @brief LU_Decomp
 * @param A
 * @param Indx
 * @return
 */
bool MATH_EXPORT LU_Decomp ( CMatrix& A, CVector& Indx );

/**
 * @brief 勒让德多项式
 * @param LL  阶数
 * @param x
 * @param PL 归一化勒让德多项式的值
 */
void MATH_EXPORT Legendre_sphPl(const int LL,const double x,double PL[]);

/**
 * @brief 缔合勒让德多项式
 * @param LL  阶数
 * @param x
 * @param PL 归一化勒让德多项式的值
 */
void MATH_EXPORT Legendre_sphPlm(const int LL,const double x,double PLM[][71]);

/**
 * @brief sin(m*x)和cos(m*x)
 * @param LL   阶数
 * @param S1X  sin(x)
 * @param C1X  cos(x)
 * @param SX   返回sin(m*x)数组
 * @param CX   返回cos(m*x)数组
 */
void MATH_EXPORT SmxCmx(const int LL,const double S1X,const double C1X,double* SX,double* CX);
////////////基础数学 end////////////
}
#endif // YTY_COMMONALGORITHM_MATH_H

