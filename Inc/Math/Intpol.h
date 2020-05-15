#ifndef YTY_INTPOL_H
#define YTY_INTPOL_H
/*****************************************
  作用：封装插值算法
  注意：此类是非线程安全的，在多线程中使用需要对变量
       进行加锁保护，否则运算结果不能保证
  备注：包含的差值算法有 Neville Lagrange Newton 三次样条插值
       目前不知道什么时候使用什么插值算法，比较好
 *****************************************/
#include "Math_global.h"
#include <vector>
namespace Math {
class CVector;
}
namespace Numerical {

using namespace Math;
using namespace std;

class MATH_EXPORT Cntpol
{
public:
    Cntpol();
    ~Cntpol();

    /**
     * @brief 内维尔插值算法(Neville)
     * @param nNum  数组的个数
     * @param pdX   自变量数组  x
     * @param pdY   因变量数组 y=f(x)
     * @param dX0   需要计算的自变量的值
     * @return      需要计算的自变量对应的因变量
     */
    static double ItNeville(int nNum,const double *pdX,const double *pdY,double dX0);

    /**
     * @brief 拉格郎日插值 (Lagrange)
     * @param nNum  数组的个数
     * @param pdX   自变量数组  x
     * @param pdY   因变量数组 y=f(x)
     * @param dX0   需要计算的自变量的值
     * @return      需要计算的自变量对应的因变量
     */
    static double ItLagrange(int nNum,const double *pdX,const double *pdY, double dX0);

    /**
     * @brief 牛顿插值 (Newton)
     * @param nNum  数组的个数
     * @param pdX   自变量数组  x
     * @param pdY   因变量数组 y=f(x)
     * @param dX0   需要计算的自变量的值
     * @return      需要计算的自变量对应的因变量
     */
    static double ItNewton(int nNum,const double *pdX,const double *pdY, double dX0);

    /**
     * @brief 三次样条插值
     * @param nNum 数组的个数
     * @param pdX  自变量数组  x
     * @param pdY  因变量数组 y=f(x)
     * @param dX0  需要计算的自变量的值
     * @return     需要计算的自变量对应的因变量
     */
    static double ItCubicSpline(int nNum, const double*pdX, const double*pdY,double dX0);

    /**
     * @brief 三次样条插值
     * @param vX  自变量数组X
     * @param vY  因变量数组Y
     * @param dX  需要计算的自变量的值
     * @return    通过三次样条插值，计算因变量。
     */
    static CVector ItCubicSpline(const vector<double>& vX, const vector<CVector>& vY, double dX);
};
}
#endif // YTY_INTPOL_H
