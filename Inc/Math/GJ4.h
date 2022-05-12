#ifndef YTY_GJ4_H
#define YTY_GJ4_H
/*****************************************
  作用：4阶Gauss-Jackson数值积分(定步长多步法)
  注意：此类是非线程安全的，在多线程中使用需要对变量
       进行加锁保护，否则运算结果不能保证
  备注：完全不懂如何进行积分的，还需要知道的人进行注释
 *****************************************/
#include <Math/Vector.h>

namespace Numerical{
typedef void (*GJ4Pfunct)(double t, const Math::CVector& r
                          , const Math::CVector& v, Math::CVector& a
                          , void* pAux);

class MATH_EXPORT CGJ4
{
public:
    /**
     * @brief 构造函数
     * @param pfGJ4 微分方程
     * @param nEqn 维数
     * @param pAux 积分需要的相关参数列表
     */
    CGJ4(GJ4Pfunct pfGJ4, int nEqn, void* pAux);

    /// 初始化
    void Init(double dT0, const Math::CVector& vecR0, const Math::CVector& vecV0, double dH);

    /// 积分步长设置
    void Step(double& dT,Math::CVector& vecR,Math::CVector& vecV);

private:

    // 4th order Runge-Kutta step
    void RK4(double& dT,Math::CVector& vecR,Math::CVector& vecV, double dH);

    /// 成员变量
    int         m_nEqn;       // Dimension
    GJ4Pfunct   m_pfGJ4;           // Differential equation
    double      m_dH;           // Step size
    void*       m_pAux;        // Pointer to auxiliary data requird by f
    Math::CVector      m_vecS2,m_vecS1;       // First and second sum of acceleration
    Math::CVector      m_vecD[4];        // Backward differences of acceleration at t
    Math::CVector      m_vecd[4];        // Backward differences of acceleration at t+h
    Math::CVector      m_vecRp,m_vecVp;     // Predictor
};
}
#endif // YTY_GJ4_H
