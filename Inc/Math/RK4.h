#ifndef YTY_RK4_H
#define YTY_RK4_H
/*****************************************
  作用：四阶龙格库塔数值积分法
  注意：此类是非线程安全的，在多线程中使用需要对变量
       进行加锁保护，否则运算结果不能保证
 *****************************************/

#include "Math_global.h"
namespace Math{
class CVector;
}

namespace Numerical{
using namespace Math;
/**
 * 定义微分方程的函数指针
 */
typedef void (*RK4funct)(double x, const CVector& y
                       , CVector& yp,void* pAux);

class MATH_EXPORT CRK4
{
public:
    /**
     * @brief 构造函数
     * @param pfRK4 微分方程
     * @param nEqn 维数
     * @param pAux 积分需要的相关参数列表
     */
    CRK4(RK4funct pfRK4, int nEqn, void* pAux );
    ~CRK4();

    /// 积分步骤
    virtual void Step (double& dt, CVector& vY, double dh);
protected:

  /// 成员变量
  RK4funct       m_pfRK4;                 /// 龙格库塔微分方程
  int            m_nEqn;                  /// 维数
  void*          m_pAux;                  /// 附加参数
  CVector    m_vK1,m_vK2,m_vK3,m_vK4; /// 函数的四个斜率
};
}
#endif // YTY_RK4_H
