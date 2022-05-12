#ifndef YTY_RKF78_H
#define YTY_RKF78_H
/*****************************************
  作用：7阶龙格库塔数值积分法
  注意：此类是非线程安全的，在多线程中使用需要对变量
       进行加锁保护，否则运算结果不能保证
 *****************************************/
#include <Math/RK4.h>

namespace Numerical{
class MATH_EXPORT CRKF78:public CRK4
{
public:
    /**
     * @brief 构造函数
     * @param pfRK4 微分方程
     * @param nEqn 维数
     * @param pAux 积分需要的相关参数列表
     */
    CRKF78(RK4funct pfRK4, int nEqn, void* pAux );
    ~CRKF78();

    virtual void Step (double& dt, Math::CVector& vY, double dh);
private:
    Math::CVector m_vK[13];

    const static double m_dAK[13];
    const static double m_dW[13];

    const static double m_dB1;
    const static double m_dB2[2];
    const static double m_dB3[3];
    const static double m_dB4[4];
    const static double m_dB5[5];
    const static double m_dB6[6];
    const static double m_dB7[7];
    const static double m_dB8[8];
    const static double m_dB9[9];
    const static double m_dB10[10];
    const static double m_dB11[11];
    const static double m_dB12[12];
};
}
#endif // YTY_RKF78_H
