#ifndef YTY_DE_H
#define YTY_DE_H
/*****************************************
  作用：多步法数值积分(Shampine & Gordon)(变阶变步长多步法)
  注意：此类是非线程安全的，在多线程中使用需要对变量
       进行加锁保护，否则运算结果不能保证
  备注：完全不懂如何进行积分的，还需要知道的人进行注释
 *****************************************/
#include "Math_global.h"

namespace Math{
class CVector;
class CMatrix;
}

namespace Numerical{

using namespace Math;

enum DE_STATE {
    DE_INIT     = 1,   // Restart integration
    DE_DONE     = 2,   // Successful step
    DE_BADACC   = 3,   // Accuracy requirement could not be achieved
    DE_NUMSTEPS = 4,   // Permitted number of steps exceeded
    DE_STIFF    = 5,   // Stiff problem suspected
    DE_INVPARAM = 6    // Invalid input parameters
};

typedef void (*DEfunct)(double x, const CVector& y
                        ,CVector& yp,void* pAux);

class MATH_EXPORT CDE
{
public:

    /**
     * @brief 构造函数
     * @param pfDE 微分方程
     * @param nEqn 维数
     * @param pAux 积分需要的相关参数列表
     */
    CDE(DEfunct pfDE, int nEqn, void* pAux);

    /// 初始化属性列表
    void Define(DEfunct pfDE, int nEqn, void* pAux);

    /**
     * @brief 积分初始化
     * @param dT0  开始积分时间
     * @param dRel 相对误差
     * @param dAbs 绝对误差
     */
    void Init(double dT0, double dRel, double dAbs);

    /**
     * @brief Integ 积分
     * @param dTout 结束时间
     * @param vecY  积分结果
     */
    void Integ(double dTout, CVector& vecY);

    /// 插入
    void Intrp(double dTout, CVector& vecY);
private:
    // Elementary integration step
    void Step(double& dX, CVector& vecY, double& dEps, bool& bCrash);

    // Interpolation
    void Intrp( double dXout, CVector& vecYout, CVector& vecYpout );

    // Integration (with full control of warnings and error status codes)
    void Integ(double& dt, double dTout, CVector& vecY);

    /// 成员变量
public:
    double       m_dRelerr;      /// 相对误差
    double       m_dAbserr;      /// 绝对误差
    DE_STATE     m_emState;      /// 计算状态
    bool         m_bPermitTOUT;  // Flag for integrating past tout
    double       m_dt;           // Value of independent variable
private:
    DEfunct  m_pfDE;
    int      m_nEqn;
    void*    m_pAux;
    CVector   m_vecYy,m_vecWt,m_vecP,m_vecYp,m_vecYpout;
    CMatrix   m_matPhi;
    double   m_dAlpha[13],m_dBeta[13],m_dV[13],m_dW[13],m_dPsi[13];
    double   m_dSig[14],m_dG[14];
    double   m_dX,m_dH,m_dHold,m_dTold,m_dDelsgn;
    int      m_nNs,m_nK,m_nKold;
    bool     m_bOldPermit, m_bPhase1,m_bStart,m_bNornd;
    bool     m_bInit;
};
}
#endif // YTY_DE_H
