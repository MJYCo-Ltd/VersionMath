#ifndef YTY_J2_H
#define YTY_J2_H

#include <Satellite/SAT_global.h>
#include <Math/Quaternion.h>
#include <Math/VecMat.h>


namespace Satellite{
class ALGORITHM_EXPORT CJ2
{
public:
    CJ2(double dA, double dE, double dI, double dRAAN, double dMA, double dAP);

    const CVector& CalPV(double dT);
private:
    double m_dA,m_dE,m_dI,m_dRAAN,m_dMA,m_dAP,m_dN,m_dFac,m_dV;
    double m_dOmega,m_dW,m_domega;
    CVector m_vPV;
};
}
#endif // YTY_J2_H
