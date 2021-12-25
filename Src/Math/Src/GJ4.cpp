#include <Math/VecMat.h>
#include <Math/GJ4.h>
using namespace Numerical;
using namespace Math;

CGJ4::CGJ4(GJ4Pfunct pfGJ4, int nEqn, void *pAux):
    m_nEqn(nEqn), m_pfGJ4(pfGJ4), m_pAux(pAux)
{
}

//
// 4th order Runge-Kutta step for 2nd order differential equation
//

void CGJ4::RK4(double& dT, CVector& vecR, CVector& vecV, double dH)
{
    CVector v_1, v_2, v_3, v_4;
    CVector a_1, a_2, a_3, a_4;

    v_1 = vecV;               m_pfGJ4( dT      , vecR              , v_1, a_1, m_pAux );
    v_2 = vecV+(dH/2.0)*a_1;  m_pfGJ4( dT+dH/2.0, vecR+(dH/2.0)*v_1, v_2, a_2, m_pAux );
    v_3 = vecV+(dH/2.0)*a_2;  m_pfGJ4( dT+dH/2.0, vecR+(dH/2.0)*v_2, v_3, a_3, m_pAux );
    v_4 = vecV+dH*a_3;        m_pfGJ4( dT+dH    , vecR+dH*v_3      , v_4, a_4, m_pAux );

    dT = dT + dH;
    vecR = vecR + (dH/6.0)*( v_1 + 2.0*v_2 + 2.0*v_3 + v_4 );
    vecV = vecV + (dH/6.0)*( a_1 + 2.0*a_2 + 2.0*a_3 + a_4 );

}


//
// Initialization of backwards differences from initial conditions
//

void CGJ4::Init(double dT0, const CVector& vecR0, const CVector& vecV0, double dH)
{
    // Order of method

    const int m = 4;

    // Coefficients gamma/delta of 1st/2nd order Moulton/Cowell corrector method

    const double gc[m+1] = {+1.0, -1/2.0, -1/12.0, -1/24.0, -19/720.0 };
    const double dc[m+2] = {+1.0,   -1.0, +1/12.0,     0.0,  -1/240.0, -1/240.0 };

    int       i,j;
    double    t = dT0;
    CVector    r = vecR0;
    CVector    v = vecV0;

    // Save step size

    m_dH = dH;

    // Create table of accelerations at past times t-3h, t-2h, and t-h using
    // RK4 steps

    m_pfGJ4(t,r,v,m_vecD[0],m_pAux);     // D[i]=a(t-ih)
    for (i=1;i<=m-1;++i)
    {
        RK4(t,r,v,-m_dH);
        m_pfGJ4(t,r,v,m_vecD[i],m_pAux);
    }

    // Compute backwards differences

    for (i=1;i<=m-1;++i)
    {
        for (j=m-1;j>=i;--j)
        {
            m_vecD[j] = m_vecD[j-1]-m_vecD[j];
        }
    }

    // Initialize backwards sums using 4th order GJ corrector

    m_vecS1 = vecV0/m_dH;

    for (i=1;i<=m  ;++i)
    {
        m_vecS1 -= gc[i]*m_vecD[i-1];
    }

    m_vecS2 = vecR0/(m_dH*m_dH)-dc[1]*m_vecS1;

    for (i=2;i<=m+1;++i)
    {
        m_vecS2 -= dc[i]*m_vecD[i-2];
    }

}


//
// Step from t to t+h
//

void CGJ4::Step(double& dT, CVector& vecR, CVector& vecV)
{
    // Order of method

    const int m = 4;

    // Coefficients gamma/delta of 1st/2nd order Bashforth/Stoermr predictor

    const double gp[m+1] = {+1.0, +1/2.0, +5/12.0,  +3/8.0, +251/720.0 };
    const double dp[m+2] = {+1.0,    0.0, +1/12.0, +1/12.0,  +19/240.0,  +3/40.0 };

    int i;

    // 4th order predictor

    m_vecRp = dp[0]*m_vecS2;

    for(i=2;i<=m+1;++i)
    {
        m_vecRp += dp[i]*m_vecD[i-2];
    }

    m_vecRp = (m_dH*m_dH)*m_vecRp;

    m_vecVp = gp[0]*m_vecS1;

    for(i=1;i<=m  ;++i)
    {
        m_vecVp += gp[i]*m_vecD[i-1];
    }

    m_vecVp =     m_dH*m_vecVp;

    // Update backwards difference table

    m_pfGJ4 ( dT+m_dH, m_vecRp,m_vecVp, m_vecd[0], m_pAux );               // Acceleration at t+h

    for (i=1;i<=m-1;++i) // New differences at t+h
    {
        m_vecd[i]=m_vecd[i-1]-m_vecD[i-1];
    }

    for (i=0;i<=m-1;++i)  // Update differences
    {
        m_vecD[i]=m_vecd[i];
    }

    m_vecS1 += m_vecd[0];  m_vecS2 += m_vecS1;                        // Update sums

    // Update independent variable and solution

    dT = dT + m_dH;
    vecR = m_vecRp;
    vecV = m_vecVp;

}
