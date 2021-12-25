#include <cmath>
#include <Math/VecMat.h>
#include <Math/RKF78.h>
using namespace Numerical;
using namespace Math;

const double CRKF78::m_dAK[13] =
{
    0, 2.0 / 27.0, 1.0 / 9.0, 1.0 / 6.0, 5.0 / 12.0, 1.0 / 2.0,
    5.0 / 6.0,  1.0 / 6.0, 2.0 / 3.0, 1.0 / 3.0, 1.0,		0.0, 1.0
};

const double CRKF78::m_dB1 = 2.0 / 27.0;
const double CRKF78::m_dB2[2] = {1.0/36.0, 1.0/12.0};
const double CRKF78::m_dB3[3] = { 1.0/24.0, 0.0, 0.125 };
const double CRKF78::m_dB4[4] = { 5.0/12.0, 0.0, -1.5625, 1.5625 };
const double CRKF78::m_dB5[5] = { 0.05,	0.0, 0.0, 0.25, 0.2 };
const double CRKF78::m_dB6[6] = { -25.0/108.0, 0.0, 0.0, 125.0/108.0, -65.0/27.0, 125.0/54.0 };
const double CRKF78::m_dB7[7] = {  31.0/300.0, 0.0, 0.0, 0.0, 61.0/225.0, -2.0/9.0, 13.0/900.0 };
const double CRKF78::m_dB8[8] = { 2.0, 0.0, 0.0, -53.0/6.0, 704.0/45.0, -107.0/9.0, 67.0/90.0, 3.0 };
const double CRKF78::m_dB9[9] = { -91.0/108.0, 0.0, 0.0, 23.0/108.0, -976.0/135.0, 311.0/54.0, -19.0/60.0, 17.0/6.0, -1.0/12.0 };
const double CRKF78::m_dB10[10] = { 2383.0/4100.0, 0.0, 0.0, -341.0/164.0, 4496.0/1025.0, -301.0/82.0, 2133.0/4100.0, 45.0/82.0, 45.0/164.0, 18.0/41.0 };
const double CRKF78::m_dB11[11] = { 3.0/205.0, 0.0, 0.0, 0.0, 0.0, -6.0/41.0, -3.0/205.0, -3.0/41.0, 3.0/41.0, 6.0/41.0, 0.0 };
const double CRKF78::m_dB12[12] = { -1777.0/4100.0, 0.0, 0.0, -341.0/164.0, 4496.0/1025.0, -289.0/82.0, 2193.0/4100.0, 51.0/82.0, 33.0/164.0, 12.0/41.0, 0.0, 1.0 };
const double CRKF78::m_dW[13] = { 0.0, 0.0, 0.0, 0.0, 0.0, 34.0/105.0, 9.0/35.0, 9.0/35.0, 9.0/280.0, 9.0/280.0, 0.0, 41.0/840.0, 41.0/840.0 };


CRKF78::CRKF78(RK4funct pfRK4, int nEqn, void *pAux):CRK4(pfRK4,nEqn,pAux)

{
}

CRKF78::~CRKF78()
{
}

void CRKF78::Step(double &dt, CVector &vY, double dh)
{
    double dt1;
    dt1 = dt;
    m_pfRK4(dt1, vY, m_vK[0], m_pAux );

    dt1 = dt + dh*m_dAK[1];

    m_pfRK4(dt1, vY+dh*m_dB1*m_vK[0], m_vK[1], m_pAux );

    dt1 = dt + dh*m_dAK[2];
    m_pfRK4( dt1 , vY+dh*(m_dB2[0]*m_vK[0] + m_dB2[1]*m_vK[1]), m_vK[2], m_pAux );

    dt1 = dt + dh*m_dAK[3];
    m_pfRK4( dt1 , vY+dh*(m_dB3[0]*m_vK[0] + m_dB3[2]*m_vK[2]), m_vK[3], m_pAux );

    dt1 = dt + dh*m_dAK[4];
    m_pfRK4( dt1 , vY+dh*(m_dB4[0]*m_vK[0] + m_dB4[2]*m_vK[2] + m_dB4[3]*m_vK[3]), m_vK[4], m_pAux );

    dt1 = dt + dh*m_dAK[5];
    m_pfRK4( dt1 , vY+dh*(m_dB5[0]*m_vK[0] + m_dB5[3]*m_vK[3] + m_dB5[4]*m_vK[4]) , m_vK[5], m_pAux );

    dt1 = dt + dh*m_dAK[6];
    m_pfRK4( dt1 , vY+dh*(m_dB6[0]*m_vK[0] + m_dB6[3]*m_vK[3] + m_dB6[4]*m_vK[4] + m_dB6[5]*m_vK[5]) , m_vK[6], m_pAux );

    dt1 = dt + dh*m_dAK[7];
    m_pfRK4( dt1 , vY+dh*(m_dB7[0]*m_vK[0] + m_dB7[4]*m_vK[4] + m_dB7[5]*m_vK[5] + m_dB7[6]*m_vK[6]), m_vK[7], m_pAux );

    dt1 = dt + dh*m_dAK[8];
    m_pfRK4( dt1 , vY+dh*(m_dB8[0]*m_vK[0] + m_dB8[3]*m_vK[3] + m_dB8[4]*m_vK[4] + m_dB8[5]*m_vK[5] + m_dB8[6]*m_vK[6] + m_dB8[7]*m_vK[7]), m_vK[8], m_pAux );

    dt1 = dt + dh*m_dAK[9];
    m_pfRK4( dt1 , vY+dh*(m_dB9[0]*m_vK[0] + m_dB9[3]*m_vK[3] + m_dB9[4]*m_vK[4] + m_dB9[5]*m_vK[5] + m_dB9[6]*m_vK[6] + m_dB9[7]*m_vK[7] + m_dB9[8]*m_vK[8]), m_vK[9], m_pAux );

    dt1 = dt + dh;  // ak[10] = 1.0
    m_pfRK4( dt1 , vY+dh*(m_dB10[0]*m_vK[0] + m_dB10[3]*m_vK[3] + m_dB10[4]*m_vK[4] + m_dB10[5]*m_vK[5] + m_dB10[6]*m_vK[6] + m_dB10[7]*m_vK[7] + m_dB10[8]*m_vK[8] + m_dB10[9]*m_vK[9]), m_vK[10], m_pAux );

    dt1 = dt;     // ak[11] = 0.0
    m_pfRK4( dt1 , vY+dh*(m_dB11[0]*m_vK[0] + m_dB11[5]*m_vK[5] + m_dB11[6]*m_vK[6] + m_dB11[7]*m_vK[7] + m_dB11[8]*m_vK[8] + m_dB11[9]*m_vK[9]), m_vK[11], m_pAux );

    dt1 = dt + dh;   // ak[12] = 1.0
    m_pfRK4( dt1 , vY+dh*(m_dB12[0]*m_vK[0] + m_dB12[3]*m_vK[3] + m_dB12[4]*m_vK[4] + m_dB12[5]*m_vK[5] + m_dB12[6]*m_vK[6] + m_dB12[7]*m_vK[7] + m_dB12[8]*m_vK[8] + m_dB12[9]*m_vK[9] + m_dB12[11]*m_vK[11]), m_vK[12], m_pAux );

    double EE=0.0;
    for(int ii=0; ii<6; ++ii)
    {
        EE += fabs((m_vK[0](ii) + m_vK[10](ii) - m_vK[11](ii) - m_vK[12](ii))*dh*41.0/840.0);
    }

    vY += dh*(m_dW[5]*m_vK[5] + m_dW[6]*m_vK[6] + m_dW[7]*m_vK[7] + m_dW[8]*m_vK[8] + m_dW[9]*m_vK[9] + m_dW[11]*m_vK[11] + m_dW[12]*m_vK[12]);
}

