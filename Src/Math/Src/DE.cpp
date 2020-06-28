#include <cmath>
#include <iostream>
#include <limits>
using namespace std;
#include "VecMat.h"
#include "DE.h"
#include "MathCommonAlgorithm.h"

using namespace Numerical;
using namespace Math;

/// 最大步数
static const int    MAXNUM = 500;

static const double umach = std::numeric_limits<double>::epsilon();
static const double twou   = 2.0*umach;
static const double fouru  = 4.0*umach;

#if _MSC_VER < 1700
double fmax(double fa1,double fa2)
{
    return(fa1>fa2 ? fa1 : fa2);
}

double fmin(double fa1,double fa2)
{
    return(fa1<fa2 ? fa1 : fa2);
}
#endif

CDE::CDE(DEfunct pfDE,int nEqn,void* pAux)
{
    Define(pfDE,nEqn,pAux);
}

void CDE::Define(DEfunct pfDE,int nEqn, void* pAux)
{
    m_nEqn          = nEqn;
    m_pfDE          = pfDE;
    m_pAux          = pAux;
    m_vecYy         = CVector(m_nEqn);   // Allocate vectors with proper dimension
    m_vecWt         = CVector(m_nEqn);
    m_vecP          = CVector(m_nEqn);
    m_vecYp         = CVector(m_nEqn);
    m_vecYpout      = CVector(m_nEqn);
    m_matPhi        = CMatrix(m_nEqn,17);
    m_emState       = DE_INVPARAM;     // Status flag
    m_bPermitTOUT   = true;            // Allow integration past tout by default
    m_dt            = 0.0;
    m_dRelerr       = 0.0;             // Accuracy requirements
    m_dAbserr       = 0.0;
}

/// Integration step
void CDE::Step (double& dX, CVector& vecY, double& dEps, bool& bCrash)
{

    // Constants

    // Powers of two (two(n)=2**n)
    static const double two[14] =
    { 1.0, 2.0, 4.0, 8.0, 16.0, 32.0, 64.0, 128.0,
      256.0, 512.0, 1024.0, 2048.0, 4096.0, 8192.0 };

    static const double gstr[14] =
    {1.0, 0.5, 0.0833, 0.0417, 0.0264, 0.0188,
     0.0143, 0.0114, 0.00936, 0.00789, 0.00679,
     0.00592, 0.00524, 0.00468 };

    // Variables

    bool    success;
    int     i,ifail, im1, ip1, iq, j, km1, km2, knew, kp1, kp2;
    int     l, limit1, limit2, nsm2, nsp1, nsp2;
    double  absh, erk, erkm1, erkm2, erkp1, err, hnew;
    double  p5eps, r, reali, realns, rho, round, sum, tau;
    double  temp1, temp2, temp3, temp4, temp5, temp6, xold;

    if (fabs(m_dH) < fouru*fabs(dX))
    {
        m_dH = sign(fouru*fabs(dX),m_dH);
        bCrash = true;
        return;                      // Exit
    }

    p5eps  = 0.5*dEps;
    bCrash  = false;
    m_dG[1]   = 1.0;
    m_dG[2]   = 0.5;
    m_dSig[1] = 1.0;

    ifail = 0;

    round = 0.0;

    for (l=0;l<m_nEqn;++l)
    {
        round += (vecY(l)*vecY(l))/(m_vecWt(l)*m_vecWt(l));
    }

    round = twou*sqrt(round);

    if (p5eps<round)
    {
        dEps = 2.0*round*(1.0+fouru);
        bCrash = true;
        return;
    }


    if (m_bStart)
    {
        // Initialize. Compute appropriate step size for first step.
        m_pfDE(dX, vecY, m_vecYp, m_pAux);
        sum = 0.0;

        for (l=0;l<m_nEqn;++l)
        {
            m_matPhi(l,1) = m_vecYp(l);
            m_matPhi(l,2) = 0.0;
            sum += (m_vecYp(l)*m_vecYp(l))/(m_vecWt(l)*m_vecWt(l));
        }

        sum  = sqrt(sum);
        absh = fabs(m_dH);

        if (dEps<16.0*sum*m_dH*m_dH)
        {
            absh=0.25*sqrt(dEps/sum);
        }

        m_dH    = sign(fmax(absh, fouru*fabs(dX)), m_dH);
        m_dHold = 0.0;
        hnew = 0.0;
        m_nK    = 1;
        m_nKold = 0;
        m_bStart  = false;
        m_bPhase1 = true;
        m_bNornd  = true;
        if (p5eps<=100.0*round)
        {
            m_bNornd = false;
            for (l=0;l<m_nEqn;++l)
            {
                m_matPhi(l,15)=0.0;
            }
        }
    }


    do
    {
        kp1 = m_nK+1;
        kp2 = m_nK+2;
        km1 = m_nK-1;
        km2 = m_nK-2;

        // ns is the number of steps taken with size h, including the
        // current one. When k<ns, no coefficients change.

        if (m_dH !=m_dHold)
        {
            m_nNs=0;
        }

        if (m_nNs<=m_nKold)
        {
            m_nNs=m_nNs+1;
        }

        nsp1 = m_nNs+1;

        if (m_nK>=m_nNs)
        {

            // Compute those components of alpha[*],beta[*],psi[*],sig[*]
            // which are changed
            m_dBeta[m_nNs] = 1.0;
            realns = m_nNs;
            m_dAlpha[m_nNs] = 1.0/realns;
            temp1 = m_dH*realns;
            m_dSig[nsp1] = 1.0;
            if (m_nK>=nsp1)
            {
                for (i=nsp1;i<=m_nK;++i)
                {
                    im1   = i-1;
                    temp2 = m_dPsi[im1];
                    m_dPsi[im1] = temp1;
                    m_dBeta[i]  = m_dBeta[im1]*m_dPsi[im1]/temp2;
                    temp1    = temp2 + m_dH;
                    m_dAlpha[i] = m_dH/temp1;
                    reali = i;
                    m_dSig[i+1] = reali*m_dAlpha[i]*m_dSig[i];
                }
            }
            m_dPsi[m_nK] = temp1;

            // Compute coefficients g[*]; initialize v[*] and set w[*].
            if (m_nNs>1)
            {
                // If order was raised, update diagonal part of v[*]
                if (m_nK>m_nKold)
                {
                    temp4 = m_nK*kp1;
                    m_dV[m_nK] = 1.0/temp4;
                    nsm2 = m_nNs-2;

                    for (j=1;j<=nsm2;++j)
                    {
                        i = m_nK-j;
                        m_dV[i] = m_dV[i] - m_dAlpha[j+1]*m_dV[i+1];
                    }
                }

                // Update V[*] and set W[*]
                limit1 = kp1 - m_nNs;
                temp5  = m_dAlpha[m_nNs];

                for (iq=1;iq<=limit1;++iq)
                {
                    m_dV[iq] = m_dV[iq] - temp5*m_dV[iq+1];
                    m_dW[iq] = m_dV[iq];
                }

                m_dG[nsp1] = m_dW[1];
            }
            else
            {
                for (iq=1;iq<=m_nK;++iq)
                {
                    temp3 = iq*(iq+1);
                    m_dV[iq] = 1.0/temp3;
                    m_dW[iq] = m_dV[iq];
                }
            }

            // Compute the g[*] in the work vector w[*]
            nsp2 = m_nNs + 2;
            if (kp1>=nsp2)
            {
                for (i=nsp2;i<=kp1;++i)
                {
                    limit2 = kp2 - i;
                    temp6  = m_dAlpha[i-1];

                    for (iq=1;iq<=limit2;++iq)
                    {
                        m_dW[iq] = m_dW[iq] - temp6*m_dW[iq+1];
                    }

                    m_dG[i] = m_dW[1];
                }
            }

        } // if K>=NS

        //
        // End block 1
        //

        //
        // Begin block 2
        //
        // Predict a solution p[*], evaluate derivatives using predicted
        // solution, estimate local error at order k and errors at orders
        // k, k-1, k-2 as if constant step size were used.
        //

        // Change phi to phi star
        if (m_nK>=nsp1)
        {
            for (i=nsp1;i<=m_nK;++i)
            {
                temp1 = m_dBeta[i];

                for (l=0;l<m_nEqn;++l)
                {
                    m_matPhi(l,i) = temp1 * m_matPhi(l,i);
                }
            }
        }

        // Predict solution and differences
        for (l=0;l<m_nEqn;++l)
        {
            m_matPhi(l,kp2) = m_matPhi(l,kp1);
            m_matPhi(l,kp1) = 0.0;
            m_vecP(l)       = 0.0;
        }

        for (j=1;j<=m_nK;++j)
        {
            i     = kp1 - j;
            ip1   = i+1;
            temp2 = m_dG[i];

            for (l=0; l<m_nEqn;++l)
            {
                m_vecP(l)     = m_vecP(l) + temp2*m_matPhi(l,i);
                m_matPhi(l,i) = m_matPhi(l,i) + m_matPhi(l,ip1);
            }
        }

        if (m_bNornd)
        {
            m_vecP = vecY + m_dH*m_vecP;
        }
        else
        {
            for (l=0;l<m_nEqn;++l)
            {
                tau = m_dH*m_vecP(l) - m_matPhi(l,15);
                m_vecP(l) = vecY(l) + tau;
                m_matPhi(l,16) = (m_vecP(l) - vecY(l)) - tau;
            }
        }

        xold = dX;
        dX = dX + m_dH;
        absh = fabs(m_dH);
        m_pfDE(dX, m_vecP, m_vecYp, m_pAux);

        // Estimate errors at orders k, k-1, k-2
        erkm2 = 0.0;
        erkm1 = 0.0;
        erk = 0.0;

        for (l=0;l<m_nEqn;++l)
        {
            temp3 = 1.0/m_vecWt(l);
            temp4 = m_vecYp(l) - m_matPhi(l,1);

            if (km2> 0)
            {
                erkm2 = erkm2 + ((m_matPhi(l,km1)+temp4)*temp3)
                       *((m_matPhi(l,km1)+temp4)*temp3);
            }

            if (km2>=0)
            {
                erkm1 = erkm1 + ((m_matPhi(l,m_nK)+temp4)*temp3)
                       *((m_matPhi(l,m_nK)+temp4)*temp3);
            }
            erk = erk + (temp4*temp3)*(temp4*temp3);
        }

        if (km2> 0)
        {
            erkm2 = absh*m_dSig[km1]*gstr[km2]*sqrt(erkm2);
        }

        if (km2>=0)
        {
            erkm1 = absh*m_dSig[m_nK]*gstr[km1]*sqrt(erkm1);
        }

        temp5 = absh*sqrt(erk);
        err = temp5*(m_dG[m_nK]-m_dG[kp1]);
        erk = temp5*m_dSig[kp1]*gstr[m_nK];
        knew = m_nK;

        // Test if order should be lowered
        if (km2 >0)
        {
            if (fmax(erkm1,erkm2)<=erk)
            {
                knew=km1;
            }
        }

        if (km2==0)
        {
            if (erkm1<=0.5*erk)
            {
                knew=km1;
            }
        }

        //
        // End block 2
        //


        //
        // If step is successful continue with block 4, otherwise repeat
        // blocks 1 and 2 after executing block 3
        //

        success = (err<=dEps);

        if (!success)
        {

            //
            // Begin block 3
            //

            // The step is unsuccessful. Restore x, phi[*,*], psi[*]. If
            // 3rd consecutive failure, set order to 1. If step fails more
            // than 3 times, consider an optimal step size. Double error
            // tolerance and return if estimated step size is too small
            // for machine precision.
            //

            // Restore x, phi[*,*] and psi[*]
            m_bPhase1 = false;
            dX = xold;

            for (i=1;i<=m_nK;++i)
            {
                temp1 = 1.0/m_dBeta[i];
                ip1 = i+1;
                for (l=0;l<m_nEqn;++l)
                {
                    m_matPhi(l,i)=temp1*(m_matPhi(l,i)-m_matPhi(l,ip1));
                }
            }

            if (m_nK>=2)
            {
                for (i=2;i<=m_nK;++i)
                {
                    m_dPsi[i-1] = m_dPsi[i] - m_dH;
                }
            }

            // On third failure, set order to one.
            // Thereafter, use optimal step size
            ++ifail;
            temp2 = 0.5;
            if (ifail>3)
            {
                if (p5eps < 0.25*erk)
                {
                    temp2 = sqrt(p5eps/erk);
                }
            }
            if (ifail>=3)
            {
                knew = 1;
            }
            m_dH = temp2*m_dH;
            m_nK = knew;
            if (fabs(m_dH)<fouru*fabs(dX))
            {
                bCrash = true;
                m_dH = sign(fouru*fabs(dX), m_dH);
                dEps *= 2.0;
                return;                     // Exit
            }

            //
            // End block 3, return to start of block 1
            //

        }  // end if(success)

    }
    while (!success);


    //
    // Begin block 4
    //
    // The step is successful. Correct the predicted solution, evaluate
    // the derivatives using the corrected solution and update the
    // differences. Determine best order and step size for next step.
    //

    m_nKold = m_nK;
    m_dHold = m_dH;


    // Correct and evaluate
    temp1 = m_dH*m_dG[kp1];
    if (m_bNornd)
    {
        for (l=0;l<m_nEqn;++l)
        {
            vecY(l) = m_vecP(l) + temp1*(m_vecYp(l) - m_matPhi(l,1));
        }
    }
    else
    {
        for (l=0;l<m_nEqn;++l)
        {
            rho = temp1*(m_vecYp(l) - m_matPhi(l,1)) - m_matPhi(l,16);
            vecY(l) = m_vecP(l) + rho;
            m_matPhi(l,15) = (vecY(l) - m_vecP(l)) - rho;
        }
    }

    m_pfDE(dX,vecY,m_vecYp,m_pAux);


    // Update differences for next step
    for (l=0;l<m_nEqn;++l)
    {
        m_matPhi(l,kp1) = m_vecYp(l) - m_matPhi(l,1);
        m_matPhi(l,kp2) = m_matPhi(l,kp1) - m_matPhi(l,kp2);
    }

    for (i=1;i<=m_nK;++i)
    {
        for (l=0;l<m_nEqn;++l)
        {
            m_matPhi(l,i) = m_matPhi(l,i) + m_matPhi(l,kp1);
        }
    }


    // Estimate error at order k+1 unless
    // - in first phase when always raise order,
    // - already decided to lower order,
    // - step size not constant so estimate unreliable
    erkp1 = 0.0;
    if ( (knew==km1) || (m_nK==12) )
    {
        m_bPhase1=false;
    }

    if (m_bPhase1)
    {
        m_nK = kp1;
        erk = erkp1;
    }
    else
    {
        if (knew==km1)
        {
            // lower order
            m_nK = km1;
            erk = erkm1;
        }
        else
        {

            if (kp1<=m_nNs)
            {
                for (l=0;l<m_nEqn;++l)
                {
                    erkp1 = erkp1 + (m_matPhi(l,kp2)/m_vecWt(l))*(m_matPhi(l,kp2)/m_vecWt(l));
                }

                erkp1 = absh*gstr[kp1]*sqrt(erkp1);

                // Using estimated error at order k+1, determine
                // appropriate order for next step
                if (m_nK>1)
                {
                    if ( erkm1<=fmin(erk,erkp1))
                    {
                        // lower order
                        m_nK=km1; erk=erkm1;
                    }
                    else
                    {
                        if ( (erkp1<erk) && (m_nK!=12) )
                        {
                            // raise order
                            m_nK=kp1; erk=erkp1;
                        }
                    }
                }
                else
                {
                    if (erkp1<0.5*erk)
                    {
                        // raise order
                        // Here erkp1 < erk < max(erkm1,ermk2) else
                        // order would have been lowered in block 2.
                        // Thus order is to be raised
                        m_nK = kp1;
                        erk = erkp1;
                    }
                }

            } // end if kp1<=ns

        } // end if knew!=km1

    } // end if !phase1


    // With new order determine appropriate step size for next step
    if ( m_bPhase1 || (p5eps>=erk*two[m_nK+1]) )
    {
        hnew = 2.0*m_dH;
    }
    else
    {
        if (p5eps<erk)
        {
            temp2 = m_nK+1;
            r = pow(p5eps/erk, 1.0/temp2);
            hnew = absh*fmax(0.5, fmin(0.9,r));
            hnew = sign(fmax(hnew, fouru*fabs(dX)), m_dH);
        }
        else hnew = m_dH;
    }

    m_dH = hnew;

    //
    // End block 4
    //
}

/// Interpolation
void CDE::Intrp (double dXout, CVector& vecYout, CVector& vecYpout )
{

    // Variables

    int     i, j, ki;
    double  eta, gamma, hi, psijm1;
    double  temp1, term;
    double  g[14], rho[14], w[14];


    g[1]   = 1.0;
    rho[1] = 1.0;

    hi = dXout - m_dX;
    ki = m_nKold + 1;

    // Initialize w[*] for computing g[*]
    for (i=1;i<=ki;++i)
    {
        temp1 = i;
        w[i] = 1.0/temp1;
    }

    // Compute g[*]
    term = 0.0;

    for (j=2;j<=ki;++j)
    {
        psijm1 = m_dPsi[j-1];
        gamma = (hi + term)/psijm1;
        eta = hi/psijm1;
        for (i=1;i<=ki+1-j;++i) w[i] = gamma*w[i] - eta*w[i+1];
        g[j] = w[1];
        rho[j] = gamma*rho[j-1];
        term = psijm1;
    }

    // Interpolate for the solution yout and for
    // the derivative of the solution ypout
    vecYpout = 0.0;
    vecYout  = 0.0;

    for (j=1;j<=ki;++j)
    {
        i = ki+1-j;
        vecYout  = vecYout  + g[i]  *m_matPhi.GetCol(i);
        vecYpout = vecYpout + rho[i]*m_matPhi.GetCol(i);
    }

    vecYout = m_vecYy + hi*vecYout;

}

/// DE integration
/// (with full control of warnings and errros status codes)
void CDE::Integ (double& dt, double dTout, CVector& vecY)
{

    // Variables

    bool    stiff, crash;           // Flags
    int     nostep;                 // Step count
    int     kle4 = 0;
    double  releps, abseps, tend;
    double  absdel, del, eps;


    // Return, if output time equals input time

    if (dt==dTout) // No integration
    {
        return;
    }


    // Test for improper parameters

    eps   = fmax(m_dRelerr,m_dAbserr);

    if ( ( m_dRelerr <  0.0         ) ||      // Negative relative error bound
         ( m_dAbserr <  0.0         ) ||      // Negative absolute error bound
         ( eps       <= 0.0         ) ||      // Both error bounds are non-positive
         ( m_emState >  DE_INVPARAM ) ||      // Invalid status flag
         ( (m_emState!= DE_INIT)  &&
           (dt != m_dTold)           ) )
    {
        m_emState = DE_INVPARAM;                 // Set error code
        return;                              // Exit
    }


    // On each call set interval of integration and counter for
    // number of steps. Adjust input error tolerances to define
    // weight vector for subroutine STEP.

    del    = dTout - dt;
    absdel = fabs(del);

    tend   = dt + 100.0*del;

    if (!m_bPermitTOUT)
    {
        tend = dTout;
    }

    nostep = 0;
    kle4   = 0;
    stiff  = false;
    releps = m_dRelerr/eps;
    abseps = m_dAbserr/eps;

    if( (m_emState==DE_INIT) || (!m_bOldPermit) || (m_dDelsgn*del<=0.0) )
    {
        // On start and restart also set the work variables x and yy(*),
        // store the direction of integration and initialize the step size
        m_bStart  = true;
        m_dX      = dt;
        m_vecYy   = vecY;
        m_dDelsgn = sign(1.0, del);
        m_dH      = sign( fmax(fouru*fabs(m_dX), fabs(dTout-m_dX)), dTout-m_dX );
    }

    while (true)
    {  // Start step loop

        // If already past output point, interpolate solution and return
        if (fabs(m_dX-dt) >= absdel)
        {
            Intrp (dTout, vecY, m_vecYpout);
            m_emState     = DE_DONE;          // Set return code
            dt         = dTout;             // Set independent variable
            m_dTold      = dt;                // Store independent variable
            m_bOldPermit = m_bPermitTOUT;
            return;                       // Normal exit
        }

        // If cannot go past output point and sufficiently close,
        // extrapolate and return
        if ( !m_bPermitTOUT && ( fabs(dTout-m_dX) < fouru*fabs(m_dX) ) )
        {
            m_dH = dTout - m_dX;
            m_pfDE(m_dX,m_vecYy,m_vecYp,m_pAux);              // Compute derivative yp(x)
            vecY = m_vecYy + m_dH*m_vecYp;                // Extrapolate vector from x to tout
            m_emState     = DE_DONE;          // Set return code
            dt         = dTout;             // Set independent variable
            m_dTold      = dt;                // Store independent variable
            m_bOldPermit = m_bPermitTOUT;
            return;                       // Normal exit
        }

        // Test for too much work
        if (nostep >= MAXNUM)
        {
            m_emState = DE_NUMSTEPS;          // Too many steps
            if (stiff) m_emState = DE_STIFF;  // Stiffness suspected
            vecY         = m_vecYy;               // Copy last step
            dt         = m_dX;
            m_dTold      = dt;
            m_bOldPermit = true;
            return;                       // Weak failure exit
        }

        // Limit step size, set weight vector and take a step
        m_dH  = sign(fmin(fabs(m_dH), fabs(tend-m_dX)), m_dH);

        for (int l=0; l<m_nEqn; ++l)
        {
            m_vecWt(l) = releps*fabs(m_vecYy(l)) + abseps;
        }

        Step ( m_dX, m_vecYy, eps, crash );


        // Test for too small tolerances
        if (crash)
        {
            m_emState     = DE_BADACC;
            m_dRelerr    = eps*releps;       // Modify relative and absolute
            m_dAbserr    = eps*abseps;       // accuracy requirements
            vecY         = m_vecYy;               // Copy last step
            dt         = m_dX;
            m_dTold      = dt;
            m_bOldPermit = true;
            return;                       // Weak failure exit
        }

        ++nostep;  // Count total number of steps

        // Count number of consecutive steps taken with the order of
        // the method being less or equal to four and test for stiffness
        ++kle4;
        if (m_nKold>  4)
        {
            kle4=0;
        }

        if (kle4>=50)
        {
            stiff=true;
        }

    } // End step loop
}

/// 初始化
void CDE::Init (double dT0, double dRel, double dAbs)
{
    m_dt      = dT0;
    m_dRelerr = dRel;
    m_dAbserr = dAbs;
    m_emState = DE_INIT;
}

/// DE integration with simplified state code handling
/// (skips over warnings, aborts in case of error)
void CDE::Integ (double dTout, CVector& vecY)
{
    do
    {
        Integ (m_dt,dTout,vecY);

        if ( m_emState==DE_INVPARAM )
        {
            std::cerr << "ERROR: invalid parameters in DE::Integ"
                      << std::endl;
            return;
        }
        else if ( m_emState==DE_BADACC )
        {
            std::cerr << "WARNING: Accuracy requirement not achieved in DE::Integ"
                      << std::endl;
        }
        else if ( m_emState==DE_STIFF )
        {
            std::cerr << "WARNING: Stiff problem suspected in DE::Integ"
                      << std::endl;
        }
    }
    while ( m_emState > DE_DONE );
}

/// Interpolation
void CDE::Intrp(double dTout, CVector& vecY)
{
    Intrp( dTout, vecY, m_vecYpout );     // Interpolate and discard interpolated
}

