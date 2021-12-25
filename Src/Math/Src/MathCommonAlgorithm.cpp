#include <cmath>
#include <cstdio>
using namespace std;

#include <VersionMathCommon.h>
#include <Math/VecMat.h>
#include <Math/MathCommonAlgorithm.h>

using namespace Math;

/// 根据 b的正负设置 a 的正负
double Math::sign(double a, double b)
{
    return (b>=0.0) ? fabs(a) : - fabs(a);
}

/// x的小数部分
double Math::Frac(double x)
{
    return x-int(x);
}

/// x%y
double Math::Modulo(double x, double y)
{
    return (y*Frac(x/y));
}

double Math::acose(double arg)
{
   if( arg >= 1.)
      return( 0.);
   if( arg <= -1.)
      return( DPI);
   return( acos( arg));
}

double Math::asine( const double arg)
{
   if( arg >= 1.)
      return( DPI / 2);
   if( arg <= -1.)
      return( -DPI / 2.);
   return( asin( arg));
}
//------------------------------------------------------------------------------
//
// LU_BackSub
//
// Purpose:
//
//   LU Backsubstitution
//
//   Solves the set of n linear equations Ax=b. Here A is input, not as the
//   matrix A but rather as its LU decomposition, determined by the function
//   LU_Decomp. b is input as the right-hand side vector b, and returns with
//   the solution vector x. A and Indx are not modified by this function and
//   can be left in place for successive calls with different right-hand
//   sides b. This routine takes into account the posssibility that B will
//   begin with many zero elements, so it is efficient for use in matrix
//   inversions.
//
// Input/output:
//
//   A       LU decomposition of permutation of A
//   Indx    Permutation index vector
//   b       Right-hand side vector b; replaced by solution x of Ax=b on output
//
//------------------------------------------------------------------------------
void Math::LU_BackSub( CMatrix& A, CVector& Indx, CVector& b )
{

    // Constants

    const int  nRow = A.Row();

    // Local variables

    int     ii,i,ll,j;
    double  Sum;

    //
    // Start
    //

    ii = -1;                        // When ii is set to a nonegative value, it will
    // become the first nonvanishing element of B.
    for (i=0; i<nRow; ++i)          // We now do the forward substitution.
    {
        ll = (int) Indx(i);         // The only wrinkle is to unscramble the
        Sum = b(ll);                // permutation as we go.
        b(ll) = b(i);
        if (ii != -1)
        {
            for (j=ii; j<i; ++j)
            {
                Sum -= A(i,j)*b(j);
            }
        }
        else
        {
            if (Sum != 0.0)          // A nonzero element was encountered, so from
            {
                ii = i;
            }
        };                           // now on we will have to do the sums in the
        b(i) = Sum;                  // loop above.
    };

    for (i=nRow-1; i>=0; i--)        // Now we do the backsubstitution, eqn 2.3.7.
    {
        Sum=b(i);
        if (i<nRow-1)
        {
            for (j=i+1;j<nRow;++j)
            {
                Sum = Sum-A(i,j)*b(j);
            };
        };
        b(i) = Sum/A(i,i);          // Store a component of the solution vector X.
    };

}

//------------------------------------------------------------------------------
//
// LU_Decomp
//
// Purpose:
//
//   LU-Decomposition.
//
//   Given an nxn matrix A, this routine replaces it by the LU decomposition
//   of a rowwise permutation of itself. A is output, arranged as in
//   equation (2.3.14) of Press et al. (1986); Indx is an ouput vector which
//   records the row permutation effected by partial pivoting. This routine is
//   used in combination with LU_BackSub to solve linear equations or invert
//   a matrix.
//
// Input/output:
//
//   A       Square matrix; replaced by LU decomposition of permutation of A
//           on output
//   Indx    Permutation index vector
//
// Note:
//
//   Adapted from LUDCMP of Press et al. (1986).
//
//------------------------------------------------------------------------------
bool Math::LU_Decomp ( CMatrix& A, CVector& Indx )
{

    // Constants

    const int    nRow    = A.Row();
    const double tiny = 1.0e-20;       // A small number

    // Variables

    int     imax=0;
    int     i,j,k;
    double  aAmax, Sum, Dum;
    CVector  V(nRow);

    // Loop over rows to get scaling information

    for (i=0; i<nRow; ++i)
    {
        aAmax = 0.0;
        for (j=0;j<nRow;++j)
        {
            if (fabs(A(i,j)) > aAmax )
            {
                aAmax=fabs(A(i,j));
            }
        }

        if (aAmax==0.0)
        {
            // No nonzero largest element
            cerr << "ERROR: Singular matrix A in LU_Decomp";
            return(false);
        };
        V(i) = 1.0/aAmax;           // V stores the implicit scaling of each row
    };

    // Loop over columns of Crout's method

    for ( j=0; j<nRow; ++j )
    {

        if (j > 0)
        {
            for ( i=0; i<j; ++i )
            {   // This is equation 2.3.12 except for i=j
                Sum = A(i,j);
                if (i>0)
                {
                    for ( k=0; k<i; ++k )
                    {
                        Sum -= A(i,k)*A(k,j);
                    }
                    A(i,j) = Sum;
                };
            };
        };

        aAmax=0.0;                  // Initialize for the search of the largest
        // pivot element

        for ( i=j; i<nRow; ++i )
        {     // This is i=j of equation 2.3.12 and
            Sum = A(i,j);             // i=j+1..N of equation 2.3.13
            if (j > 0)
            {
                for ( k=0; k<j; ++k )
                {
                    Sum -= A(i,k)*A(k,j);
                }
                A(i,j) = Sum;
            };
            Dum = V(i)*fabs(Sum);     // Figure of merit for the pivot
            if (Dum >= aAmax)         // Is it better than the best so far ?
            {
                imax  = i;
                aAmax = Dum;
            };
        };

        if (j != imax)                // Do we need to interchange rows?
        {
            for ( k=0; k<nRow; ++k)      // Yes, do so ...
            {
                Dum = A(imax,k);
                A(imax,k) = A(j,k);
                A(j,k) = Dum;
            }
            V(imax) = V(j);           // Also interchange the scale factor
        };

        Indx(j) = imax;

        if (j != nRow-1)                 // Now finally devide by the pivot element
        {
            if (A(j,j) == 0.0)        // If the pivot element is zero the matrix
            {
                A(j,j) = tiny;        // is singular (at least to the precision of
            };                        // the algorithm). For some applications on
            Dum=1.0/A(j,j);           // singular matrices, it is desirable to
            for (i=j+1;i<nRow;++i)
            {                         // substitude tiny for zero.
                A(i,j)=A(i,j)*Dum;
            };
        };

    };   // Go back for the next column in the reduction

    if (A(nRow-1,nRow-1)==0.0)
    {
        A(nRow-1,nRow-1)=tiny;
    }

    return (true);
}

/// 将值转换到 [0~2pi]
double Math::zero_to_two_pi( double angle_in_radians)
{
    angle_in_radians = fmod( angle_in_radians, D2PI);
    if( angle_in_radians < 0.)
    {
        angle_in_radians += D2PI;
    }
    return( angle_in_radians);
}

/// 归一化勒让德多项式
void Math::Legendre_sphPl(const int LL,const double x,double PL[])
{
    double W1,W2,W3,W4;
    PL[0] = 1;
    PL[1] = sqrt(3.0) * x;
//    PL[2] = sqrt(5.0) * (1.5*x*x-0.5);
    if(LL<2) return;
    for(int L=2;L<=LL;++L)
    {
        W1 = 1.0/L;
        W2 = 2.0*L;
        W3 = sqrt((W2-1.0)/(W2-3.0)); // sqrt(2l-1)/(2L-3)
        W4 = sqrt((W2+1.0)/(W2-1.0)); // sqrt(2l+1)/(2L-1)
        PL[L] = W4*( (2.0-W1)*x*PL[L-1]-W3*(1.0-W1)*PL[L-2] );
    }
}

/// 归一化勒让德多项式
void Math::Legendre_sphPlm(const int LL,const double x,double PLM[][71])
{
    double G3 = sqrt(3.0);
    double GX2 = sqrt(1-x*x);
    double W1,W2,WL1,WL2,WL3,WL0,WM1,WM2,WM3,WM12,WLM1,WLM2;
    PLM[1][1] = G3*GX2;
    for(int L=2;L<=LL;++L)
    {
        W1 = 2.0*L;         // 2L
        W2 = sqrt(W1);      // sqrt(2L)
        WL1 = sqrt(W1+1.0); // sqrt(2L+1)
        WL2 = sqrt(W1-1.0); // sqrt(2L-1)
        WL3 = sqrt(W1-3.0); // sqrt(2L-3)
        WL0 = WL1/W2;       // sqrt(2L+1)/2L
        PLM[L][L] = WL0 * GX2 * PLM[L-1][L-1];
        for(int M=1;M<=L-1;++M)
        {
            WM1 = L + M;       // L+M
            WM2 = L - M;       // L-M
            WM3 = sqrt((WM1-1.0)*(WM2-1.0)); // sqrt(L+M-1)*(L-M-1)
            WM12 = sqrt(WM1*WM2);    // sqrt(L+M)*(L-M)
            WLM1 = WL1*WL2/WM12;     // sqrt(2L+1)*(2L-1)/(L+M)/(L-M)
            WLM2 = WM3*WL1/(WM12*WL3); // sqrt(2L+1)*(L+M-1)*(L-M-1)/(2L-3)/(L+M)/(L-M)
            PLM[L][M] = WLM1 * x * PLM[L-1][M] - WLM2 * PLM[L-2][M];
        }
    }
}

/// 用递推方法计算sin(m*x)和cos(m*x)
void Math::SmxCmx(const int LL,const double S1X,const double C1X,double* SX,double* CX)
{
    SX[0] = 0;    // sin(0)
    CX[0] = 1;    // cos(0)
    SX[1] = S1X;  // sin(x)
    CX[1] = C1X;  // cos(x)
    SX[2] = 2.0*SX[1]*CX[1]; // 2sin(x)cos(x)
    CX[2] = 2.0*CX[1]*CX[1] - 1.0; // 2cos(x)^2 -1
    if(LL<3) return;
    double CX2 = 2.0 * CX[1]; // 2cos(x)
    for(int MM=3;MM<=LL;++MM)
    {
        SX[MM] = CX2*SX[MM-1] - SX[MM-2]; // 2cos(x)sin[(m-1)x] - sin[(m-2)x]
        CX[MM] = CX2*CX[MM-1] - CX[MM-2]; // 2cos(x)cos[(m-1)x] - cos[(m-2)x]
    }
}
