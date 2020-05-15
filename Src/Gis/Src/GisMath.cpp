#include "GisMath.h"
#include "sofam.h"
#include "Vector.h"
#include <cmath>

using namespace std;
using namespace Math;

CGisMath::CGisMath(double dRequator,double dRPolar)
{
    SetEllipsePara(dRequator,dRPolar);
}

/// 设置椭球参数
void CGisMath::SetEllipsePara(double dREquator, double dRPolar)
{
    m_dREquator = dREquator;
    m_dRPolar   = dRPolar;
    Init();
}

double CGisMath::CalArea(const vector<CVector>& v2DPoint, bool bGeo)
{
    int nCount = v2DPoint.size();
    double dMtotalArea = 0;

    if (nCount > 2)
    {

        if(2 == nCount && fabs(v2DPoint[nCount-1](0) - v2DPoint[0](0)) < 1e-9
                      && fabs(v2DPoint[nCount-1](1) - v2DPoint[0](1)) < 1e-9)
        {
            return(dMtotalArea);
        }

        /// 经纬度坐标下的球面多边形
        if (bGeo)
        {
            double LowX = 0.0;
            double LowY = 0.0;
            double MiddleX = 0.0;
            double MiddleY = 0.0;
            double HighX = 0.0;
            double HighY = 0.0;

            double AM = 0.0;
            double BM = 0.0;
            double CM = 0.0;

            double AL = 0.0;
            double BL = 0.0;
            double CL = 0.0;

            double AH = 0.0;
            double BH = 0.0;
            double CH = 0.0;

            double CoefficientL = 0.0;
            double CoefficientH = 0.0;

            double ALtangent = 0.0;
            double BLtangent = 0.0;
            double CLtangent = 0.0;

            double AHtangent = 0.0;
            double BHtangent = 0.0;
            double CHtangent = 0.0;

            double ANormalLine = 0.0;
            double BNormalLine = 0.0;
            double CNormalLine = 0.0;

            double OrientationValue = 0.0;

            double AngleCos = 0.0;

            double Sum1 = 0.0;
            double Sum2 = 0.0;
            double Count2 = 0;
            double Count1 = 0;


            double Sum = 0.0;
            double Radius = 6378000;

            for (int i = 0; i < nCount; i++)
            {
                if (i == 0)
                {
                    LowX = v2DPoint[nCount - 1](0);
                    LowY = v2DPoint[nCount - 1](1);
                    MiddleX = v2DPoint[0](0);
                    MiddleY = v2DPoint[0](1);
                    HighX = v2DPoint[1](0);
                    HighY = v2DPoint[1](0);
                }
                else if (i == nCount - 1)
                {
                    LowX = v2DPoint[nCount - 2](0);
                    LowY = v2DPoint[nCount - 2](1);
                    MiddleX = v2DPoint[nCount - 1](0);
                    MiddleY = v2DPoint[nCount - 1](1);
                    HighX = v2DPoint[0](0);
                    HighY = v2DPoint[0](1);
                }
                else
                {
                    LowX = v2DPoint[i - 1](0);
                    LowY = v2DPoint[i - 1](1);
                    MiddleX = v2DPoint[i](0);
                    MiddleY = v2DPoint[i](1);
                    HighX = v2DPoint[i + 1](0);
                    HighY = v2DPoint[i + 1](1);
                }

                AM = cos(MiddleY) * cos(MiddleX);
                BM = cos(MiddleY) * sin(MiddleX);
                CM = sin(MiddleY);
                AL = cos(LowY) * cos(LowX);
                BL = cos(LowY) * sin(LowX);
                CL = sin(LowY);
                AH = cos(HighY) * cos(HighX);
                BH = cos(HighY) * sin(HighX);
                CH = sin(HighY);


                CoefficientL = (AM * AM + BM * BM + CM * CM) / (AM * AL + BM * BL + CM * CL);
                CoefficientH = (AM * AM + BM * BM + CM * CM) / (AM * AH + BM * BH + CM * CH);

                ALtangent = CoefficientL * AL - AM;
                BLtangent = CoefficientL * BL - BM;
                CLtangent = CoefficientL * CL - CM;
                AHtangent = CoefficientH * AH - AM;
                BHtangent = CoefficientH * BH - BM;
                CHtangent = CoefficientH * CH - CM;


                AngleCos = (AHtangent * ALtangent + BHtangent * BLtangent + CHtangent * CLtangent) / (sqrt(AHtangent * AHtangent + BHtangent * BHtangent + CHtangent * CHtangent) * sqrt(ALtangent * ALtangent + BLtangent * BLtangent + CLtangent * CLtangent));

                AngleCos = acos(AngleCos);

                ANormalLine = BHtangent * CLtangent - CHtangent * BLtangent;
                BNormalLine = 0 - (AHtangent * CLtangent - CHtangent * ALtangent);
                CNormalLine = AHtangent * BLtangent - BHtangent * ALtangent;

                if (AM != 0)
                    OrientationValue = ANormalLine / AM;
                else if (BM != 0)
                    OrientationValue = BNormalLine / BM;
                else
                    OrientationValue = CNormalLine / CM;

                if (OrientationValue > 0)
                {
                    Sum1 += AngleCos;
                    Count1++;

                }
                else
                {
                    Sum2 += AngleCos;
                    Count2++;
                    //Sum +=2*Math.PI-AngleCos;
                }

            }

            if (Sum1 > Sum2)
            {
                Sum = Sum1 + (D2PI * Count2 - Sum2);
            }
            else
            {
                Sum = (D2PI * Count1 - Sum1) + Sum2;
            }

            //平方米
            dMtotalArea = (Sum - (nCount - 2) * DPI) * Radius * Radius;
        }
        else
        { //非经纬度坐标下的平面多边形

            int i, j;
            //double j;
            double p1x, p1y;
            double p2x, p2y;
            for ( i = nCount - 1, j = 0; j < nCount; i = j, j++)
            {

                p1x = v2DPoint[i](0);
                p1y = v2DPoint[i](1);

                p2x = v2DPoint[j](0);
                p2y = v2DPoint[j](1);

                dMtotalArea += p1x * p2y - p2x * p1y;
            }
            dMtotalArea /= 2.0;
        }
        return dMtotalArea;
    }
    return 0;
}

double CGisMath::ComputePolygonArea(const vector<CVector> &v2DPoint, bool bGeo)
{
    double x1, y1, dx, dy;
    double Qbar1, Qbar2;

    int nCount = v2DPoint.size();

    double area = 0.0;

    if (nCount < 3)
    {
        return(area);
    }



    double x2 = v2DPoint[nCount-1](0);
    double y2 = v2DPoint[nCount-1](1);
    Qbar2 = GetQbar( y2 );

    for ( int i = 0; i < nCount; i++ )
    {
        x1 = x2;
        y1 = y2;
        Qbar1 = Qbar2;

        x2 = v2DPoint[i](0);
        y2 = v2DPoint[i](1);
        Qbar2 = GetQbar( y2 );

        if ( x1 > x2 )
            while ( x1 - x2 > DPI )
                x2 += D2PI;
        else if ( x2 > x1 )
            while ( x2 - x1 > DPI )
                x1 += D2PI;

        dx = x2 - x1;
        area += dx * ( m_Qp - GetQ( y2 ) );

        if (( dy = y2 - y1 ) != 0.0 )
            area += dx * GetQ( y2 ) - ( dx / dy ) * ( Qbar2 - Qbar1 );
    }
    if (( area *= m_AE ) < 0.0 )
        area = -area;

    if ( area > m_E )
        area = m_E;
    if ( area > m_E / 2 )
        area = m_E - area;

    return (area);
}

/// 初始化数据
void CGisMath::Init()
{
    double a2 = ( m_dREquator * m_dREquator );
    double e2 = 1 - ( a2 / ( m_dRPolar * m_dRPolar ) );
    double e4, e6;

    e4 = e2 * e2;
    e6 = e4 * e2;

    m_AE = a2 * ( 1 - e2 );

    m_QA = ( 2.0 / 3.0 ) * e2;
    m_QB = ( 3.0 / 5.0 ) * e4;
    m_QC = ( 4.0 / 7.0 ) * e6;

    m_QbarA = -1.0 - ( 2.0 / 3.0 ) * e2 - ( 3.0 / 5.0 ) * e4  - ( 4.0 / 7.0 ) * e6;
    m_QbarB = ( 2.0 / 9.0 ) * e2 + ( 2.0 / 5.0 ) * e4  + ( 4.0 / 7.0 ) * e6;
    m_QbarC =                     - ( 3.0 / 25.0 ) * e4 - ( 12.0 / 35.0 ) * e6;
    m_QbarD = ( 4.0 / 49.0 ) * e6;

    m_Qp = GetQ( DPI / 2 );
    m_E  = 4 * DPI * m_Qp * m_AE;
    if ( m_E < 0.0 )
    m_E = -m_E;
}

double CGisMath::GetQ(double dX )
{
    double sinx, sinx2;

    sinx = sin( dX );
    sinx2 = sinx * sinx;

    return sinx *( 1 + sinx2 *( m_QA + sinx2 *( m_QB + sinx2 * m_QC ) ) );
}


double CGisMath::GetQbar(double dX )
{
    double cosx, cosx2;

    cosx = cos( dX );
    cosx2 = cosx * cosx;

    return cosx *( m_QbarA + cosx2 *( m_QbarB + cosx2 *( m_QbarC + cosx2 * m_QbarD ) ) );
}
