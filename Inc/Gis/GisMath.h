#ifndef GISMATH_H
#define GISMATH_H
#include <vector>
#include "GIS_global.h"
using namespace std;
namespace Math{
class CVector;
}

using namespace Math;
class GISSHARED_EXPORT CGisMath
{

public:
    /**
     * @brief 默认采用
     * @param 赤道半径 为 6378140.0 m
     * @param 极半径   为 6356755.0 m
     */
    CGisMath(double dRequator=6378140.0,double dRPolar=6356755.0);

    /**
     * @brief 设置椭球参数
     * @param dREquator 赤道半径  [m]
     * @param dRPolar   极半径    [m]
     */
    void SetEllipsePara(double dREquator,double dRPolar);

    /**
     * @brief 计算覆盖面积
     * @param v2DPoint 点数据
     *                 1、 经度、纬度 [rad、rad]
     *                 2、 x位置、y位置 [m、m]
     * @param bGeo     是否是经纬度
     * @return 计算出来的面积   m^2
     */
    static double CalArea(const vector<CVector>& v2DPoint, bool bGeo=true);

    /**
     * @brief 计算球面多边形面积
     * @param v2DPoint 点数据
     *                 1、 经度、纬度 [rad、rad]
     *                 2、 x位置、y位置 [m、m]
     * @param bGeo     是否是经纬度
     * @return 计算出来的面积   m^2
     */
    double ComputePolygonArea(const vector<CVector>& v2DPoint, bool bGeo=true);
private:
    /**
     * @brief 根据椭球参数初始化数据
     */
    void Init();

    double GetQ( double dX );

    double GetQbar( double dX );

private:
    double  m_dREquator;    /// 赤道半径  [m]
    double  m_dRPolar;      /// 极半径    [m]

    double m_QA, m_QB, m_QC;
    double m_QbarA, m_QbarB, m_QbarC, m_QbarD;
    double m_AE;  /* a^2(1-e^2) */
    double m_Qp;
    double m_E;
};

#endif // GISMATH_H
