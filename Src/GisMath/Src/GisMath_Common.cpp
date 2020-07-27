#include "GisMath_Common.h"

/*****************************************************
 * 东北天坐标系 笛卡尔坐标系，以物体的重心为坐标原点，物体的正东方向为X轴正方向，
 *            物体重心与地心连线的反方向(天)为Z轴正方向，
 *            根据右手定则，确定Y轴指向(指向正北)。
 *            该坐标系又称站心坐标系，因该坐标常用于地面站观测物体，因此原点又称站心。
 *****************************************************/

/// 检查是否初始化
void CheckInit()
{
    if(!S_INIT)
    {
        GisMath::InitGis(GisMath::WGS_84);
    }
}

/// 初始化地理信息系统
void GisMath::InitGis(ELLIPSOID typeEllipsoid)
{
    switch (typeEllipsoid)
    {
    case BJ_54:
    case WGS_72:
    case GRS_80:
    case CGCS_2000:
        iauEform(typeEllipsoid,&EARTH_DA,&EARTH_DF);
        EM_TYPE = typeEllipsoid;
        break;
    default:
        iauEform(WGS84,&EARTH_DA,&EARTH_DF);
        EM_TYPE = WGS_84;
        break;
    }

    EARTH_DB = EARTH_DA*(1-EARTH_DF);
    /// 初始化系数
    geod_init(&PJ_GEOD,EARTH_DA,EARTH_DF);
    S_INIT = true;
}

/// 求射线与椭球的交点
bool GisMath::CalLineInterEllipsoid(const CVector &pt, const CVector &stDir, CVector &rInsertPos)
{
    CheckInit();
    /// 椭球参数
    double a  = EARTH_DA;
    double b  = EARTH_DB;
    CVector vDir(stDir);
    vDir.Normalize();

    double a2 = a*a;
    double b2 = b*b;
    double dx = vDir.GetX();
    double dy = vDir.GetY();
    double dz = vDir.GetZ();
    double ex = pt.GetX();
    double ey = pt.GetY();
    double ez = pt.GetZ();

    double A   = b2*dx*dx + b2*dy*dy + a2*dz*dz;
    double B   = 2 * (b2*dx*ex + b2*dy*ey+ a2*dz*ez);
    double C   = b2*ex*ex + b2*ey*ey + a2*ez*ez - a2*b2;

    /// 求该视线与椭球的两个交点
    double delta = B*B-4*A*C;
    if(delta<0) return false; //没有交点

    double deltaS = sqrt(delta);
    double t1 = (-B + deltaS)/(2*A);
    double t2 = (-B - deltaS)/(2*A);

    /// 如果t1,t2都小于0，则认为视线与椭球没有交点
    if(t1<0 && t2<0)
    {
        return(false);
    }

    /// 取t1, t2中绝对值较小的那个
    double t = fabs(t1)<fabs(t2)? t1:t2;

    ///  根据直线参数方程：X = Dt + E，可得直线与椭球的交点坐标
    rInsertPos(0) = t * dx + ex;
    rInsertPos(1) = t * dy + ey;
    rInsertPos(2) = t * dz + ez;

    return(true);
}
