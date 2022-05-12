#include <cmath>
#include <Math/Transit.h>

using namespace std;
using namespace Math;

CTransit::CTransit()
{

}

CTransit::~CTransit()
{
}


void polygon_isect_line(const CVector &p1, const CVector &p2, const CVector &pos,int &winding)
{
    double x1 = p1.GetX();
    double y1 = p1.GetY();
    double x2 = p2.GetX();
    double y2 = p2.GetY();
    double y = pos.GetY();
    int dir = 1;
    if (fabs(y1-y2)<1e-7)
    {
        /// 忽略了垂直线和比较接近的值
        return;
    }
    /// 计算两个点的Y值，如果在下面就交换
    else if (y2 < y1)
    {
        double x_tmp = x2; x2 = x1; x1 = x_tmp;
        double y_tmp = y2; y2 = y1; y1 = y_tmp;
        dir = -1;
    }
    /// 确保检测点在线段之内
    if (y >= y1 && y < y2)
    {
        /// 射线的参数化
        double x = x1 + ((x2 - x1) / (y2 - y1)) * (y - y1);
        /// 检测位于线段哪一方左方表示已经相交，右侧表示没有
        if (x<=pos.GetX())
        {
            winding += dir;
        }
    }
}

bool CTransit::ContainsPoint(const vector<CVector> &AllPoint, const CVector &pt)
{
    /// 采用射线法
    if (AllPoint.empty())
    {
        return false;
    }

    int winding_number = 0;
    CVector last_pt = AllPoint[0];
    CVector last_start = AllPoint[0];
    for (vector<CVector>::size_type i = 1; i < AllPoint.size(); ++i)
    {
        const CVector &e = AllPoint[i];
        polygon_isect_line(last_pt, e, pt, winding_number);
        last_pt = e;
    }

    if (last_pt != last_start)
    {
        polygon_isect_line(last_pt, last_start, pt, winding_number);
    }
    //	return (?(winding_number != 0):((winding_number % 2) != 0));
    //	return false;
    if(winding_number%2!=0)
    {
        return true;
    }
    else
    {
        return false;
    }
}
