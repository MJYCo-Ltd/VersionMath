#ifndef YTY_OPTICALALG_H
#define YTY_OPTICALALG_H
#include <vector>
using namespace std;
#include "SAT_global.h"
namespace Math{
class CVector;
}

namespace Satellite{
using namespace Math;
class ALGORITHM_EXPORT COpticalAlg
{
public:
    COpticalAlg();
    ~COpticalAlg();

    /**
     * @brief 阴影计算
     * @param vecPos    卫星位置 (ECI) 单位(m,m,m)
     * @param vecSun    太阳位置 (ECI) 单位(m,m,m)
     * @param vecCbPos  中心天体位置 (ECI) 单位(m,m,m)
     * @param dRcb      中心天体半径    单位(m)
     * @return 可视因子(本影区为0,半影区在0~1之间,阳照区为1)
     * @attention 圆锥形模型
     */
    static double Shadow(const CVector& vecPos,const CVector& vecSun,
                         const CVector& vecCbPos,const double& dRcb );
    /**
     * @brief 对流层对光(电磁波)的折射
     * @param dPa   大气压 单位 mb
     * @param dfH   空气湿度
     * @param dT    气温  单位 K (开尔文)
     * @param dElev 仰角  单位 弧度
     * @return 仰角折射修正角 单位 弧度
     */
    static double TroposphericRef(double dPa, double dfH, double dT, double dElev);

    /**
     * @brief 下行链路的光行迭代(卫星到地面站)
     * @param vSatPos      卫星位置   (ECI坐标 单位 m)
     * @param vecStaion    地面站位置 (ECI坐标 单位 m)
     * @param nPos         卫星的第几个位置需要与地面站通信
     * @param dMJd         卫星位置向量的 第一个值的开始时间
     * @param dStep        卫星每个数据的间隔步长
     * @param vecSat       进行链路通信时刻卫星的位置 (ECI坐标 单位 m)
     * @return             通信到达地面站需要的时间 单位 s
     */
    static double DownlegLightTime(const vector<CVector>& vSatPos, const CVector& vecStaion,
                                   int nPos,double dMJd, double dStep, CVector& vecSat);

    /**
     * @brief 上行链路的光行迭代(地面站到卫星)
     * @param vSatPos     卫星位置   (ECI坐标 单位 m)
     * @param vecStaion   地面站位置 (ECI坐标 单位 m)
     * @param nPos        卫星的第几个位置需要与地面站通信
     * @param dMJd        卫星位置向量的 第一个值的开始时间
     * @param dStep       卫星每个数据的间隔步长
     * @param vecSat      进行链路通信时刻卫星的位置 (ECI坐标 单位 m)
     * @return            通信到达卫星需要的时间 单位 s
     */
    static double UplegLightTime(const vector<CVector>& vSatPos, const CVector& vecStaion,
                                 int nPos, double dMJd, double dStep, CVector& vecSat);
};
}
#endif // YTY_OPTICALALG_H
