#ifndef YTY_KEPLER_H
#define YTY_KEPLER_H

#include <Satellite/SAT_global.h>

#include <Math/VecMat.h>
#include <Math/Quaternion.h>

namespace Satellite{

class ALGORITHM_EXPORT CKepler
{
public:
    CKepler(double dA, double dE, double dI, double dRAAN, double dMA, double dAP);
    ~CKepler();

    const Math::CVector& CalPV(double dT);
    /**
     * @brief 解开普勒方程
     * @param dGM  中心天体引力系数   [m^3/s^2]
     * @param vKep 开普勒六根数    (半长轴,离心率,轨道倾角,升交点赤经,近地点幅角,平近点角)
     *                             [m, ,rad,rad,rad,rad]
     * @param dT   相对于历元的时间       [s]
     * @return 位置速度(x,y,z,vx,vy,vz)  [m,m/s][ECI]
     */
    static Math::CVector State(const double& dGM, const Math::CVector& vKep, double dT=0.0);

    /**
     * @brief 通过位置速度解算开普勒六根数
     * @param dGM    中心天体引力系数   [m^3/s^2]
     * @param vY     位置速度(x,y,z,vx,vy,vz) [m,m/s][ECI]
     * @return 开普勒六根数(半长轴,离心率,轨道倾角,升交点赤经,近地点幅角,平近点角)
     *                     [m, ,rad,rad,rad,rad]
     * @attention 该方法仅适用于椭圆轨道
     */
    static Math::CVector Elements(const double& dGM, const Math::CVector& vY );

    /**
     * @brief 通过位置速度解算开普勒六根数
     * @param dGM    中心天体引力系数          [m^3/s^2]
     * @param vY     位置速度(x,y,z,vx,vy,vz) [m,m/s][ECI]
     * @return 开普勒六根数(半长轴,离心率,轨道倾角,升交点赤经,近地点幅角,平近点角)
     *                    [m, ,rad,rad,rad,rad]
     */
    static Math::CVector ClassicalElements(const double& dGM, const Math::CVector& vY);

    /**
     * @brief 通过两组位置解算开普勒六根数
     * @param dGM   中心天体引力系数   [m^3/s^2]
     * @param dMjd1 位置1的时间       [MJD]
     * @param dMjd2 位置2的时间       [MJD]
     * @param vPos1 位置1(x,y,z)     [m][ECI]
     * @param vPos2 位置2(x,y,z)     [m][ECI]
     * @return dMjd1时的开普勒六根数(半长轴,离心率,轨道倾角,升交点赤经,近地点幅角,平近点角)
     *                             [m, ,rad,rad,rad,rad]
     */
    static Math::CVector Elements(const double& dGM, const double& dMjd1, const double& dMjd2,
                            const Math::CVector& vPos1, const Math::CVector& vPos2);

    /**
     * @brief 将经典六根数转化成 TLE数据
     * @param vKep     dJDEpoch 对应的开普勒六根数(半长轴,离心率,轨道倾角,升交点赤经,近地点幅角,平近点角)
     *                                         [m, ,rad,rad,rad,rad]
     * @param dJDEpoch 历元时间        [JD]
     * @param nNorad   卫星的NORAD编号
     * @param str      存储两行根数的数组
     * @return 成功 返回true 失败返回 false
     */
    static bool Classical2TLE(const Math::CVector& vKep, const double& dJDEpoch,
                              int nNorad, char str[2][73]);

    /**
     * @brief 根据两个瞬时根数计算 相对位置速度
     * @param vKepChaser 追踪根数 (半长轴,离心率,轨道倾角,升交点赤经,近地点幅角,平近点角)
     *                               [m, ,rad,rad,rad,rad]
     * @param vKepTarget 目标根数 (半长轴,离心率,轨道倾角,升交点赤经,近地点幅角,平近点角)
     *                              [m, ,rad,rad,rad,rad]
     * @param dGM        中心天体引力系数 [m^3/s^2]
     * @return 追踪相对于目标的位置、速度 (x,y,z,vx,vy,vz)[m,m/s]
     */
    static Math::CVector RIC(const Math::CVector& vKepChaser,const Math::CVector& vKepTarget,
                       const double& dGM);

    /**
     * @brief 根据目标轨道的瞬时根数 以及另一轨道相对于目标星的位置速度 计算轨道根数
     * @param vKepTarget 目标轨道瞬时根数 (半长轴,离心率,轨道倾角,升交点赤经,近地点幅角,平近点角)
     *                                   [m, ,rad,rad,rad,rad]
     * @param pv         相对于目标卫星的位置速度(x,y,z,vx,vy,vz) [m,m/s]
     * @param dGM        中心天体引力系数 [m^3/s^2]
     * @return 轨道瞬时根数 (半长轴,离心率,轨道倾角,升交点赤经,近地点幅角,平近点角)
     *                      [m,  ,rad,rad,rad,rad]
     */
    static Math::CVector CIR(const Math::CVector& vKepTarget,const Math::CVector& vPV,
                       const double& dGM);
    ///////////////////根数转换////////////////////////////////
    /**
     * 开普勒六根数(半长轴,离心率,轨道倾角,升交点赤经,近地点幅角,平近点角)[m,,rad,rad,rad,rad]
     */

    /**
     * @brief 将瞬时根数转换成平根数
     * @param rIKep 瞬时根数 (半长轴,离心率,轨道倾角,升交点赤经,近地点幅角,平近点角)
     *                        [m, ,rad,rad,rad,rad]
     * @return 平根数        (半长轴,离心率,轨道倾角,升交点赤经,近地点幅角,平近点角)
     *                         [m, ,rad,rad,rad,rad]
     */
    static Math::CVector Instant2Mean(const Math::CVector &rIKep);

    /**
     * @brief 将平根数转换成瞬时根数
     * @param rMKep 平根数  (半长轴,离心率,轨道倾角,升交点赤经,近地点幅角,平近点角)
     *                        [m, ,rad,rad,rad,rad]
     * @return  瞬时根数     (半长轴,离心率,轨道倾角,升交点赤经,近地点幅角,平近点角)
     *                         [m, ,rad,rad,rad,rad]
     */
    static Math::CVector Mean2Instant(const Math::CVector &rMKep);
    ///////////////////根数转换end////////////////////////////////

    /**
     * @brief 计算真近点角
     * @param de 轨道离心率 [0,1]
     * @param dM 平近点角   [Rad]
     * @return 真近点角     [Rad]
     */
    static double TrueAnom(double de, double dM);

    /**
     * @brief 计算平近点角
     * @param de 轨道离心率 [0,1]
     * @param dT 真近点角   [Rad]
     * @return 平近点角     [Rad]
     */
    static double MeanAnom(double da,double de, double dT);

    /**
     * @brief 计算偏近点角
     * @param dM 平近点角 单位 [rad]
     * @param de 轨道偏心率    [0,1]
     * @return 偏近点角 单位   [rad]
     */
    static double EccAnom(double dM, double de);

    /**
     * @brief 计算卫星的周期
     * @param da  轨道长半轴        [m]
     * @param dGM 中心天体引力系数   [m^3/s^2]
     * @return 卫星运行一个周期的时间 [s]
     * @attention da与dGM的单位必须统一
     */
    static double T(double da, double dGM);

    /**
     * @brief 远地点
     * @param da 长半轴       [m]
     * @param de 离心率       [0~1]
     * @return 远地点据地心距离 [m]
     */
    static double Apogee(double da, double de);

    /**
     * @brief 近地点
     * @param da 长半轴      [m]
     * @param de 离心率      [0~1]
     * @return 近地点地心距离 [m]
     */
    static double Perigee(double da, double de);
private:
    double m_dA,m_dE,m_dI,m_dRAAN,m_dMA,m_dAP,m_dN,m_dFac,m_dV;
    Math::CQuaternion m_quatPQW;
    Math::CVector m_vPV;
};
}
#endif // YTY_KEPLER_H
