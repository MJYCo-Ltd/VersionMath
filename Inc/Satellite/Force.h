#ifndef YTY_FORCE_H
#define YTY_FORCE_H
/*****************************************
  作用：本类进行卫星受力情况的计算，力学模型有
       光压、大气阻力、非球形引力
  注意：此类是非线程安全的，在多线程中使用需要对变量
       进行加锁保护，否则运算结果不能保证
采用的模型： 大气模型 参见CAtmosphere
           引力模型 参见CSTKGraveModel
 *****************************************/
#include "SAT_global.h"
#include "Atmosphere.h"

namespace Math {
class CVector;
class CMatrix;
}

namespace Physical{
using namespace Math;
class ALGORITHM_EXPORT CForce
{
public:
    CForce();
    ~CForce();

    /**
     * @brief 地球引力(如果有其他星球的重力场系数，可以替换第三体引力)
     * @param vecR    卫星在惯性系下的位置 (x,y,z)[m]
     * @param matT    ECI到ECF的旋转矩阵
     * @param dGM     地球引力系数        [m^3/s^2]
     * @param dR      地球参考半径        [m]
     * @param matCS   未归一的重力势系数
     * @param nN_Max  阶数
     * @param nM_Max  次数
     * @return 摄动引起的加速度          [m/s^2]
     * @attention nM_Max<=nN_Max
     *            nM_Max 和 nN_Max 小于等于 matCS 最大的阶次
     */
    static CVector AccelHarmonic (const CVector& vecR, const CMatrix& matT,
                                  double dGM, double dR, const CMatrix& matCS,
                                  int nN_Max, int nM_Max );

    /**
     * @brief 带谐项摄动
     * @param vecR   卫星在惯性系下的位置 (x,y,z)[m]
     * @param matT   ECI到ECF的旋转矩阵
     * @param dGM    地球引力系数        [m^3/s^2]
     * @param dR     地球参考半径        [m]
     * @param matCS  归一化的重力势系数
     * @param nN_Max 阶数
     * @return 摄动引起的加速度         [m/s^2]
     */
    static CVector ZonalHarmonic(const CVector& vecR, const CMatrix& matT,
                                 double dGM, double dR, const CMatrix& matCS,
                                 int nN_Max);

    /**
     * @brief 田谐项摄动
     * @param vecR   卫星在惯性系下的位置 (x,y,z)[m]
     * @param matT   ECI到ECF的旋转矩阵
     * @param dGM    地球引力系数        [m^3/s^2]
     * @param dR     地球参考半径        [m]
     * @param matCS  归一化的重力势系数
     * @param nN_Max 阶数
     * @return 摄动引起的加速度           [m/s^2]
     */
    static CVector TesseralHarmonic(const CVector& vecR, const CMatrix& matT,
                                    double dGM, double dR, const CMatrix& matCS,
                                    int nN_Max);
    /**
     * @brief 太阳光压(默认认为太阳光与太阳能帆板垂直)
     * @param vecR   卫星在惯性系下的位置 (x,y,z)[m]
     * @param vecRSun太阳在惯性系下的位置 (x,y,z)[m]
     * @param dArea  卫星面积           [m^2]
     * @param dMass  卫星质量           [kg]
     * @param dCR    太阳光压系数(由卫星材质决定)
     * @param dP0    地球受太阳光辐射压   [N/m^-2]
     * @param dAU    地球到太阳平均距离   [m]
     * @return 摄动引起的加速度 单位 m/s^2
     */
    static CVector AccelSolrad(const CVector& vecR, const CVector& vecRSun,
                               double dArea, double dMass, double dCR,
                               double dP0, double dAU );

    /**
     * @brief 大气阻力
     * @param vecRSun 太阳在惯性系下的位置 (x,y,z)[m]
     * @param vecR    卫星在惯性系下的位置 (x,y,z)[m]
     * @param vecV    卫星在惯性系下的速度 (vx,vy,vz)[m/s]
     * @param matT    ECI到ECF的旋转乘以地球旋转的逆矩阵得到的矩阵
     * @param dArea   卫星面积           [m^2]
     * @param dMass   卫星质量           [kg]
     * @param dCD     阻尼系数(由卫星材质决定)
     * @param typeAtmos 大气模型
     * @return 摄动引起的加速度 [m/s^2]
     */
    static CVector AccelDrag(const CVector& vecRSun, const CVector& vecR,
                             const CVector& vecV,const CMatrix& matT, double dMJD,
                             double dArea, double dMass, double dCD, AtmosModel typeAtmos);

    /**
     * @brief 第三体引力(将第三体作为质点考虑)
     * @param vecR      卫星在惯性系下的位置   (x,y,z)[m]
     * @param vecPoint  第三体在惯性系下的位置 (x,y,z)[m]
     * @param dGM       第三体引力系数        [m^3/s^2]
     * @return 摄动引起的加速度                 [m/s^2]
     */
    static CVector AccelPointMass(const CVector& vecR, const CVector& vecPoint, double dGM);

    /**
     * @brief 后牛顿效应
     * @param vecR     卫星在惯性系下的位置   (x,y,z)[m]
     * @param vecV     卫星在惯性系下的速度(vx,vy,vz)[m/s]
     * @param dGM      第三体引力系数               [m^3/s^2]
     * @return 摄动引起的加速度                     [m/s^2]
     */
    static CVector AccelPostNewton(const CVector& vecR, const CVector& vecV, double dGM);
};
}
#endif // YTY_FORCE_H
