#ifndef YTY_COORSYS_H
#define YTY_COORSYS_H

/*****************************************
  作用：完成各种坐标系的相互转换
       本类采用的ECI为ICRF；ECF为ITRF
  注意：此类是非线程安全的，在多线程中使用需要对变量
       进行加锁保护，否则运算结果不能保证
名词解释：
                           地固系[ECF]
        GTRS(Geocentric Terrestrial Reference System)    地心地球参考系
        ITRS(International Terrestrial Reference System) 国际地球参考系
        ITRF(International Terrestrial Reference Frame)  国际地球参考架
        VVLH()                                           轨道参考系

                           地惯系[ECI]
        GCRS(Geocentric Celestial Reference System)      地心天球参考系
        ICRS(International Celestial Reference System)   国际天球参考系
        ICRF(International Celestial Reference Frame)    国际天球参考架
        CIRS(Celestial Intermediate Reference System)    天球中间参考系
        J2000(Mean Equator and Mean Equinox at 2000)     平赤道平春分点坐标系
        TEME(True Equator and Mean Equinox)              真赤道平春分点坐标系

                         其他坐标系
        ENU(East-North-Zenith coordinate System)         东北天坐标系(站心坐标系)
 *****************************************/
#include "SAT_global.h"

namespace Math{
class CVector;
class CMatrix;
}

namespace Aerospace{
using namespace Math;

class ALGORITHM_EXPORT CCoorSys
{
public:
    /************************ 不同层次的转换 ************************/
    /////////////////// ICRS 与 ITRS 坐标转换///////////////////////////
    /// @brief  地惯系 到 地固系的坐标转换,考虑了岁差、章动、极移
    /// @param  dMJD   约简儒略日
    /// @param  vEci   地惯系下的坐标 (适合数据为3的倍数)
    /// @param  VEcf   地固系下的坐标 (同上)
    /// @return 成功返回true 失败返回false
    /// @attention 支持批量处理
    static bool ECI2ECF(const double& dMJD, const CVector &vEci,CVector &vEcf);
    static bool ECF2ECI(const double& dMJD, const CVector& vEcf, CVector& vEci);
    /// @return 地惯系 到 地固系 的旋转矩阵
    static CMatrix ECI2ECF(const double& dMJD);
    static CMatrix ECF2ECI(const double& dMJD);

    /////////////////// TEME 与 ICRS 坐标转换///////////////////////////
    /// @brief  TEME到ICRS的坐标转换
    /// @param  dJD     约简儒略日
    /// @param  VTeme   VTeme下的坐标 (同 ECI2ECF)
    /// @param  vEci    地惯系下的坐标 (同上)
    /// @return 成功返回true 失败返回false
    /// @attention 支持批量处理
    static bool TEME2ECI(const double& dMJD, const CVector &vTeme, CVector &vEci);
    static bool ECI2TEME(const double& dMJD, const CVector &vEci, CVector &vTeme);
    /// @return TEME 转 ICRS 的旋转矩阵
    static CMatrix TEME2ECI(const double& dMJD);
    static CMatrix ECI2TEME(const double& dMJD);

    ////////////////// TEME 与 ECF 坐标转换 ///////////////////////////
    /// @brief TEME 到 ECF的坐标转换
    /// @param dMJD
    /// @param vTeme
    /// @param vEcf
    /// @return
    /// @attention 支持批处理
    static bool TEME2ECF(const double& dMJD, const CVector &vTeme, CVector &vEcf);
    static bool ECF2TEME(const double& dMJD, const CVector &vEcf, CVector &vTeme);
    /// @return TEME 转 ECF 的旋转矩阵
    static CMatrix TEME2ECF(const double& dMJD);
    static CMatrix ECF2TEME(const double& dMJD);

    /////////////////// J2000 与 ITRS 坐标转换///////////////////////////
    /// @brief  J2000 到 地固系的坐标转换,考虑了岁差、章动、极移
    /// @param  dMJD   约简儒略日
    /// @param  vJ2000 J2000下的坐标 (适合数据为3的倍数)
    /// @param  VEcf   地固系下的坐标 (同上)
    /// @return 成功返回true 失败返回false
    /// @attention 支持批量处理
    static bool J20002ECF(const double& dMJD, const CVector& vJ2000, CVector& vEcf);
    static bool ECF2J2000(const double& dMJD, const CVector& vEcf, CVector& vJ2000);
    /// @return J2000坐标系 到 地固系 的旋转矩阵
    static CMatrix J20002ECF(const double& dMJD);
    static CMatrix ECF2J2000(const double& dMJD);

    /////////////////// VVLH 与 ECI 坐标转换///////////////////////////
    /// @brief VVLH2ECI
    /// @param vKep
    /// @param vVVLH
    /// @param vECI
    /// @return 成功返回true 失败返回false
    ///
    static bool VVLH2ECI(const CVector &vKep, const CVector& vVVLH, CVector& vECI);
    static bool ECI2VVLH(const CVector &vKep, const CVector& vECI, CVector& vVVLH);
    /// ECI到轨道坐标系的旋转矩阵
    static CMatrix VVLH2ECI(const CVector& vKep);
    static CMatrix ECI2VVLH(const CVector& vKep);
private:
    CCoorSys();
};
}
#endif // YTY_COORSYS_H
