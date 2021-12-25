#ifndef YTY_COORSYS_H
#define YTY_COORSYS_H

/*****************************************
  作用：完成各种坐标系的相互转换
       本类采用的ECI为GCRF；ECF为ITRF
  注意：此类是非线程安全的，在多线程中使用需要对变量
       进行加锁保护，否则运算结果不能保证
名词解释：
                           地固系[ECF]
        GTRS(Geocentric Terrestrial Reference System)    地心地球参考系
        ITRS(International Terrestrial Reference System) 国际地球参考系
        ITRF(International Terrestrial Reference Frame)  国际地球参考架
        VVLH()                                           轨道参考系

                           地惯系[ECI]
        GCRF(Geocentric Celestial Reference Frame)       地心天球参考架
        GCRS(Geocentric Celestial Reference System)      地心天球参考系
        ICRS(International Celestial Reference System)   国际天球参考系
        ICRF(International Celestial Reference Frame)    国际天球参考架
        CIRS(Celestial Intermediate Reference System)    天球中间参考系
        J2000(Mean Equator and Mean Equinox of 2000)     J2000的平赤道平春分点坐标系
        MOD(Mean Equator and Mean Equinox of Date)       平赤道平春分点坐标系
        TOD(True Equator and True Equinox of Date)       真赤道真春分点坐标
        TEME(True Equator and Mean Equinox)              真赤道平春分点坐标系

                         其他坐标系
        ENU(East-North-Zenith coordinate System)         东北天坐标系(站心坐标系)
 *****************************************/
#include <Satellite/SAT_global.h>

namespace Math{
class CVector;
class CMatrix;
}

namespace Aerospace{
using namespace Math;

class ALGORITHM_EXPORT CCoorSys
{
public:

    /////////////////// TEME 与 J2000 坐标转换///////////////////////////
    /// @brief  TEME到ICRS的坐标转换
    /// @param  dJD     约简儒略日
    /// @param  VTeme   TEMEOfDate下的坐标
    /// @param  vJ2000  J2000下的坐标
    /// @return 成功返回true 失败返回false
    /// @attention 支持批量处理
    static bool TEME2J2000(const double& dMJD, const CVector &vTeme, CVector &vJ2000);
    static bool J20002TEME(const double& dMJD, const CVector &vJ2000, CVector &vTeme);
    /// @return TEME 转 J2000 的旋转矩阵
    static CMatrix TEME2J2000(const double& dMJD);
    static CMatrix J20002TEME(const double& dMJD);

    /////////////////// J2000 与 其他惯性系 坐标转换///////////////////////////
    static CMatrix J20002GCRF();
    static CMatrix GCRF2J2000();

    static CMatrix J20002MOD(const double& dMJD);
    static CMatrix MOD2J2000(const double& dMJD);

    static CMatrix J20002TOD(const double& dMJD);
    static CMatrix TOD2J2000(const double& dMJD);

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
    static bool VVLH2J2000(const CVector &vKep, const CVector& vVVLH, CVector& vECI);
    static bool J20002VVLH(const CVector &vKep, const CVector& vECI, CVector& vVVLH);
    /// ECI到轨道坐标系的旋转矩阵
    static CMatrix VVLH2J2000(const CVector& vKep);
    static CMatrix J20002VVLH(const CVector& vKep);
};
}
#endif // YTY_COORSYS_H
