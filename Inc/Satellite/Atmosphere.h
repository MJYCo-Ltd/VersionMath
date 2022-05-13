#ifndef YTY_ATMOSPHERE_H
#define YTY_ATMOSPHERE_H
/*****************************************
  作用：大气密度模型
  注意：此类是非线程安全的，在多线程中使用需要对变量
       进行加锁保护，否则运算结果不能保证
  描述：HarrisPriester 模型范围 100km~1000km
       NRLMSISE2000   模型范围  >-10km
       STANDARD1976   模型范围 0km~1000km
 *****************************************/

#include <string>
#include <Satellite/SAT_global.h>
#include <Math/Vector.h>


enum AtmosModel
{
    HarrisPriester,NRLMSISE2000,STANDARD1976
};

namespace Physical{
class ALGORITHM_EXPORT CAtmosphere
{
public:
    CAtmosphere(AtmosModel type=HarrisPriester);
    ~CAtmosphere();

    //    /**
    //     * @brief 用文件初始化大气模型
    //     * @param csFileName
    //     * @return
    //     */
    //    bool Init(const char* csFileName);

    /**
     * @brief 判断大气模型是否初始化
     * @return
     */
    bool IsInit(){return(m_bInit);}

    /**
     * @brief 获取大气模型名称
     * @return
     */
    const std::string& GetModelName(){return (m_strModelName);}

    /**
     * @brief 获取大气模型的类型
     * @return
     */
    AtmosModel GetModelType(){return(m_typeAtmos);}

    /**
     * @brief 获得大气密度值
     * @param dMJD                  [MJD]
     * @param rSat 卫星ECF位置 (x,y,z)[m]
     * @param rSun 太阳ECI位置 (x,y,z)[m]
     * @return 大气密度 [kg/m^3]
     */
    double GetDensity(const double& dMJD,const Math::CVector& rSat,const Math::CVector& rSun);

    /**
     * @brief 获得大气密度值
     * @param dMJD                 [MJD]
     * @param rSat 卫星ECF位置 (x,y,z)[m]
     * @param rSun 太阳ECI位置 (x,y,z)[m]
     * @param typeAtmos 大气类型
     * @return大气密度 [kg/m^3]
     */
    static double GetDensity(const double &dMJD, const Math::CVector &rSat, const Math::CVector &rSun,
                             AtmosModel typeAtmos);

protected:
    void SetModelName();

private:
    AtmosModel m_typeAtmos;   /// 大气模型类型
    std::string m_strModelName;/// 大气名称
    bool       m_bInit;       /// 是否初始化成功
};
}
#endif // YTY_ATMOSPHERE_H
