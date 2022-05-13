#ifndef YTY_SGP4_H
#define YTY_SGP4_H

/*****************************************
  作用：进行SGP4/SDP4的轨道运算
  注意：所有的操作不成功会返回一个空的向量，
       使用时需要注意，如果调用函数时不能保证输入
       数值的正确性，需要检测对象是否可为空，或者
       直接判断对象是否合理。例如：
       CSGP4 tleSGP4;
       …
       if(tleSGP4)
       {
         代码;
       }
       …
       TLE计算结果为TEME of DATA 坐标系下的坐标，
       如果转换成 J2000或者 ECF请参考 CCoorSys

       没有考虑多线程，如果在多线程下使用，需要自己进行
       数据同步控制
 *****************************************/

#include <string>
#include <Satellite/SAT_global.h>
#include <Math/Vector.h>

struct elsetrec;

const int TLELENGTH(73);
const unsigned int FIRSTVALUELENGTH(61);
const unsigned int SECONDVALUELENGTH(63);

namespace Satellite{

class ALGORITHM_EXPORT CSGP4
{
public:
    CSGP4(const std::string& strLine1, const std::string& strLine2);
    CSGP4(const char strLine1[73],const char strLine2[73]);
    ~CSGP4();

    /**
     * @brief 设置两行数据
     * @param strLine1  第一行数据
     * @param strLine2  第二行数据
     * @return 如果数据有效则返回 true 无效则返回 false
     */
    bool SetTLE(const std::string& strLine1, const std::string& strLine2);

    /**
     * @brief 计算TEME坐标系下的位置速度
     * @param dMJD  约简儒略日
     * @return 错误则向量为空 正确则返回 位置(x,y,z)速度(vx,vy,vz)
     *                                  [m,m/s][TEME]
     */
    Math::CVector CalPV(const double& dMJD);

    /**
     * @brief 重载bool类型，用于判断是否有效
     */
    inline operator bool()const{return(m_bValid);}

    /**
     * @brief 获取TLE中包含的经典六根数
     * @return 开普勒六根数  (半长轴,离心率,轨道倾角,升交点赤经,近地点幅角,平近点角)
     *                     [m, ,rad,rad,rad,rad]
     * @attention 只有数据有效才能调用时值才有效
     */
    const Math::CVector& ClassicalElements(){return (m_Oribit);}

    /**
     * @brief 通过位置速度计算 经典六根数
     * @param vRV(x,y,z,vx,vy,vz) [m,m/s]
     * @return 开普勒六根数  (半长轴,离心率,轨道倾角,升交点赤经,近地点幅角,平近点角)
     *                         [m, ,rad,rad,rad,rad]
     * @attention 如果计算失败会返回一个空的向量
     */
    static Math::CVector ClassicalElements(const Math::CVector& vRV);

    /**
     * @brief 获取TLE的历元时间
     * @return 约简儒略日
     * @attention 只有数据有效调用时值才有效
     */
    const double& GetTLEEpoch(){return(m_dEpochMJD);}

    /**
     * @brief 获取tle数据
     * @return
     */
    const elsetrec* GetTleInfo(){return(m_pTle);}

private:
    char   m_strTLE[2][TLELENGTH]; /// 两行根数数据
    bool   m_bValid;               /// 数据是否有效
    Math::CVector m_Oribit;        /// 经典六根数信息
    double m_dEpochMJD;            /// TLE对应的历元时间
    elsetrec* m_pTle;              /// 读取到的TLE数据
};
}
#endif // YTY_SGP4_H
