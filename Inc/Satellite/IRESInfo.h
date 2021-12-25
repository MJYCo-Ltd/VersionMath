#ifndef YTY_IRESINFO_H
#define YTY_IRESINFO_H
/*****************************************
  作用：本文件负责读取finals2000A.data文件，
       并解析文件内容
  注意：此类是非线程安全的，在多线程中使用需要对变量
       进行加锁保护，否则运算结果不能保证
  结构：采用单例设计模式进行设计
数据获取：finals2000A.data 下载地址 www.iers.org
        如地址不可达则直接搜索 IERS 进行数据下载。
 *****************************************/
#include "SAT_global.h"
#include <string>
#include <vector>

using namespace std;

namespace Aerospace{
struct BulletinB
{
    double dPX{};           /// 极移x    角度单位 秒
    double dPY{};           /// 极移y    角度单位 秒
    double dUT1_UTC{};      /// UT1-UTC 时间单位 秒
};

class ALGORITHM_EXPORT CIRESInfo
{
public:
    /// 获取实例
    static CIRESInfo  *GetInstance();

    /// 通过文件初始化
    bool Init(const string& strFileName);
    /// 判断是否初始化
    bool IsInit(){return m_bInit;}

    /// 根据简约儒略日获取布告B中的数据
    bool GetBulletinB(double dMJD, BulletinB &bulletinB) const;
    /// 获取
    double      GetUT2_UT1(double dMJD) const;
private:
    /// 清除数据，在再次初始化时使用
    void ClearData();
    CIRESInfo();
    ~CIRESInfo();
private:
    vector<double>     m_vMJD;
    vector<double>     m_vPx;
    vector<double>     m_vPy;
    vector<double>     m_vUt1_Utc;
    unsigned int       m_unCount{};
    bool               m_bInit{};                 //是否初始化
};
}
#endif // YTY_IRESINFO_H
