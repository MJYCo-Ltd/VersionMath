﻿#ifndef YTY_STKGRAVEMODEL_H
#define YTY_STKGRAVEMODEL_H
/*****************************************
  作用：本文件负责读取STK中的重力场系数文件，
       并解析文件内容
  注意：此类是非线程安全的，在多线程中使用需要对变量
       进行加锁保护，否则运算结果不能保证
数据获取：STK安装目录下 STKData\CentralBodies\Earth
        所有的 grv 文件即为重力场系数模型文件
 *****************************************/
#include <string>
using namespace std;
#include "SAT_global.h"

namespace Math{
class CMatrix;
}

namespace Phycical{
using namespace Math;
class ALGORITHM_EXPORT CSTKGraveModel
{
public:
    CSTKGraveModel();
    CSTKGraveModel(const string& sFileName);
    virtual ~CSTKGraveModel(){}

    /**
     * 设置文件名字
     */
    bool SetFile(const string& sFileName);

    /**
     * 获取模型名称
     */
    const string& GetModelName(){return m_sModelName;}

    /**
     *
     * 获取阶数
     */
    const int& GetDegree(){return m_nDegree;}

    /**
     * 获取次数
     */
    const int& GetOrder(){return m_nOrder;}

    /**
     * 获取地球的引力常数
     */
    const double& GetEarthGM(){return m_dGM;}

    /**
     * 获取地球半径
     */
    const double& GetEarthRref(){return m_dR;}

    /**
     * 获取归一化重力系数
     */
    const CMatrix& GetCoefficients(){return m_mtInfo;}

    /**
     * 获取非归一化重力系数
     */
    const CMatrix& GetUnNormalizeCoefficients(){return m_unNorInfo;}

    /**
     * 查看是否准备好
     */
    bool GetReady(){return m_bReady;}
private:
    /**
     * 读取文件
     */
    bool ReadFile();

    /**
     * 初始化缓存
     */
    void InitBuffer(const streamsize &nSize);

    /**
     * 处理缓存数据
     */
    void DealBuffer();
    bool ReadCS();
    bool ReadStr(string& rsName,const char* cstrSign);
    bool ReadDouble(double& dDouble, const char* cstrSign);
    bool ReadInt(int &nInt, const char *cstrSign);
    double UnNormaliz(int nN, int nM);
private:
    string m_sFileName;      /// 文件名称
    string m_sModelName;     /// 地球引力模型名称

    int    m_nDegree;        /// 地球引力模型的阶数
    int    m_nOrder;         /// 地球引力模型的次数

    double m_dGM;            /// 地球万有引力常数
    double m_dR;             /// 地球的半径

    CMatrix m_mtInfo;    /// 重力场系数
    CMatrix m_unNorInfo; /// 非归一化重力场系数

    bool   m_bReady;         /// 数据是否处理完毕
    bool   m_bNormal;        /// 是否是归一化的系数

    char*  m_pBuffer;        /// 将文件读取到内存中的缓存
};
}
#endif // YTY_STKGRAVEMODEL_H
