#ifndef YTY_VECTOR_H
#define YTY_VECTOR_H

/*****************************************
  作用：进行向量的运算，包括点乘、叉乘、组合、子集
       、数乘、求模等操作
  注意：此类是非线程安全的，在多线程中使用需要对变量
       进行加锁保护，否则运算结果不能保证
 空向量：维数为0的向量
 *****************************************/
#include <Math/Math_global.h>
namespace Math {

class MATH_EXPORT CVector
{
public:
    const static CVector X_AXIS;
    const static CVector Y_AXIS;
    const static CVector Z_AXIS;
    const static CVector NULL_VECTOR;
    /// 构造函数
    CVector (){}
    CVector (unsigned int nSize);                       /// 构造一个全为零的向量
    CVector (const CVector& rV);           /// 复制构造
    CVector (CVector&& rV);
    CVector (const double* pdArry,unsigned int nDim);  /// 通过数组构造
    CVector (double dX, double dY, double dZ,bool bIsNormal=false); /// 3维向量
    CVector (double dx, double dy, double dz
             ,double dX, double dY, double dZ); /// 6维向量

    /// 析构函数
    ~CVector();

    /// 判断向量是否为空
    inline bool IsEmpty()const {return (0 == m_unDim);}

    /// Set方法
    void Set(double xT,double yT, double zT);
    inline void SetX(double x){m_pdV[0] = x;}
    inline void SetY(double y){m_pdV[1] = y;}
    inline void SetZ(double z){m_pdV[2] = z;}

    /// Get方法
    inline double GetX()const {return m_pdV[0];}
    inline double GetY()const {return m_pdV[1];}
    inline double GetZ()const {return m_pdV[2];}

    /// 大小
    inline unsigned int Size() const { return m_unDim; }
    CVector& Resize(unsigned int nSize);

    /// 获取向量的子集
    CVector slice(unsigned int nFirst, unsigned int nLast) const;

    /// 向量长度
    double Length() const;

    /// 将向量归一化
    void Normalize();
    inline bool IsNormal()const{return(m_bIsNormal);}

    /// 获取向量的开平方向量
    CVector Sqrt() const;

    /// 转换操作符
    inline operator bool()const {return !IsEmpty();}
    /// 赋值操作符
    CVector& operator=(const double value);
    CVector& operator=(const CVector& V);
    CVector& operator=(CVector&& rV);

    /// 获取向量中的值
    /// 注意没有进行越界检查，使用时按照规范使用
    inline double  operator () (int i) const { return m_pdV[i]; }
    inline double& operator () (int i)       { return m_pdV[i]; }

    /// 对向量进行加、减、相等的判断
    void operator += (const CVector& V);
    void operator -= (const CVector& V);
    bool operator == (const CVector& V)const;
    bool operator != (const CVector& V)const;

    /// 对向量进行缩放
    CVector& operator /= (const double value);
    CVector& operator *= (const double value);
private:
    inline void InitAttr();
private:
    /// 成员变量
    unsigned int m_unDim{};  /// 向量维数
    bool   m_bIsNormal{false};
    double *m_pdV{};         /// 向量存储空间
};
}
#endif // YTY_VECTOR_H
