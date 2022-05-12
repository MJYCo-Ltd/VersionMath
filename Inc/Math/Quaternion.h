#ifndef YTY_QUATERNION_H
#define YTY_QUATERNION_H
/*****************************************
  作用：四元数
  注意：此类是非线程安全的，在多线程中使用需要对变量
       进行加锁保护，否则运算结果不能保证
 *****************************************/

#include <Math/VecMat.h>
namespace Math{
class MATH_EXPORT CQuaternion
{
public:
    /**
     * @brief 默认构造函数
     */
    CQuaternion();
    ~CQuaternion();

    CQuaternion(const CQuaternion& rOther);             /// 复制构造
    CQuaternion(double dx, double dy, double dz, double ds);

    /**
     * @brief 通过向量与旋转角度构建四元数
     * @param vVec 旋转绕的轴 (vVec的维数必须大于2)
     * @param dphi 旋转角度 单位(弧度)
     */
    CQuaternion(const CVector& vVec, double dphi);

    /**
     * @brief 通过向量构建四元数
     * @param vVec 向量 (vVec的维数必须大于3)
     */
    CQuaternion(const CVector& vVec);

    /**
     * @brief 通过旋转矩阵构造
     * @param vMat 旋转矩阵
     */
    CQuaternion(const CMatrix& vMat);

    /**
     * @brief 将四元数归一化
     */
    void Normalize();

    /**
     * @brief 四元数的模
     * @return
     */
    double Norm() const;

    /**
     * @brief 获取四元数对应的旋转矩阵
     * @return
     */
    CMatrix GetMatrix();

    /**
     * @brief 以向量的形式获取值
     * @return
     */
    CVector GetVector()const;

    /**
     * @brief 获取旋转信息
     * @param rVec   旋转轴
     * @param dAngle 旋转角度 单位(弧度)
     */
    void GetRotate(CVector& rVec, double& dAngle) const;

    /// Set方法
    void SetX(double dX){m_dX = dX;}
    void SetY(double dY){m_dY = dY;}
    void SetZ(double dZ){m_dZ = dZ;}
    void SetS(double dS){m_dS = dS;}

    /// Get方法
    double GetX()const{return(m_dX);}
    double GetY()const{return(m_dY);}
    double GetZ()const{return(m_dZ);}
    double GetS()const{return(m_dS);}


    /// -号操作符 结果与 共轭一样
    CQuaternion& operator -();

    ///
    CQuaternion& operator +=(const CQuaternion& rOther);
    CQuaternion& operator -=(const CQuaternion& rOther);

    /// *= 操作符
    CQuaternion& operator *=(double dScal);
    CQuaternion& operator *=(const CQuaternion& rOther);
    CQuaternion& operator /=(const CQuaternion& rOther);
    CQuaternion& operator /=(double dScal);

    /// 赋值操作符
    CQuaternion& operator =(const CMatrix& rMat);
    CQuaternion& operator =(const CQuaternion& rOther);

    /// 重载等号运算符
    bool operator ==(const CQuaternion& rOther)const;
    bool operator ==(const CMatrix& rMat)const;

    /// 重载bool类型
    operator bool()const;

    bool Rebuild(const CMatrix &vMat);
    bool Rebuild(const CVector &vVec, double dphi);

    bool Translate(const CVector &vIn, CVector &vOut);
private:
    void InitQuat();
    void InitCoefficient();

private:
    double m_dX;
    double m_dY;
    double m_dZ;
    double m_dS;
    double m_dX1,m_dY1,m_dZ1;
    double m_dX2,m_dY2,m_dZ2;
    double m_dX3,m_dY3,m_dZ3;
    bool   m_bNeedInitCoefficient;
};
}
/// 计算rOther的逆四元数
MATH_EXPORT Math::CQuaternion Inv(const Math::CQuaternion& rOther);

/// 计算rOther的共轭四元数
MATH_EXPORT Math::CQuaternion Conj(const Math::CQuaternion& rOther);

/// 计算从 rFrom 到 rTo 的中间位置的四元数
/// dPos 取值范围 (0~1) 超过的取余，为负数则直接乘-1
MATH_EXPORT Math::CQuaternion Slerp(const Math::CQuaternion& rqFrom,const Math::CQuaternion& rqTo,double dPos);

/// 实现运算符
MATH_EXPORT Math::CQuaternion operator + (const Math::CQuaternion& qOne,const Math::CQuaternion& qTwo);
MATH_EXPORT Math::CQuaternion operator - (const Math::CQuaternion& qOne,const Math::CQuaternion& qTwo);
MATH_EXPORT Math::CQuaternion operator * (const Math::CQuaternion& q,double dScale);
MATH_EXPORT Math::CQuaternion operator * (double dScale,const Math::CQuaternion& q);

MATH_EXPORT Math::CVector     operator * (const Math::CQuaternion& q,const Math::CVector& v);
MATH_EXPORT Math::CQuaternion operator * (const Math::CQuaternion& rOne, const Math::CQuaternion& rTwo);
MATH_EXPORT Math::CQuaternion operator / (const Math::CQuaternion& rOne, const Math::CQuaternion& rTwo);

#endif // YTY_QUATERNION_H
