#ifndef YTY_EULERANGLE_H
#define YTY_EULERANGLE_H
/*****************************************
  作用：欧拉角旋转
  注意：此类是非线程安全的，在多线程中使用需要对变量
       进行加锁保护，否则运算结果不能保证
 *****************************************/
#include "Math_global.h"

namespace Math {
/**
 * @brief 欧拉角的旋转类型
 */
enum EULERROTATE
{
    E121,E123,E131,E132,E212,E213
    ,E231,E232,E312,E313,E321,E323
};

/**
 * @brief 欧拉角度值
 */
struct Euler_Rotate
{
    double dA;
    double dB;
    double dC;
};

class CMatrix;

class MATH_EXPORT CEulerAngle
{
public:
    CEulerAngle(const CMatrix& rMatrix,EULERROTATE type=E123);
    CEulerAngle(double dA, double dB, double dC, EULERROTATE type=E123);
    CEulerAngle(const Euler_Rotate& stRotate, EULERROTATE type=E123);
    ~CEulerAngle();

    /**
     * @brief 通过角A、角B、角C角度值，以及旋转类型，计算旋转矩阵
     * @param dA      角A [rad]
     * @param dB      角B [rad]
     * @param dC      角C [rad]
     * @param type    旋转类型
     * @return 旋转矩阵
     */
    static CMatrix CreateMatrix(double dA,double dB, double dC, EULERROTATE type);

    /**
     * @brief 获取旋转矩阵
     * @return 旋转矩阵的值
     */
    const CMatrix& GetMatrix()const{return(m_matRote);}

    /**
     * @brief 获取旋转角度信息
     * @return
     */
    const Euler_Rotate& GetRotate()const{return(m_stRotate);}

    /**
     * @brief 获取旋转类型
     * @return
     */
    EULERROTATE GetRotateType()const{return(m_emOriginalType);}

private:
    EULERROTATE  m_emOriginalType;   // 原始旋转类型
    CMatrix  m_matRote;          // 旋转矩阵
    Euler_Rotate m_stRotate;         // 旋转角度信息
};
}
#endif // YTY_EULERANGLE_H
