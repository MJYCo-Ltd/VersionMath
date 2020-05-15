#ifndef YTY_YPRANGLE_H
#define YTY_YPRANGLE_H
/*****************************************
  作用：俯仰翻滚偏航 旋转转换
  注意：此类是非线程安全的，在多线程中使用需要对变量
       进行加锁保护，否则运算结果不能保证
 *****************************************/
#include "Math_global.h"

//////////////////////////////////////////////////
// x 对应 Roll
// y 对应 Pitch
// z 对应 yaw
///////////////////////////////////////////////////

namespace Math {

class CMatrix;

/**
 * @brief 俯仰翻滚偏航的旋转类型
 */
enum YPRROTATE
{
    XYZ,XZY,YXZ,YZX,ZXY,ZYX
    ,RPY,RYP,PRY,PYR,YRP,YPR
};

/**
 * @brief 俯仰翻滚偏航的旋转值
 */
struct YPR_Rotate
{
    union
    {
        double dX;
        double dRoll;
    };

    union
    {
        double dY;
        double dPitch;
    };

    union
    {
        double dZ;
        double dYaw;
    };
};

class MATH_EXPORT CYPRAngle
{
public:
    /// 构造函数
    CYPRAngle(const CMatrix& rMatrix,YPRROTATE type=RPY);
    CYPRAngle(double dRoll, double dPitch, double dYaw, YPRROTATE type=RPY);
    CYPRAngle(const YPR_Rotate& stRotate, YPRROTATE type=PRY);
    ~CYPRAngle();

    /**
     * @brief 通过初始旋转类型与旋转角度，计算所求的类型的旋转角度
     * @param OrigType  原始旋转类型
     * @param PurType   目标旋转类型
     * @param OrigRotate原始旋转角度 [rad]
     * @param PurRotate 目标旋转角度 [rad]
     * @return 计算成功返回 true 失败返回 false
     */
    static bool CalTransform(YPRROTATE OrigType, YPRROTATE PurType, const YPR_Rotate& OrigRotate,
                             YPR_Rotate& PurRotate);

    /**
     * @brief 根据矩阵计算，所求的类型的旋转角度
     * @param rotateMatrix 旋转矩阵
     * @param type         所求的旋转类型
     * @param stRotate     所求的旋转角度 [rad]
     */
    static void CalTransform(const CMatrix& rotateMatrix, YPRROTATE type, YPR_Rotate& stRotate);

    /**
     * @brief 通过翻滚、俯仰、偏航角度值，以及旋转类型，计算旋转矩阵
     * @param dRoll   翻滚角 [rad]
     * @param dPitch  俯仰角 [rad]
     * @param dYaw    偏航角 [rad]
     * @param type    旋转类型
     * @return 旋转矩阵
     */
    static CMatrix CreateMatrix(double dRoll,double dPitch, double dYaw, YPRROTATE type);


    /**
     * @brief 获取指定类型的角度序列值
     * @param type 要获取的指定类型
     * @return 各角度的数值 [rad]
     */
    YPR_Rotate GetRotate(YPRROTATE type) const;

    /**
     * @brief 获取存储的旋转信息
     * @return 各角度的数值 [rad]
     */
    const YPR_Rotate& GetRotate()const{return(m_stRotate);}

    /**
     * @brief 获取旋转矩阵
     * @return
     */
    const CMatrix& GetMatrix()const{return(m_matRote);}

private:
    YPRROTATE m_emOriginalType;   // 原始旋转类型
    CMatrix m_matRote;        // 旋转矩阵
    YPR_Rotate     m_stRotate;    // 旋转角度信息
};
}
#endif // YTY_YPRANGLE_H
