#ifndef YTY_MATRIX_H
#define YTY_MATRIX_H
/*****************************************
  作用：进行矩阵的数据的保存，以及简单的操作包括：
       矩阵的加减，矩阵的缩放
  注意：此类是非线程安全的，在多线程中使用需要对变量
       进行加锁保护，否则运算结果不能保证
 空矩阵：行数列数都为零的矩阵
 *****************************************/

#include "Math_global.h"
#include "Vector.h"

/// 世界坐标系 不随着本体旋转而旋转的 坐标系
/// 局部坐标系  跟物体固联的坐标系，一直保持机头为x轴，机身右侧为y轴，机身背部为z轴
/// 向量 矩阵在左边 (M*v)  得到的向量 表示 世界坐标系中这个点的位置在局部坐标中的描述
/// 向量 矩阵在右边 (v*M)  得到的向量 表示 局部坐标系中的这个点的位置在世界坐标中的描述

namespace Math {

class MATH_EXPORT CMatrix
{

public:
    const static CMatrix IDENTITY_MAT;
    const static CMatrix NULL_MATRIX;
    /// 构造函数
    CMatrix ();                                      /// 默认构造函数
    CMatrix (unsigned int dim1,unsigned int dim2);   /// Nullmatrix
    CMatrix (const CMatrix& rM);                     /// 复制构造函数
    CMatrix (const double *p,unsigned int dim1,unsigned int dim2);     /// 通过一维数组构造
    CMatrix (const double p[3][3],unsigned int dim1,unsigned int dim2);/// 通过二维数组构造

    /// 析构函数
    ~CMatrix();

    /// 判断矩阵是否为空
    inline bool IsEmpty()const{return (0 == m_unRow || 0 == m_unCol);}
    inline bool IsSquare()const{return(m_unRow == m_unCol && m_unRow > 0);}

    /// 判断矩阵是否是单位矩阵
    bool IsId() const;

    bool Translate(const CVector &vIn, CVector &vOut)const;

    /// 赋值
    /// 将所有的元素赋值成指定的值
    CMatrix& operator=(const double value);

    /// 如果矩阵为空，赋值后两个矩阵相同
    /// 如果矩阵非空，只有当两个矩阵行数与列数相同才能进行赋值操作：
    /// 1、行数与列数不同不会造成左操作数的改变
    /// 2、行数与类书相同，赋值后两个矩阵相同
    CMatrix& operator=(const CMatrix& rM);

    /// 获取大小
    inline unsigned int Row() const { return m_unRow; }
    inline unsigned int Col() const { return m_unCol; }
    /// 重置矩阵大小
    CMatrix& Resize(unsigned int nRow,unsigned int nCol);

    /// 返回指定列数据
    CVector GetCol(unsigned int j) const;

    /// 返回指定行数据
    CVector GetRow(unsigned int i) const;

    /// 获取对角线上的数据
    /// @attention 如果无法获取则返回一个空的向量
    /// 可能失败的原因：1、行数列数不相等；2、矩阵为空
    CVector Diag() const;

    /// 矩阵的迹
    double Trace() const;
    double Trace(unsigned int low,unsigned int upp) const;

    /// 获取矩阵的一块儿
    /// 如果指定的区域超出矩阵的范围则返回一个空的矩阵
    CMatrix slice(unsigned int first_row, unsigned int last_row, unsigned int first_col, unsigned int last_col);

    /** @brief     将向量设置到指定列
      * @param j   指定矩阵的列
      * @param Col 准备填充的向量
      * @return    如果成功返回true
      */
    bool SetCol(unsigned int j, const CVector& Col);

    /// 设置指定行成指定数据
    bool SetRow(unsigned int i, const CVector& Row);

    /// 转换操作符
    inline operator bool() const {return !IsEmpty();}
    /// 访问数据
    inline double  operator () (unsigned int i, unsigned int j) const { return m_ppdM[i][j]; }
    inline double& operator () (unsigned int i, unsigned int j)       { return m_ppdM[i][j]; }
    /// 矩阵的加减
    /// 如果传入的矩阵的行列不同，则保持原值不变
    bool operator += (const CMatrix& rM);
    bool operator -= (const CMatrix& rM);
    bool operator == (const CMatrix& rM)const;
    bool operator != (const CMatrix& rM)const;
    /// 对矩阵进行缩放
    CMatrix& operator /= (const double value);
    CMatrix& operator *= (const double value);

private:
    /// 将矩阵置空
    void Empty();
    /// 初始化属性
    inline void InitAttr();
private:

    /// 属性
    unsigned int      m_unRow;  // First dimension (number of rows)
    unsigned int      m_unCol;  // Second dimension (number of columns)
    double **m_ppdM;           // Matrix M(n,m)

};
}
#endif // YTY_MATRIX_H
