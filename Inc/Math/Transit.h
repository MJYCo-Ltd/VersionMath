#ifndef YTY_TRANSIT_H
#define YTY_TRANSIT_H
#include <vector>
#include <Math/Vector.h>

namespace Math{
class MATH_EXPORT CTransit
{
public:
    CTransit();
    ~CTransit();

    /**
     * @brief 验证点是否在区域内
     * @param AllPoint 区域的点
     * @param pt       被判断的点
     */
    static bool ContainsPoint(const std::vector<CVector>& AllPoint,const Math::CVector &pt);
};
}
#endif // YTY_TRANSIT_H
