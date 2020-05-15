#ifndef YTY_TRANSIT_H
#define YTY_TRANSIT_H
#include <vector>
#include "Vector.h"
using namespace std;
#include "SAT_global.h"

using namespace Math;
namespace Satellite{
class ALGORITHM_EXPORT CTransit
{
public:
    CTransit();
    ~CTransit();

    /**
     * @brief 验证点是否在区域内
     * @param AllPoint 区域的点
     * @param pt       被判断的点
     */
    static bool ContainsPoint(const vector<CVector>& AllPoint,const CVector &pt);
};
}
#endif // YTY_TRANSIT_H
