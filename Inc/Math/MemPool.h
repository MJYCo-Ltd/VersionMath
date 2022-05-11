#ifndef YTY_MEMPOOL_H
#define YTY_MEMPOOL_H

/*****************************************
  作用：统一管理内存，减少内存碎片
 *****************************************/
#include <list>
#include <mutex>
#include "Math_global.h"

class MATH_EXPORT CMemPool
{
public:
    static CMemPool* GetInstance();

    /**
     * @brief 初始化空间大小，以后将以此大小为基本单位开辟空间
     * @param nBytesSize
     */
    bool InitSize(size_t nBytesSize);

    /**
     * @brief 开辟空间
     */
    template<typename T>
    T* Create(size_t nSize);

    /**
     *@brief 移除开辟的空间
     */
    void Remove(void* pT);

    /**
     * @brief 释放空间
     */
    void Release();

protected:
    ~CMemPool();
private:
    size_t m_nBytesSize{};
    size_t m_nTotalSize{};
    void*     m_pFirstBuffer{};
    std::list<void*> m_listBuffer{};
    std::mutex m_mtxMemPool;
};

#endif // YTY_MEMPOOL_H
