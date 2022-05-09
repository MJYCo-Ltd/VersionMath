#include <cstdlib>
#include <Math/MemPool.h>
#include "Inc/tlsf.h"

static const unsigned int unDefaultSize=1024*1024;
/// 获取单例
CMemPool *CMemPool::GetInstance()
{
    static CMemPool s_memPool;
    return(&s_memPool);
}

/// 初始化空间大小
bool CMemPool::InitSize(size_t nBytesSize)
{
    if(0 == nBytesSize)
    {
        return(false);
    }


    if(nBytesSize > m_nBytesSize)
    {
        m_nBytesSize = nBytesSize;
    }

    if(nullptr != m_pFirstBuffer)
    {
        if(m_nTotalSize > nBytesSize)
        {
            return(true);
        }
        else
        {
            /// 如果总空间大小小于期望的空间大小
            while(m_nTotalSize < nBytesSize)
            {
                void* pNewArea = malloc(m_nBytesSize);
                m_nTotalSize += m_nBytesSize;
                m_listBuffer.push_back(pNewArea);
                add_new_area(pNewArea,m_nBytesSize,m_pFirstBuffer);
            }
        }
    }
    else
    {
        /// 开辟空间
        m_pFirstBuffer = malloc(m_nBytesSize);
        if( nullptr ==m_pFirstBuffer || size_t(-1) == init_memory_pool(m_nBytesSize,m_pFirstBuffer))
        {
            free(m_pFirstBuffer);
            m_pFirstBuffer = nullptr;
            return(false);
        }
        m_nTotalSize = m_nBytesSize;
    }

    return(true);
}

/// 释放空间
void CMemPool::Release()
{
    if(nullptr != m_pFirstBuffer)
    {
        destroy_memory_pool(m_pFirstBuffer);
        free(m_pFirstBuffer);
        m_pFirstBuffer = nullptr;

        for(auto one : m_listBuffer)
        {
            free(one);
        }

        m_listBuffer.clear();
    }
}

CMemPool::~CMemPool()
{
    Release();
}

template<typename T>
T *CMemPool::Create(size_t nSize)
{
    if(nullptr == m_pFirstBuffer)
    {
        InitSize(unDefaultSize);
        if(nullptr == m_pFirstBuffer)
        {
            return(nullptr);
        }
    }

    /// 尝试开辟空间
    T* pT = (T*)tlsf_malloc(sizeof(T)* nSize);

    while(nullptr == pT)
    {
        /// 开辟新的空间
        void* pNewArea = malloc(m_nBytesSize);
        if(nullptr != pNewArea)
        {
            m_nTotalSize += m_nBytesSize;
            m_listBuffer.push_back(pNewArea);
            add_new_area(pNewArea,m_nBytesSize,m_pFirstBuffer);
            pT = (T*)tlsf_malloc(sizeof(T)* nSize);
        }
        else /// 开辟空间失败则跳出循环
        {
            break;
        }
    }

    return(pT);
}

void CMemPool::Remove(void *pT)
{
    if(pT != nullptr && m_pFirstBuffer != nullptr)
    {
        tlsf_free(pT);
    }
}

template double* CMemPool::Create<double>(size_t nSize);
template char* CMemPool::Create<char>(size_t nSize);
template double** CMemPool::Create<double*>(size_t nSize);
