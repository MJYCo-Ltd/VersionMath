#include <cstring>
#include <cmath>
using namespace std;
#include <Math/Vector.h>
#include <Math/Intpol.h>
#include <Math/MemPool.h>

using namespace Numerical;
using namespace Math;

/// 内维尔插值算法
double Intpol::ItNeville(int nNum, const double *pdX, const double *pdY, double dX0)
{
    int     i,k;                    // Loop counters
    double* p = CMemPool::GetInstance()->Create<double>(nNum);  // Interpolation tableau
    double  dY0;

    /// 将值复制到新开辟的空间
    memcpy(p,pdY,sizeof(*p)*nNum);

    /// Neville 插值算法
    for (k=1;k<nNum;++k)
    {
        for (i=nNum-1;i>=k;--i)
        {
            p[i] += (dX0-pdX[i])*(p[i]-p[i-1])/(pdX[i]-pdX[i-k]);
        }
    }

    dY0 = p[nNum-1];

    /// 删除开辟的空间
    CMemPool::GetInstance()->Remove(p);

    return(dY0);
}

/// 拉格朗日插值算法
double Intpol::ItLagrange(int nNum, const double *pdX, const double *pdY, double dX0)
{
    int i, j;
    double dY0(0);

    /// Lagrange 插值算法
    for (i = 0; i < nNum; ++i)
    {
        if (fabs(dX0 - pdX[i]) < 1e-8)
        {
            return pdY[i];
        }

        double X = 1;
        for (j = 0; j < nNum; ++j)
        {
            if (i != j)
            {
                X = X * (dX0 - pdX[j]) / (pdX[i] - pdX[j]);
            }
        }
        dY0 += X * pdY[i];
    }

    return (dY0);
}

/// 牛顿插值
double Intpol::ItNewton(int nNum, const double *pdX, const double *pdY, double dX0)
{
    int i, j;
    double* p = CMemPool::GetInstance()->Create<double>(nNum);
    double  dY0;

    /// 将值复制到新开辟的空间
    memcpy(p,pdY,sizeof(*p)*nNum);

    /// Newton 插值算法
    for (j = 0; j < nNum - 1; ++j)
    {
        for (i = nNum - 1; i > j; --i)
        {
            p[i] = (p[i] - p[i - 1]) / (pdX[i] - pdX[i - j - 1]);
        }
    }

    dY0 = p[nNum - 1];

    for (i = nNum - 2; i >= 0; --i)
    {
        dY0 = dY0 * (dX0 - pdX[i]) + p[i];
    }

    /// 删除开辟的空间
    CMemPool::GetInstance()->Remove(p);

    return (dY0);
}

/// 三次样条插值
double Intpol::ItCubicSpline(int nNum, const double *pdX, const double *pdY, double dX0)
{

    if(nNum <2)
    {
        return(0.0);
    }

    /// i 作为循环索引
    /// nIndex作为保存所求所在的区域
    int i(0),nIndex(0);

    /// 防止越界
    if(pdX[0] >= dX0)
    {
        return(pdY[0]);
    }

    /// 防止越界
    if(pdX[nNum - 1] <= dX0)
    {
        return(pdY[nNum-1]);
    }

    /// 查找时间落在哪个区域中
    for(i = 1; i<nNum; ++i)
    {
        if(dX0 == pdX[i])
        {
            return(pdY[i]);
        }

        if(dX0 > pdX[i-1] && dX0<pdX[i])
        {
            nIndex = i;
            break;
        }
    }

    /// 如果值没有在边界点上，进行三次样条插值
    double* dis = CMemPool::GetInstance()->Create<double>(nNum); /// 保存当前点到第一个点的距离
    double* H   = CMemPool::GetInstance()->Create<double>(nNum); /// 保存相邻两点间距离
    double dTempX0 = dX0-pdX[0];

    /// 初始化第一个点
    dis[0] = 0.0 ;
    H[0]   = 0.0;

    //========== 遍历关键点计算相邻两点间的时间  ========================
    for(i = 0 ; i < nNum-1 ; ++i)
    {
        double dTemp = pdX[i+1]-pdX[i];

        dis[i+1] = dis[i]+dTemp;
        H[i+1]=dTemp;
    }
    //================ END =========================================

    double* A = CMemPool::GetInstance()->Create<double>(nNum); /// Miu
    double* B = CMemPool::GetInstance()->Create<double>(nNum); /// 对角线
    double* C = CMemPool::GetInstance()->Create<double>(nNum); /// Lanmida
    double* D = CMemPool::GetInstance()->Create<double>(nNum); //
    double* COEX = CMemPool::GetInstance()->Create<double>(nNum); // 保存X的
    double* BTA = CMemPool::GetInstance()->Create<double>(nNum); // 二阶导
    double* CTA = CMemPool::GetInstance()->Create<double>(nNum); // 一阶导

    //================ 已知相邻两个关键点的时间差求比率 =================================================
    for( i = 1 ; i < nNum-1 ; ++i)
    {
        A[i] = H[i] / (H[i]+H[i+1]); // 计算三点当前点的比率
        B[i] = 2 ;
        C[i] = 1 - A[i]; // 剩下的比率
    }
    A[0] = 1 ;
    B[0] = 2 ;
    C[0] = 0 ;
    A[nNum-1] = 0 ;
    B[nNum-1] = 2 ;
    C[nNum-1] = 1 ;
    //================ END  =================================================

    /// get  coef of function " x = x(s) "
    // 计算两个关键点之间关于X的曲线方程

    //===================  X的D[] ===============================================
    for( i = 1 ; i < nNum-1 ; ++i)
    {
        D[i] = 3*( C[i]*(pdY[i]-pdY[i-1])/H[i] + A[i]*(pdY[i+1]-pdY[i])/H[i+1] );
    }

    D[0]   =  3.0*(pdY[1]   - pdY[0])  /H[1];
    D[nNum-1] =  3.0*(pdY[nNum-1] - pdY[nNum-2])/H[nNum-1];
    //=================== X的D[] END ===============================================

    //=================== X的BTA[]  ===============================================
    BTA[0] = C[0] / B[0] ;
    for( i = 1 ; i < nNum-1 ; ++i)
    {
        BTA[i] = C[i] / (B[i] - A[i]*BTA[i-1]);
    }
    //=================== X的BTA[] END ===============================================

    //=================== X的CTA[]   ===============================================
    CTA[0] = D[0] / B[0] ;
    for( i = 1 ; i < nNum ; ++i)
    {
        CTA[i] = (D[i] - A[i]*CTA[i-1]) / (B[i] - A[i]*BTA[i-1]) ;
    }
    //=================== X的CTA[] END ===============================================

    //===================  X的COEX[]  ===============================================
    COEX[nNum-1] = CTA[nNum-1] ;
    for( i = nNum-2 ; i >= 0 ; --i)
    {
        COEX[i] = CTA[i] - BTA[i]*COEX[i+1];
    }
    //=================== X的COEX[] END ===============================================

    //=================== 根据方程系数计算各个段的插入点坐标 ==============================
    double dResult =  COEX[nIndex-1]*(dis[nIndex]-dTempX0)*(dis[nIndex]-dTempX0)*(dTempX0-dis[nIndex-1])/(H[nIndex]*H[nIndex])
                     -COEX[nIndex]*(dTempX0-dis[nIndex-1])*(dTempX0-dis[nIndex-1])*(dis[nIndex]-dTempX0)/(H[nIndex]*H[nIndex])
                    +pdY[nIndex-1]*(dis[nIndex]-dTempX0)*(dis[nIndex]-dTempX0)*(2*(dTempX0-dis[nIndex-1])+H[nIndex])/(H[nIndex]*H[nIndex]*H[nIndex])
                    +pdY[nIndex]*(dTempX0-dis[nIndex-1])*(dTempX0-dis[nIndex-1])*(2*(dis[nIndex]-dTempX0)+H[nIndex])/(H[nIndex]*H[nIndex]*H[nIndex]);


    /// 释放内存
    CMemPool::GetInstance()->Remove(dis);
    CMemPool::GetInstance()->Remove(H);
    CMemPool::GetInstance()->Remove(A);
    CMemPool::GetInstance()->Remove(B);
    CMemPool::GetInstance()->Remove(C);
    CMemPool::GetInstance()->Remove(D);
    CMemPool::GetInstance()->Remove(COEX);
    CMemPool::GetInstance()->Remove(BTA);
    CMemPool::GetInstance()->Remove(CTA);

    return(dResult);
}

CVector Intpol::ItCubicSpline(const vector<double> &vX, const vector<CVector> &vY, double dX)
{
    /// 定义一个临时变量，用于返回
    CVector tmp;

    int nXSize = vX.size();
    int nYSize = vY.size();
    int nDim   = vY[0].Size();
    /// 判断两个数组是否相等，不相等则直接退出
    if(nXSize != nYSize || nXSize < 2 || nDim < 1)
    {
        return (tmp);
    }

    /// i 作为循环索引
    /// nIndex作为保存所求所在的区域
    int i(0),nFind(0);

    /// 防止越界
    if(vX[0] >= dX)
    {
        return(vY[0]);
    }

    /// 防止越界
    if(vX[nXSize-1] <= dX)
    {
        return(vY[nXSize-1]);
    }

    /// 查找时间落在哪个区域中
    for(i = 1; i<nXSize; ++i)
    {
        if(dX == vX[i])
        {
            return(vY[i]);
        }

        if(dX > vX[i-1] && dX < vX[i])
        {
            nFind = i;
            break;
        }
    }

    double* dis = CMemPool::GetInstance()->Create<double>(nXSize); // 保存当前点到第一个点的距离
    double* H   = CMemPool::GetInstance()->Create<double>(nXSize); // 保存相邻两点间距离
    double dTempX = dX-vX[0];

    // 初始化第一个点
    dis[0] = 0.0 ;
    H[0]   = 0.0;

    //========== 遍历关键点计算相邻两点间的时间  ========================
    for(i = 0 ; i < nXSize-1 ; ++i)
    {
        double dTemp = vX[i+1]-vX[i];

        dis[i+1] = dis[i]+dTemp;
        H[i+1]=dTemp;
    }
    //================ END =========================================

    double* A = CMemPool::GetInstance()->Create<double>(nXSize); /// Miu
    double* B = CMemPool::GetInstance()->Create<double>(nXSize); /// 对角线
    double* C = CMemPool::GetInstance()->Create<double>(nXSize); /// Lanmida
    double* D = CMemPool::GetInstance()->Create<double>(nXSize); //
    double** COEX = CMemPool::GetInstance()->Create<double*>(nDim); // 保存X的
    for(i=0; i<nDim; ++i)
    {
        COEX[i] = CMemPool::GetInstance()->Create<double>(nXSize);
    }
    double* BTA = CMemPool::GetInstance()->Create<double>(nXSize); // 二阶导
    double* CTA = CMemPool::GetInstance()->Create<double>(nXSize); // 一阶导

    //================ 已知相邻两个关键点的时间差求比率 =================================================
    for( i = 1 ; i < nXSize-1 ; ++i)
    {
        A[i] = H[i] / (H[i]+H[i+1]); // 计算三点当前点的比率
        B[i] = 2 ;
        C[i] = 1 - A[i]; // 剩下的比率
    }
    A[0] = 1 ;
    B[0] = 2 ;
    C[0] = 0 ;
    A[nXSize-1] = 0 ;
    B[nXSize-1] = 2 ;
    C[nXSize-1] = 1 ;
    //================ END  =================================================

    //=================== BTA[]  ===============================================
    BTA[0] = C[0] / B[0] ;
    for( i = 1 ; i < nXSize-1 ; ++i)
    {
        BTA[i] = C[i] / (B[i] - A[i]*BTA[i-1]);
    }
    //=================== BTA[] END ===============================================

    /// 遍历各个维度的数据，进行计算
    for(int nIndex = 0; nIndex<nDim; ++nIndex)
    {
        //===================  D[] ===============================================
        for( i = 1 ; i < nXSize-1 ; ++i)
        {
            D[i] = 3*( C[i]*(vY[i](nIndex)-vY[i-1](nIndex))/H[i] + A[i]*(vY[i+1](nIndex)-vY[i](nIndex))/H[i+1] );
        }

        D[0]   =  3.0*(vY[1](nIndex)   - vY[0](nIndex))  /H[1];
        D[nXSize-1] =  3.0*(vY[nXSize-1](nIndex) - vY[nXSize-2](nIndex))/H[nXSize-1];
        //=================== D[] END ===============================================

        //=================== CTA[]   ===============================================
        CTA[0] = D[0] / B[0] ;
        for( i = 1 ; i < nXSize ; ++i)
        {
            CTA[i] = (D[i] - A[i]*CTA[i-1]) / (B[i] - A[i]*BTA[i-1]) ;
        }
        //=================== CTA[] END ===============================================

        //===================  COEX[]  ===============================================
        COEX[nIndex][nXSize-1] = CTA[nXSize-1] ;
        for( i = nXSize-2 ; i >= 0 ; --i)
        {
            COEX[nIndex][i] = CTA[i] - BTA[i]*COEX[nIndex][i+1];
        }
        //=================== COEX[] END ===============================================
    }

    tmp.Resize(nDim);
    for(i=0; i<nDim; ++i)
    {
        tmp(i) =  COEX[i][nFind-1]*(dis[nFind]-dTempX)*(dis[nFind]-dTempX)*(dTempX-dis[nFind-1])/(H[nFind]*H[nFind])
                     -COEX[i][nFind]*(dTempX-dis[nFind-1])*(dTempX-dis[nFind-1])*(dis[nFind]-dTempX)/(H[nFind]*H[nFind])
                    +vY[nFind-1](i)*(dis[nFind]-dTempX)*(dis[nFind]-dTempX)*(2*(dTempX-dis[nFind-1])+H[nFind])/(H[nFind]*H[nFind]*H[nFind])
                    +vY[nFind](i)*(dTempX-dis[nFind-1])*(dTempX-dis[nFind-1])*(2*(dis[nFind]-dTempX)+H[nFind])/(H[nFind]*H[nFind]*H[nFind]);
    }

    /// 释放内存
    CMemPool::GetInstance()->Remove(dis);
    CMemPool::GetInstance()->Remove(H);
    CMemPool::GetInstance()->Remove(A);
    CMemPool::GetInstance()->Remove(B);
    CMemPool::GetInstance()->Remove(C);
    CMemPool::GetInstance()->Remove(D);

    for(i=0; i<nDim; ++i)
    {
        CMemPool::GetInstance()->Remove(COEX[i]);
    }
    CMemPool::GetInstance()->Remove(COEX);

    CMemPool::GetInstance()->Remove(BTA);
    CMemPool::GetInstance()->Remove(CTA);

    return(tmp);
}
