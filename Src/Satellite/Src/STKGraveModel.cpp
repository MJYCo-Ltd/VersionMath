#include <fstream>
#include <cstring>
#include <cstdio>
#include <sstream>
#include <cctype>
#include <cmath>
#include <iomanip>

#include <Math/VecMat.h>
#include <Satellite/STKGraveModel.h>

using namespace std;
using namespace Math;
using namespace Phycical;

/// 防止重定义
#ifndef YTY_STK_GRAVE
#define YTY_STK_GRAVE 1

#define STK_MODEL      "Model"
#define STK_DEGREE     "Degree"
#define STK_ORDER      "Order"
#define STK_GM         "Gm"
#define STK_REFDITANCE "RefDistance"
#define STK_NORMAL     "Normalized"
#define STK_CS_BEGIN   "BEGIN Coefficients"
#define STK_CS_END     "END Coefficients"

#endif

CSTKGraveModel::CSTKGraveModel():m_bReady(false),m_pBuffer(0)
{
}

CSTKGraveModel::CSTKGraveModel(const string &sFileName):m_bReady(false),m_pBuffer(0)
{
    SetFile(sFileName);
}

/// 设置文件名称
bool CSTKGraveModel::SetFile(const string &sFileName)
{
    if(sFileName == m_sFileName)
    {
        return(false);
    }
    else
    {
        m_bReady = false;
        m_sFileName = sFileName;
        return ReadFile();
    }
}

bool CSTKGraveModel::ReadFile()
{
    ifstream inputFile;
    /// 打开文件
    inputFile.open(m_sFileName.c_str());
    if(inputFile.is_open())
    {
        /// 判断文件大小
        streampos sizeBegin = inputFile.tellg();
        inputFile.seekg(0,ios_base::end);
        streampos sizeEnd = inputFile.tellg();
        streamsize size = sizeEnd - sizeBegin;

        /// 一次性读取文件，可以节约时间
        if(size > 0)
        {
            inputFile.seekg(0,ios_base::beg);
            InitBuffer(size);
            inputFile.read(m_pBuffer,size);
        }

        /// 关闭文件
        inputFile.close();

        /// 处理数据
        if(size>0)
        {
            DealBuffer();
        }
    }

    return(m_bReady);
}

void CSTKGraveModel::InitBuffer(const streamsize& nSize)
{
    /// 如果存在数据则清空
    if(0 != m_pBuffer)
    {
        delete []m_pBuffer;
        m_pBuffer = 0;
    }

    /// 开辟空间
    m_pBuffer = new char[nSize + 1]();
}

void CSTKGraveModel::DealBuffer()
{
    if(   ReadStr(m_sModelName,STK_MODEL)
       && ReadInt(m_nDegree,STK_DEGREE)
       && ReadInt(m_nOrder,STK_ORDER)
       && ReadDouble(m_dGM,STK_GM)
       && ReadDouble(m_dR,STK_REFDITANCE)
       && ReadCS()
          )
    {
        m_bReady = true;
    }

}

/// 从字符串中读取一个字符串
bool CSTKGraveModel::ReadStr(string &rsName, const char *cstrSign)
{
    const char* pPos = strstr(m_pBuffer,cstrSign);
    if(NULL != pPos)
    {
        char str[255] = "";
        int n = sscanf(pPos+strlen(cstrSign),"%s",str);

        if(1 == n)
        {
            rsName = str;
            return(true);
        }
    }

    return(false);
}

/// 从字符串中读取一个 int 值
bool CSTKGraveModel::ReadInt(int& nInt, const char* cstrSign)
{
    const char* pPos = strstr(m_pBuffer,cstrSign);
    if(NULL != pPos)
    {
        int n = sscanf(pPos+strlen(cstrSign),"%d",&nInt);

        if(1 == n)
        {
            return(true);
        }
    }

    return(false);
}

/// 从字符串中读取一个 double 值
bool CSTKGraveModel::ReadDouble(double& dDouble, const char* cstrSign)
{
    const char* pPos = strstr(m_pBuffer,cstrSign);
    if(NULL != pPos)
    {
        int n = sscanf(pPos+strlen(cstrSign),"%lf",&dDouble);

        if(1 == n)
        {
            return(true);
        }
    }

    return(false);
}

double CSTKGraveModel::UnNormaliz(int nN, int nM)
{
    /*
     *    《卫星轨道模型方法和应用》
     *   (德)门斯布吕克(Oliver Montenbruck),(德)吉尔(Eberhard Gill)著
     *                 王家松        祝开建        胡小工      译
     *   国防工业出版社,2012.4第1版第1次印刷
     *   54页
     */

    double dFac=0;

    if(nM == 0)
    {
        dFac = 1.0;
    }
    else
    {
        dFac = 2.0;
    }

    /// 对应的归一化系数
    /// 乘以 (2*n + 1)
    dFac *= double((nN<<1) + 1);

    int nStart = nN - nM;
    int nEnd   = nN + nM;

    /// 乘以 (n-m)!/(n+m)!
    for(int i=nEnd; i>nStart; i--)
    {
        if(i>1)
        {
            dFac /= double(i);
        }
    }

    return sqrt(dFac);
}

bool CSTKGraveModel::ReadCS()
{
    const char* pPos = strstr(m_pBuffer,STK_CS_BEGIN);
    if(0 != pPos)
    {
        /// 跳到正确的位置
        pPos += strlen(STK_CS_BEGIN);
        while(isspace(*pPos)) ++pPos;

        /// 构建归一化系数矩阵
        m_mtInfo.Resize(m_nDegree+1,m_nDegree+1);
        m_mtInfo(0,0) = 1.0;

        /// 构建非归一化系数矩阵
        m_unNorInfo.Resize(m_nDegree+1,m_nDegree+1);
        m_unNorInfo(0,0) = 1.0;

        /// 读取数据
        int nN,nM;
        double dC,dS,dFac;
        istringstream s(pPos);
        while(strcmp(STK_CS_END,pPos) != 0 && s)
        {
            /// 读取数据
            streampos posBegin = s.tellg();
            s>>nN>>nM>>dC>>dS;

            dFac = UnNormaliz(nN,nM);
            /// 设置数据
            m_mtInfo(nN,nM) = dC;
            m_unNorInfo(nN,nM) = dC*dFac;
            if(nM > 0)
            {
                m_mtInfo(nM-1,nN) = dS;
                m_unNorInfo(nM-1,nN)=dS*dFac;
            }

            /// 移动读取位置
            pPos += s.tellg() - posBegin + 1;
        }
    }
    return(true);
}
