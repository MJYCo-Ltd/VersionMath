#include <cmath>
#include <limits>
#include <Math/VecMat.h>
#include <Math/Quaternion.h>
using namespace std;
using namespace Math;

CQuaternion::CQuaternion()
{
    InitQuat();
}

/// 复制构造
CQuaternion::CQuaternion(const CQuaternion &rOther)
{
    m_dX = rOther.m_dX;
    m_dY = rOther.m_dY;
    m_dZ = rOther.m_dZ;
    m_dS = rOther.m_dS;
}

CQuaternion::CQuaternion(double dx, double dy, double dz, double ds):
    m_dX(dx),m_dY(dy),m_dZ(dz),m_dS(ds)
{
    Normalize();
    InitCoefficient();
}

/// 通过向量与旋转角度构建四元数
CQuaternion::CQuaternion(const CVector &vVec, double dphi)
{
    Rebuild(vVec,dphi);
}

/// 通过向量构建四元数
CQuaternion::CQuaternion(const CVector &vVec)
{
    if(vVec.Size() < 4)
    {
        InitQuat();
        return;
    }

    /// 复制值
    m_dX = vVec(0);
    m_dY = vVec(1);
    m_dZ = vVec(2);
    m_dS = vVec(3);
}

/// 通过方向余弦矩阵构造
CQuaternion::CQuaternion(const CMatrix &vMat)
{
    Rebuild(vMat);
}

/// 通过矩阵构造
bool CQuaternion::Rebuild(const CMatrix &vMat)
{
    if(vMat.Row()!=3 || vMat.Col()!=3)
    {
        InitQuat();
        return(false);
    }

    /// 临时变量
    double ds;
    double dTq[4];
    int    i, j;

    // Use tq to store the largest trace
    dTq[0] = 1 + vMat(0,0)+vMat(1,1)+vMat(2,2);
    dTq[1] = 1 + vMat(0,0)-vMat(1,1)-vMat(2,2);
    dTq[2] = 1 - vMat(0,0)+vMat(1,1)-vMat(2,2);
    dTq[3] = 1 - vMat(0,0)-vMat(1,1)+vMat(2,2);

    /// 找到最大值
    j = 0;
    for(i=1;i<4;++i)
    {
        j = (dTq[i]>dTq[j])? i : j;
    }

    // check the diagonal
    if (j==0)
    {
        /* perform instant calculation */
        m_dS = dTq[0];
        m_dX = vMat(1,2)-vMat(2,1);
        m_dY = vMat(2,0)-vMat(0,2);
        m_dZ = vMat(0,1)-vMat(1,0);
    }
    else if (j==1)
    {
        m_dS = vMat(1,2)-vMat(2,1);
        m_dX = dTq[1];
        m_dY = vMat(0,1)+vMat(1,0);
        m_dZ = vMat(2,0)+vMat(0,2);
    }
    else if (j==2)
    {
        m_dS = vMat(2,0)-vMat(0,2);
        m_dX = vMat(0,1)+vMat(1,0);
        m_dY = dTq[2];
        m_dZ = vMat(1,2)+vMat(2,1);
    }
    else /* if (j==3) */
    {
        m_dS = vMat(0,1)-vMat(1,0);
        m_dX = vMat(2,0)+vMat(0,2);
        m_dY = vMat(1,2)+vMat(2,1);
        m_dZ = dTq[3];
    }

    ds = sqrt(0.25/dTq[j]);
    m_dS *= ds;
    m_dX *= ds;
    m_dY *= ds;
    m_dZ *= ds;

    m_dX1 = vMat(0,0);
    m_dY1 = vMat(0,1);
    m_dZ1 = vMat(0,2);
    m_dX2 = vMat(1,0);
    m_dY2 = vMat(1,1);
    m_dZ2 = vMat(1,2);
    m_dX3 = vMat(2,0);
    m_dY3 = vMat(2,1);
    m_dZ3 = vMat(2,2);

    m_bNeedInitCoefficient=false;
    return(true);
}

bool CQuaternion::Rebuild(const CVector &vVec, double dphi)
{
    if(vVec.Size() < 3)
    {
        InitQuat();
        return(false);
    }

    if(vVec.IsNormal())
    {
        double sphi = sin(dphi/2.0);
        m_dX = vVec.GetX() * sphi;
        m_dY = vVec.GetY() * sphi;
        m_dZ = vVec.GetZ() * sphi;
        m_dS = cos(dphi/2.0);
    }
    else
    {
        /// 复制值
        CVector vTmp = vVec;

        /// 将向量归一化
        vTmp.Normalize();

        /// 进行旋转
        double sphi = sin(dphi/2.0);
        m_dX = vTmp.GetX() * sphi;
        m_dY = vTmp.GetY() * sphi;
        m_dZ = vTmp.GetZ() * sphi;
        m_dS = cos(dphi/2.0);
    }

    m_bNeedInitCoefficient=true;
    return(true);
}

CQuaternion::~CQuaternion()
{
}

void CQuaternion::InitQuat()
{
    m_dX = m_dY = m_dZ = 0.;
    m_dS = 1.0;
}

/// 初始化系数
void CQuaternion::InitCoefficient()
{
    m_dX1=1.0 - 2.0*(m_dY*m_dY+m_dZ*m_dZ);
    m_dY1=2.0*(m_dX*m_dY+m_dZ*m_dS);
    m_dZ1=2.0*(m_dX*m_dZ-m_dY*m_dS);
    m_dX2=2.0*(m_dX*m_dY-m_dZ*m_dS);
    m_dY2=1.0 - 2.0*(m_dX*m_dX+m_dZ*m_dZ);
    m_dZ2=2.0*(m_dY*m_dZ+m_dX*m_dS);
    m_dX3=2.0*(m_dX*m_dZ+m_dY*m_dS);
    m_dY3=2.0*(m_dY*m_dZ-m_dX*m_dS);
    m_dZ3=1.0 - 2.0*(m_dX*m_dX+m_dY*m_dY);
}

/// 将四元数归一化
void CQuaternion::Normalize()
{
    double dNormal = Norm();
    if(0. == dNormal)
    {
        return;
    }
    else
    {
        double dTmp = 1.0 / dNormal;
        m_dX *= dTmp;
        m_dY *= dTmp;
        m_dZ *= dTmp;
        m_dS *= dTmp;
    }
}

/// 四元数的模
double CQuaternion::Norm() const
{
    return(sqrt(m_dX*m_dX + m_dY*m_dY + m_dZ*m_dZ + m_dS*m_dS));
}


CMatrix CQuaternion::GetMatrix()
{
    CMatrix matTmp(3,3);

    double dLength = Norm();
    if(fabs(dLength) <= numeric_limits<double>::min())
    {
        return(matTmp);
    }

    if(m_bNeedInitCoefficient)
    {
        Normalize();
        InitCoefficient();
        m_bNeedInitCoefficient = false;
    }

    /// 构建旋转矩阵
    matTmp(0,0) = m_dX1;
    matTmp(0,1) = m_dY1;
    matTmp(0,2) = m_dZ1;
    matTmp(1,0) = m_dX2;
    matTmp(1,1) = m_dY2;
    matTmp(1,2) = m_dZ2;
    matTmp(2,0) = m_dX3;
    matTmp(2,1) = m_dY3;
    matTmp(2,2) = m_dZ3;

    return (matTmp);
}

CVector CQuaternion::GetVector() const
{
    CVector vecTmp(4);

    vecTmp(0) = m_dX;
    vecTmp(1) = m_dY;
    vecTmp(2) = m_dZ;
    vecTmp(3) = m_dS;

    return(vecTmp);
}

void CQuaternion::GetRotate(CVector &rVec, double &dAngle)const
{
    /// 如果数据有效且数据不等于3，则直接返回
    if(rVec && rVec.Size() != 3)
    {
        return;
    }

    rVec.Resize(3);

    double dHalfAngle = sqrt( m_dX*m_dX + m_dY*m_dY + m_dZ*m_dZ );

    /// 计算角度
    dAngle = 2.0 * atan2( dHalfAngle, m_dS );
    if(dHalfAngle)
    {
        rVec(0) = m_dX / dHalfAngle;
        rVec(1) = m_dY / dHalfAngle;
        rVec(2) = m_dZ / dHalfAngle;
    }
    else
    {
        rVec(0) = 0.0;
        rVec(1) = 0.0;
        rVec(2) = 1.0;
    }
}

/// 完成负号操作符
CQuaternion& CQuaternion::operator -()
{
    this->operator *=(-1.);
    return(*this);
}

CQuaternion& CQuaternion::operator +=(const CQuaternion& rOther)
{
    m_dX += rOther.m_dX;
    m_dY += rOther.m_dY;
    m_dZ += rOther.m_dZ;
    m_dS += rOther.m_dS;

    return(*this);
}

CQuaternion& CQuaternion::operator -=(const CQuaternion& rOther)
{
    m_dX -= rOther.m_dX;
    m_dY -= rOther.m_dY;
    m_dZ -= rOther.m_dZ;
    m_dS -= rOther.m_dS;

    return(*this);
}

CQuaternion& CQuaternion::operator *=(double dScal)
{
    m_dX *= dScal;
    m_dY *= dScal;
    m_dZ *= dScal;
    m_dS *= dScal;

    return(*this);
}

/// *=操作符
CQuaternion& CQuaternion::operator *=(const CQuaternion& rOther)
{
/*
    /// 另一种计算方法
    CVector v1(m_dX,m_dY,m_dZ);
    CVector v2(rOther.m_dX,rOther.m_dY,rOther.m_dZ);
    CVector vOut;

    double w = (m_dS * rOther.m_dS) - CVecMat::Dot(v1,v2);

    vOut = CVecMat::Cross(v1, v2);

    v1 *= w;
    v2 *= w;

    vOut += v1;
    vOut += v2;

    m_dX = vOut.GetX();
    m_dY = vOut.GetY();
    m_dZ = vOut.GetZ();
    m_dS = w;

    return(*this);
*/
    /// 计算乘后的结果
    double dX = rOther.m_dS*m_dX + rOther.m_dX*m_dS + rOther.m_dY*m_dZ - rOther.m_dZ*m_dY;
    double dY = rOther.m_dS*m_dY - rOther.m_dX*m_dZ + rOther.m_dY*m_dS + rOther.m_dZ*m_dX;
    double dZ = rOther.m_dS*m_dZ + rOther.m_dX*m_dY - rOther.m_dY*m_dX + rOther.m_dZ*m_dS;
    double dS = rOther.m_dS*m_dS - rOther.m_dX*m_dX - rOther.m_dY*m_dY - rOther.m_dZ*m_dZ;

    m_dX = dX;
    m_dY = dY;
    m_dZ = dZ;
    m_dS = dS;

    m_bNeedInitCoefficient = true;

    return(*this);
}

/// /= 操作符
CQuaternion& CQuaternion::operator /=(const CQuaternion& rOther)
{
    this->operator *= (Inv(rOther));
    m_bNeedInitCoefficient = true;
    return (*this);
}

CQuaternion& CQuaternion::operator /=(double dScal)
{
    /// 防止除0
    if(0.0 == dScal)
    {
        InitQuat();
    }
    else
    {
        dScal = 1.0 / dScal;
        this->operator *=(dScal);
    }
    m_bNeedInitCoefficient = true;

    return(*this);
}

/// 赋值操作符
CQuaternion& CQuaternion::operator =(const CQuaternion& rOther)
{
    if(this == &rOther)
    {
        return(*this);
    }

    /// 赋值值到this指针
    this->m_dX = rOther.m_dX;
    this->m_dY = rOther.m_dY;
    this->m_dZ = rOther.m_dZ;
    this->m_dS = rOther.m_dS;

    m_bNeedInitCoefficient = true;

    return(*this);
}

CQuaternion& CQuaternion::operator =(const CMatrix& rMat)
{
    this->Rebuild(rMat);
    return(*this);
}

/// 判等 操作符
bool CQuaternion::operator ==(const CQuaternion& rOther) const
{

    /// 判断地址是否相同
    if(this == &rOther)
    {
        return(true);
    }

    /// 判断值是否相同
    if(fabs(m_dX - rOther.m_dX) > numeric_limits<double>::min())
    {
        return(false);
    }
    else if(fabs(m_dY - rOther.m_dY) > numeric_limits<double>::min())
    {
        return(false);
    }
    else if(fabs(m_dZ - rOther.m_dZ) > numeric_limits<double>::min())
    {
        return(false);
    }
    else if(fabs(m_dS - rOther.m_dS) > numeric_limits<double>::min())
    {
        return(false);
    }

    /// 值全部相等则返回true
    return(true);
}

/// 判断是否与旋转矩阵是否相等
bool CQuaternion::operator ==(const CMatrix& rMat)const
{
    CQuaternion tmp(rMat);
    return(this->operator ==(tmp));
}

bool CQuaternion::Translate(const CVector &vIn, CVector &vOut)
{
    if(vIn.Size() > 2)
    {
        if(m_bNeedInitCoefficient)
        {
            Normalize();
            InitCoefficient();
            m_bNeedInitCoefficient = false;
        }

        unsigned int nCount = vIn.Size() / 3;
        double dvX,dvY,dvZ;

        /// 如果输出大小不足以存储则扩展空间
        if(vOut.Size() < nCount*3)
        {
            vOut.Resize(nCount*3);
        }

        for(unsigned int i=0,j=0;i<nCount;++i,j+=3)
        {
            dvX = vIn(j);
            dvY = vIn(j+1);
            dvZ = vIn(j+2);
            vOut(j)   = m_dX1*dvX + m_dY1*dvY + m_dZ1*dvZ;
            vOut(j+1) = m_dX2*dvX + m_dY2*dvY + m_dZ2*dvZ;
            vOut(j+2) = m_dX3*dvX + m_dY3*dvY + m_dZ3*dvZ;
        }

        return(true);
    }
    else
    {
        return(false);
    }
}

CQuaternion::operator bool()const
{
    return(m_dX == 0. ? m_dY == 0. ? m_dZ == 0. ? m_dS == 0. ? false : true : true : true : true);
}

////////////////////////以下非成员函数///////////////////////////////////////
CQuaternion operator +(const CQuaternion& qOne, const CQuaternion& qTwo)
{
    CQuaternion tmp(qOne);
    tmp += qTwo;
    return(tmp);
}

CQuaternion operator -(const CQuaternion& qOne, const CQuaternion& qTwo)
{
    CQuaternion tmp(qOne);
    tmp -= qTwo;
    return(tmp);
}

CQuaternion operator * (const CQuaternion& q, double dScal)
{
    CQuaternion tmp(q);
    tmp *= dScal;
    return(tmp);
}

CQuaternion operator * (double dScal ,const CQuaternion& q)
{
    CQuaternion tmp(q);
    tmp *= dScal;
    return(tmp);
}

CQuaternion operator * (const CQuaternion& rOne, const CQuaternion& rTwo)
{
    CQuaternion q(rOne);
    q *= rTwo;
    return(q);
}

CQuaternion operator / (const CQuaternion& rOne, const CQuaternion& rTwo)
{
    CQuaternion q(rOne);
    q /= rTwo;
    return(q);
}

CVector operator * (const CQuaternion& q,const CVector& v)
{
    /*
    /// 另一种算法 来自Nvidia SDK
    CVector uv, uuv;
    CVector qvec(q.GetX(), q.GetY(), q.GetZ());
    uv = CVecMat::Cross(qvec,v);
    uuv = CVecMat::Cross(qvec,uv);
    uv *= ( 2.0 * q.GetS() );
    uuv *= 2.0;
    return v + uv + uuv;
    */
    double dqX = q.GetX(),
           dqY = q.GetY(),
           dqZ = q.GetZ(),
           dqS = q.GetS();

    double dvX = v.GetX(),
           dvY = v.GetY(),
           dvZ = v.GetZ();
    return CVector( (1.0 - 2.0*(dqY*dqY+dqZ*dqZ))*dvX + 2.0*(dqX*dqY+dqZ*dqS)*dvY + 2.0*(dqX*dqZ-dqY*dqS)*dvZ,
                        2.0*(dqX*dqY-dqZ*dqS)*dvX + (1.0 - 2.0*(dqX*dqX+dqZ*dqZ))*dvY + 2.0*(dqY*dqZ+dqX*dqS)*dvZ,
                        2.0*(dqX*dqZ+dqY*dqS)*dvX + 2.0*(dqY*dqZ-dqX*dqS)*dvY + (1.0 - 2.0*(dqX*dqX+dqY*dqY))*dvZ);
}

CQuaternion Inv(const CQuaternion& rOther)
{
    CQuaternion q(rOther);
    double dNomal = rOther.Norm();
    dNomal = dNomal * dNomal;

    if(dNomal <= numeric_limits<double>::min())
    {
        return(q);
    }

    q /=(-dNomal);

    return(q);
}

CQuaternion Conj(const CQuaternion& rOther)
{
    CQuaternion q(rOther);
    -q;
    return(q);
}

CQuaternion Slerp(const CQuaternion &rqFrom, const CQuaternion &rqTo, double dPos)
{
    /// 防止为负数
    if(dPos < 0.)
    {
        dPos = -dPos;
    }

    /// 防止越界
    if(dPos > 1.0)
    {
        dPos -= (int)dPos;
    }

    CQuaternion q;
    const double epsilon = 0.00001;
    double omega, cosomega, sinomega, scale_from, scale_to ;

    CQuaternion quatTo(rqTo);

    CVector vFrom = rqFrom.GetVector();
    CVector vTo   = rqTo.GetVector();
    cosomega = CVecMat::Dot(vFrom,vTo);

    /// 如果两个向量夹角大于90度
    if ( cosomega <0.0 )
    {
        cosomega = -cosomega;
        -quatTo;
    }

    /// 如果两个差别很大
    if( (1.0 - cosomega) > epsilon )
    {
        omega= acos(cosomega) ;
        sinomega = sin(omega) ;

        scale_from = sin((1.0-dPos)*omega)/sinomega ;
        scale_to = sin(dPos*omega)/sinomega ;
    }
    else
    {
       /// 如果两个向量基本上重合
       scale_from = 1.0 - dPos ;
       scale_to = dPos ;
    }

    q = (rqFrom*scale_from) + (quatTo*scale_to);
    return(q);
}
