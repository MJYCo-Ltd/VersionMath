#include "VecMat.h"
#include "RK4.h"

using namespace Math;
using namespace Numerical;
CRK4::CRK4(RK4funct pfRK4, int nEqn, void* pAux):
    m_pfRK4(pfRK4),m_nEqn(nEqn),m_pAux(pAux)
{

}

CRK4::~CRK4()
{

}

void CRK4::Step(double& dt, CVector& vY, double dh)
{

    /**
     * 《卫星轨道模型方法和应用》
     * (德)门斯布吕克(Oliver Montenbruck),(德)吉尔(Eberhard Gill)著
     *                   王家松        祝开建        胡小工      译
     * 国防工业出版社,2012.4第1版第1次印刷
     * 113页
     */

    m_pfRK4( dt       , vY               , m_vK1, m_pAux );
    m_pfRK4( dt+dh/2.0, vY+(dh/2.0)*m_vK1, m_vK2, m_pAux );
    m_pfRK4( dt+dh/2.0, vY+(dh/2.0)*m_vK2, m_vK3, m_pAux );
    m_pfRK4( dt+dh    , vY+dh*m_vK3      , m_vK4, m_pAux );

    /// 函数增量
    vY = vY + (dh/6.0)*( m_vK1 + 2.0*m_vK2 + 2.0*m_vK3 + m_vK4 );

    /// 更新自变量
    dt = dt + dh;
}
