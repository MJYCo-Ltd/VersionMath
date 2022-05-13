#include <cmath>
#include <GisMath/GisMath.h>
#include <SatelliteToolKit/SatelliteToolKit.h>

/// 计算两个卫星和地面站的夹角
std::vector<double> Cal2SatGroundAngle(const std::vector<PV>& vAx,
                                  const std::vector<PV>& vBx,
                                  const Pos& rGroundStation)
{
    std::vector<double> vTemp;

    if(vAx.size() != vBx.size())
    {
        return(vTemp);
    }

    vTemp.resize(vAx.size());

    Math::CVector v3D(3);

    GisMath::LBH2XYZ(rGroundStation.dX,rGroundStation.dY,rGroundStation.dZ,
                     v3D(0),v3D(1),v3D(2));
    Math::CMatrix g2l = GisMath::GLOBAL2LOCAL(rGroundStation.dX,rGroundStation.dY);

    int nSize = vAx.size();
    Math::CVector vSatA(3),vSatB(3);
    double dElevA,dAzimA,dElevB,dAzimB,dDotAB;
    for(int nIndex=0; nIndex<nSize; ++nIndex)
    {
        const Pos &rSatA = vAx.at(nIndex).stP;
        const Pos &rSatB = vBx.at(nIndex).stP;

        vSatA.Set(rSatA.dX - v3D(0),
                  rSatA.dY - v3D(1),
                  rSatA.dZ - v3D(2));

        vSatA = g2l * vSatA;

        vSatB.Set(rSatB.dX - v3D(0),
                  rSatB.dY - v3D(1),
                  rSatB.dZ - v3D(2));

        vSatB = g2l * vSatB;

        Math::CVecMat::AzEl(vSatA,dAzimA,dElevA);
        Math::CVecMat::AzEl(vSatB,dAzimB,dElevB);

        /// 只有当两颗星都可见的时候才返回夹角
        if(dElevA >0 && dElevB>0)
        {
            dDotAB = Math::CVecMat::Dot(vSatA,vSatB);
            vTemp[nIndex] = acos(dDotAB/vSatA.Length()/vSatB.Length());
        }
        else
        {
            vTemp[nIndex]=-1;
        }
    }

    return(vTemp);
}

/// 计算某个卫星和同一个地面站的夹角
std::vector<double>CalSat2GroundAngle(const std::vector<PV>& vAx,
                                 const Pos& rGoundA,
                                 const Pos& rGoundB)
{
    std::vector<double> vTemp;
    if(vAx.size() < 1)
    {
        return (vTemp);
    }

    vTemp.resize(vAx.size());

    Math::CVector v3DA(3),v3DB(3);

    GisMath::LBH2XYZ(rGoundA.dX,rGoundA.dY,rGoundA.dZ,
                     v3DA(0),v3DA(1),v3DA(2));
    GisMath::LBH2XYZ(rGoundB.dX,rGoundB.dY,rGoundB.dZ,
                     v3DB(0),v3DB(1),v3DB(2));
    Math::CMatrix g2lA = GisMath::GLOBAL2LOCAL(rGoundA.dX,rGoundA.dY);
    Math::CMatrix g2lB = GisMath::GLOBAL2LOCAL(rGoundB.dX,rGoundB.dY);

    int nSize = vAx.size();
    Math::CVector vGroundA(3),vGroundB(3),vSat_GroudA(3),vSat_GroundB(3);
    double dElevA,dAzimA,dElevB,dAzimB,dDotAB;
    for(int nIndex=0; nIndex<nSize; ++nIndex)
    {
        const Pos &rSat = vAx.at(nIndex).stP;

        vGroundA.Set(rSat.dX - v3DA(0),
                     rSat.dY - v3DA(1),
                     rSat.dZ - v3DA(2));

        vGroundA = g2lA * vGroundA;

        vGroundB.Set(rSat.dX - v3DB(0),
                     rSat.dY - v3DB(1),
                     rSat.dZ - v3DB(2));

        vGroundB = g2lB * vGroundB;

        Math::CVecMat::AzEl(vGroundA,dAzimA,dElevA);
        Math::CVecMat::AzEl(vGroundB,dAzimB,dElevB);

        /// 只有当两颗星都可见的时候才返回夹角
        if(dElevA >0 && dElevB>0)
        {
            vSat_GroudA.Set(v3DA(0)-rSat.dX,
                            v3DA(1)-rSat.dY,
                            v3DA(2)-rSat.dZ);

            vSat_GroundB.Set(v3DB(0)-rSat.dX,
                             v3DB(1)-rSat.dY,
                             v3DB(2)-rSat.dZ);
            dDotAB = Math::CVecMat::Dot(vSat_GroudA,vSat_GroundB);
            vTemp[nIndex] = acos(dDotAB/vSat_GroudA.Length()/vSat_GroundB.Length());
        }
        else
        {
            vTemp[nIndex]=-1;
        }
    }

    return(vTemp);
}
