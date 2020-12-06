#include <cmath>
#include <VersionMathCommon.h>
#include <Math/VecMat.h>
#include <GisMath/GisMath.h>
#include <SatelliteToolKit/SatelliteToolKit.h>

/// 计算
vector<double> CalPDOP(const vector<Satellite_Element>& vSatellite,
                      const BJTime& stStartTime,
                      const BJTime& stEndTime,
                      double dDeltaTime,
                      const Pos& rPos)
{
    vector<SatellitePos> vPos;
    vPos.resize(vSatellite.size());

    /// 计算卫星的位置
    int i(0),nSize(vSatellite.size());

    for(;i<nSize;++i)
    {
        const Satellite_Element& one = vSatellite[i];
        switch(one.elemType)
        {
        case SAT_TLE:
            SGP4(stStartTime,stEndTime,dDeltaTime,one.stTLE.sLine1,one.stTLE.sLine2,vPos[i]);
            break;
        case SAT_TWOBODY:
            TwoBody(stStartTime,stEndTime,dDeltaTime,one.stKepler,vPos[i]);
            break;
        case SAT_PV:
            break;
        }
    }

    Pos v3D;
    GisMath::LBH2XYZ(rPos.dX*DD2R, rPos.dY*DD2R, rPos.dZ, v3D.dX, v3D.dY, v3D.dZ);

    vector<double> vPDOP;
    i = vPos[0].vECF.size();
    vPDOP.resize(i);

    Pos vSat3D;

    Math::CVector vStation3D(v3D.dX,v3D.dY,v3D.dZ),v4(4),v3(3);
    v4(3) = -1;

    vector<Math::CVector> vAllVisibleSat;
    vAllVisibleSat.reserve(vPos.size());

    for(int nIndex=0;nIndex<i;++nIndex)
    {
        vAllVisibleSat.clear();
        for(int j=0;j<vPos.size();++j)
        {
            vSat3D = vPos[j].vECF[nIndex].stP;

            /// 只有对地面位置可见的星才参与pdop计算
            if(IsVisible(vSat3D,v3D,5*DD2R))
            {
                Math::CVector vSat(vSat3D.dX,vSat3D.dY,vSat3D.dZ);
                v3 = vSat-vStation3D;
                v3.Normalize();
                v4.Set(v3(0),v3(1),v3(2));
                vAllVisibleSat.push_back(v4);
            }
        }

        /// 如果有超过一个星可见
        if(vAllVisibleSat.size() > 0)
        {
            Math::CMatrix matA(vAllVisibleSat.size(),4);
            for(int j=vAllVisibleSat.size()-1; j>-1; --j)
            {
                matA.SetRow(j,vAllVisibleSat[j]);
            }

            /// 算法参见 https://en.wikipedia.org/wiki/Dilution_of_precision_(navigation)
            Math::CMatrix matResult = Math::CVecMat::Inv(Math::CVecMat::Transp(matA) * matA);

            vPDOP[nIndex] = sqrt(matResult(0,0)+matResult(1,1)+matResult(2,2));
        }
        else
        {
            vPDOP[nIndex]=0;
        }
    }
    return(vPDOP);
}
