#include <iostream>
#include <fstream>
#include <iomanip>

#include <VersionMathCommon.h>
#include <Math/Vector.h>
#include <Math/Quaternion.h>

#include <GisMath/GisMath.h>
#include <Satellite/SGP4.h>
#include <Satellite/Date.h>
#include <Satellite/Kepler.h>
#include <SatelliteToolKit/SatelliteToolKit.h>

using namespace std;

int main()
{
    string sErrInfo;

    /// 初始化运算库
    if(!InitSatelliteToolKit(".",sErrInfo))
    {
        cout<<sErrInfo<<endl;
    }

//    CQuaternion cTest(CVector(1,0,0),0);
//    cout<<setprecision(12)<<cTest.GetX()<<'\t'<<cTest.GetY()<<'\t'<<cTest.GetZ()<<'\t'<<cTest.GetS()<<endl;
//    CQuaternion cTest1(CVector(0,1,0),0);
//    cout<<setprecision(12)<<cTest1.GetX()<<'\t'<<cTest1.GetY()<<'\t'<<cTest1.GetZ()<<'\t'<<cTest1.GetS()<<endl;
//    CQuaternion cTest2(CVector(0,0,1),0);
//    cout<<setprecision(12)<<cTest2.GetX()<<'\t'<<cTest2.GetY()<<'\t'<<cTest2.GetZ()<<'\t'<<cTest2.GetS()<<endl;

//    BJTime startTime;
//    startTime.nYear = 2020;
//    startTime.uMonth = 7;
//    startTime.uDay = 1;
//    startTime.uHour = 16;
//    startTime.uMinute=0;
//    startTime.uSecond = 0;
//    startTime.uMSecond=0;
//    BJTime endTime;

//    endTime.nYear = 2020;
//    endTime.uMonth = 7;
//    endTime.uDay = 2;
//    endTime.uHour = 0;
//    endTime.uMinute=0;
//    endTime.uSecond = 0;
//    endTime.uMSecond=0;

//    Pos goundPos={116.3,39.9,0.0};
//    Pos satPYR={0,0,0};

    Satellite_Element tle;
    tle.elemType = SAT_TLE;
    tle.stTLE.sLine1 = "1 90004U 10046A   20061.16666667  .00000000  00000-0  00000-0 0 00000";
    tle.stTLE.sLine2 = "2 90004 123.0013 001.3967 0006000 135.4818 255.2335 13.41413563139433";
//    vector<Period> tmpPeriod = VisiblePeriod(startTime,endTime,tle,goundPos,Rota_PRY,satPYR,
//                                             45*DD2R,45*DD2R,eRectangle);
//    int nIndex(0),nSize(tmpPeriod.size());

//    for(;nIndex < nSize;++nIndex)
//    {
//        const Period& one = tmpPeriod[nIndex];
//        cout<<setw(4)
//           <<nIndex
//           <<" Start:"
//           << one.stStart.nYear
//           <<'-'
//          <<setw(2)
//         <<setfill('0')
//        <<(int)one.stStart.uMonth
//        <<'-'
//        <<setw(2)
//        <<setfill('0')
//        <<(int)one.stStart.uDay
//        <<' '
//        <<setw(2)
//        <<setfill('0')
//        <<(int)one.stStart.uHour
//        <<':'
//        <<setw(2)
//        <<setfill('0')
//        <<(int)one.stStart.uMinute
//        <<':'
//        <<setw(2)
//        <<setfill('0')
//        <<(int)one.stStart.uSecond
//        <<'.'
//        <<setw(3)
//        <<setfill('0')
//        <<one.stStart.uMSecond;
//        cout<<" End:"
//           << one.stEnd.nYear
//           <<'-'
//          <<setw(2)
//          <<setfill('0')
//         <<(int)one.stEnd.uMonth
//        <<'-'
//        <<setw(2)
//        <<setfill('0')
//        <<(int)one.stEnd.uDay
//        <<' '
//        <<setw(2)
//        <<setfill('0')
//        <<(int)one.stEnd.uHour
//        <<':'
//        <<setw(2)
//        <<setfill('0')
//        <<(int)one.stEnd.uMinute
//        <<':'
//        <<setw(2)
//        <<setfill('0')
//        <<(int)one.stEnd.uSecond
//        <<'.'<<setw(3)
//        <<setfill('0')
//        <<one.stEnd.uMSecond;
//        cout<<" Duration:"<<setprecision(6)<<one.dDurationTime<<endl;
//    }


    /// 根据SGP4生成六根数
    Satellite::CSGP4 tmpSGP4(tle.stTLE.sLine1,tle.stTLE.sLine2);

    Math::CVector vKepler = tmpSGP4.ClassicalElements();
//    Satellite_Element tmpElement;
//    tmpElement.elemType = SAT_TWOBODY;
//    tmpElement.stKepler.dA = vKepler(0);
//    tmpElement.stKepler.dE = vKepler(1);
//    tmpElement.stKepler.dI = vKepler(2);
//    tmpElement.stKepler.dRAAN = vKepler(3);
//    tmpElement.stKepler.dW = vKepler(4);
//    tmpElement.stKepler.dMA = vKepler(5);
//    vector<Satellite_Element> vOut = CreateConstellatory(tmpElement,5,6);

//    CDate dataBJ(tmpSGP4.GetTLEEpoch(),BJ);

//    int nYear,nMonth,nDay,nHour,nMinute;
//    double dSeconds;
//    dataBJ.GetDate(nYear,nMonth,nDay,nHour,nMinute,dSeconds);
//    cout<<dataBJ<<endl;
//    for(int nIndex=0,nSize=vOut.size();nIndex<nSize;++nIndex)
//    {
//        const Satellite_Element& one = vOut[nIndex];
//        cout<<"---------------------------\ndA:\t"<<one.stKepler.dA<<'\t'
//            <<"dE:\t"<<one.stKepler.dE<<'\t'
//            <<"dI:\t"<<one.stKepler.dI*DR2D<<'\n'
//            <<"dRAAN:\t"<<one.stKepler.dRAAN*DR2D<<'\t'
//            <<"dW:\t"<<one.stKepler.dW*DR2D<<'\t'
//            <<"dMA:\t"<<one.stKepler.dMA*DR2D<<endl;
//    }


    cout<<"=============================\ndA:\t"<<vKepler(0)<<'\t'
        <<"dE:\t"<<vKepler(1)<<'\t'
        <<"dI:\t"<<vKepler(2)*DR2D<<'\n'
        <<"dRAAN:\t"<<vKepler(3)*DR2D<<'\t'
        <<"dW:\t"<<vKepler(4)*DR2D<<'\t'
        <<"dMA:\t"<<vKepler(5)*DR2D<<endl;
    cout<<"DMJD:\t"<<tmpSGP4.GetTLEEpoch()<<endl;

    char tleOut[2][73];
    Satellite::CKepler::Classical2TLE(vKepler,tmpSGP4.GetTLEEpoch()+DJM0,99,tleOut);
    cout<<tleOut[0]<<tleOut[1]<<endl;

    tmpSGP4.SetTLE(tleOut[0],tleOut[1]);

    vKepler = tmpSGP4.ClassicalElements();
    cout<<"=============================\ndA:\t"<<vKepler(0)<<'\t'
        <<"dE:\t"<<vKepler(1)<<'\t'
        <<"dI:\t"<<vKepler(2)*DR2D<<'\n'
        <<"dRAAN:\t"<<vKepler(3)*DR2D<<'\t'
        <<"dW:\t"<<vKepler(4)*DR2D<<'\t'
        <<"dMA:\t"<<vKepler(5)*DR2D<<endl;
    cout<<"DMJD:\t"<<tmpSGP4.GetTLEEpoch()<<endl;

//     ifstream inputFile("gps-ops.txt");

//     char sBuffer[75];

//     vector<Satellite_Element> vAllSatellite;
//     Satellite_Element tmpSatellite;
//     tmpSatellite.elemType = SAT_TLE;

//     while(inputFile)
//     {
//         inputFile.getline(sBuffer,sizeof(sBuffer));
//         inputFile.getline(sBuffer,sizeof(sBuffer));
//         tmpSatellite.stTLE.sLine1 = sBuffer;
//         inputFile.getline(sBuffer,sizeof (sBuffer));
//         tmpSatellite.stTLE.sLine2 = sBuffer;
//         if(!tmpSatellite.stTLE.sLine1.empty() && !tmpSatellite.stTLE.sLine2.empty())
//         {
//             vAllSatellite.push_back(tmpSatellite);
//         }
//     }

//     vector<double> vResultPDOP = CalPDOP(vAllSatellite,startTime,endTime,60000,goundPos);

//     for(nIndex=0,nSize=vResultPDOP.size();nIndex<nSize;++nIndex)
//     {
//         cout<<"PDOP:"<<vResultPDOP[nIndex]<<endl;
//     }

//    CVector vPos(121.20115022222222*DD2R,23.15952908814338*DD2R,536.6220480057327);
//    CVector vStation(121.20115022222222*DD2R,23.078261805555556*DD2R,289);

//    CVector vPos3D(3),vStation3D(3),vLocal3D(3);
//    GisMath::LBH2XYZ(vPos,vPos3D);
//    GisMath::LBH2XYZ(vStation,vStation3D);

//    CVector vGlobal(vPos3D-vStation3D);
//    GisMath::GLOBAL2LOCAL(vStation(0),vStation(1),vGlobal,vLocal3D);

//    cout<<setprecision(6)<<setw(12)<<vLocal3D;

//    CDate mjBein(2020,3,2,1,0,0,UTC);
//    Satellite::CSGP4 spg41("1 91001U          20061.66666667 -.00000001  00000-0 -13106-2 0 00008",
//                          "2 91001 045.0073 000.0048 0004655 268.5152 091.4846 07.15404217000017"),
//            sgp42("1 91004U          20061.66666667 -.00000001  00000-0 -28120-2 0 00003",
//                  "2 91004 045.0073 180.0049 0004655 268.5153 091.4845 07.15404212000015");

//    double dMJD = mjBein.GetMJD();
//    bool bIsVisible=false;
//    bool bChanged=false;

//    Pos pv1,pv2;

//    CVector vPV;
//    int nIndex=1;
//    for(int i=0;i<86400;++i)
//    {
//        vPV = spg41.CalPV(dMJD+i*SECDAY);
//        pv1.dX = vPV(0);
//        pv1.dY = vPV(1);
//        pv1.dZ = vPV(2);

//        vPV = sgp42.CalPV(dMJD+i*SECDAY);
//        pv2.dX = vPV(0);
//        pv2.dY = vPV(1);
//        pv2.dZ = vPV(2);

//        if(bChanged != (bIsVisible = !InsertEarth(pv1,pv2)))
//        {
//            bChanged = bIsVisible;
//            if(bChanged)
//            {
//                cout<<nIndex<<'\t'<<CDate(dMJD+i*SECDAY,UTC)<<'\t';
//            }
//            else
//            {
//                cout<<CDate(dMJD+i*SECDAY,UTC)<<endl;
//                ++nIndex;
//            }
//        }
//    }
    GisMath::InitGis(GisMath::WGS_84);

    Math::CVector vEye(3);
    Math::CVector vDir(3);
    Math::CVector vOut(3);

    vDir(0) = tan(1*DD2R);
    vDir(1) = vDir(0);
    vDir(2) = 1;
    if(GisMath::CalLineInterEarth(vEye,vDir,vOut))
    {
        cout<<setw(20)<<vOut<<endl;
    }
    double dEarth_a(6378137);
    double dEarth_b(6356752.314245179);
    vDir(0) = 0.;
    vDir(1) = 0.;
    vDir(2) = 1.;
    cout<<dEarth_b<<endl;
    if(GisMath::CalLineInterEarth(vEye,vDir,vOut))
    {
        cout<<setw(20)<<vOut<<endl;
    }



    /// 释放资源
    CloseSatelliteToolKit();
    return 0;
}
