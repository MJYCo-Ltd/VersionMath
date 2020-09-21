#include <iostream>
#include <fstream>
#include <iomanip>
#include "SatelliteToolKit.h"
#include "sofa.h"
#include "Vector.h"
#include "SGP4.h"
#include "Date.h"
#include "Quaternion.h"

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

//    Satellite_Element tle;
//    tle.elemType = SAT_TLE;
//    tle.stTLE.sLine1 = "1 90004U 10046A   20061.16666667  .00000000  00000-0  00000-0 0 00000";
//    tle.stTLE.sLine2 = "2 90004 123.0013 001.3967 0006000 135.4818 255.2335 13.41413563139433";
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


//    /// 根据SGP4生成六根数
//    Satellite::CSGP4 tmpSGP4(tle.stTLE.sLine1,tle.stTLE.sLine2);

//    Math::CVector vKepler = tmpSGP4.ClassicalElements();
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
//    for(nIndex=0,nSize=vOut.size();nIndex<nSize;++nIndex)
//    {
//        const Satellite_Element& one = vOut[nIndex];
//        cout<<"---------------------------\ndA:\t"<<one.stKepler.dA<<'\t'
//            <<"dE:\t"<<one.stKepler.dE<<'\t'
//            <<"dI:\t"<<one.stKepler.dI*DR2D<<'\n'
//            <<"dRAAN:\t"<<one.stKepler.dRAAN*DR2D<<'\t'
//            <<"dW:\t"<<one.stKepler.dW*DR2D<<'\t'
//            <<"dMA:\t"<<one.stKepler.dMA*DR2D<<endl;
//    }

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

    CDate mjBein(2020,3,2,1,0,0,UTC);
    Satellite::CSGP4 spg41("1 91007U          20061.66666667 -.00000102  00000-0 -65164+3 0 00007",
                          "2 91007 000.9785 124.5021 0000033 253.1263 107.2755 01.00272001000012"),
            sgp42("1 90002U 09055A   20061.16666667  .00000136  00000-0  44938-4 0 00008",
                  "2 90002 098.5033 141.9267 0001830 317.3362 227.2961 14.38039522545503");

    double dMJD = mjBein.GetMJD();
    bool bIsVisible=false;
    bool bChanged=false;

    PV pv1,pv2;

    CVector vPV;
    int nIndex=1;
    for(int i=0;i<86400;++i)
    {
        vPV = spg41.CalPV(dMJD+i*SECDAY);
        pv1.stP.dX = vPV(0);
        pv1.stP.dY = vPV(1);
        pv1.stP.dZ = vPV(2);
        pv1.stV.dX = vPV(3);
        pv1.stV.dY = vPV(4);
        pv1.stV.dZ = vPV(5);

        vPV = sgp42.CalPV(dMJD+i*SECDAY);
        pv2.stP.dX = vPV(0);
        pv2.stP.dY = vPV(1);
        pv2.stP.dZ = vPV(2);
        pv2.stV.dX = vPV(3);
        pv2.stV.dY = vPV(4);
        pv2.stV.dZ = vPV(5);

        if(bChanged != (bIsVisible = IsVisible(pv1,pv2,0.26179938779914943653855361527329)))
        {
            bChanged = bIsVisible;
            if(bChanged)
            {
                cout<<nIndex<<'\t'<<CDate(dMJD+i*SECDAY,UTC)<<'\t';
            }
            else
            {
                cout<<CDate(dMJD+i*SECDAY,UTC)<<endl;
                ++nIndex;
            }
        }
    }
    /// 释放资源
    CloseSatelliteToolKit();
    return 0;
}
