#include <iostream>
#include <fstream>
#include <iomanip>
#include "SatelliteToolKit.h"
#include "sofa.h"
#include "Vector.h"
#include "SGP4.h"
#include "Date.h"

using namespace std;

int main()
{
    string sErrInfo;

    /// 初始化运算库
    if(!InitSatelliteToolKit(sErrInfo))
    {
        cout<<sErrInfo<<endl;
    }

    BJTime startTime;
    startTime.nYear = 2020;
    startTime.uMonth = 7;
    startTime.uDay = 1;
    startTime.uHour = 16;
    startTime.uMinute=0;
    startTime.uSecond = 0;
    startTime.uMSecond=0;
    BJTime endTime;

    endTime.nYear = 2020;
    endTime.uMonth = 7;
    endTime.uDay = 2;
    endTime.uHour = 0;
    endTime.uMinute=0;
    endTime.uSecond = 0;
    endTime.uMSecond=0;

    Pos goundPos={116.3,39.9,0.0};
    Pos satPYR={0,0,0};

    Satellite_Element tle;
    tle.elemType = SAT_TLE;
    tle.stTLE.sLine1 = "1 90004U 10046A   20061.16666667  .00000000  00000-0  00000-0 0 00000";
    tle.stTLE.sLine2 = "2 90004 123.0013 001.3967 0006000 135.4818 255.2335 13.41413563139433";
    vector<Period> tmpPeriod = VisiblePeriod(startTime,endTime,tle,goundPos,Rota_PRY,satPYR,
                                             45*DD2R,45*DD2R,eRectangle);
    int nIndex=0;

    for(auto one : tmpPeriod)
    {
        cout<<setw(4)
           <<++nIndex
           <<" Start:"
           << one.stStart.nYear
           <<'-'
          <<setw(2)
         <<setfill('0')
        <<(int)one.stStart.uMonth
        <<'-'
        <<setw(2)
        <<setfill('0')
        <<(int)one.stStart.uDay
        <<' '
        <<setw(2)
        <<setfill('0')
        <<(int)one.stStart.uHour
        <<':'
        <<setw(2)
        <<setfill('0')
        <<(int)one.stStart.uMinute
        <<':'
        <<setw(2)
        <<setfill('0')
        <<(int)one.stStart.uSecond
        <<'.'
        <<setw(3)
        <<setfill('0')
        <<one.stStart.uMSecond;
        cout<<" End:"
           << one.stEnd.nYear
           <<'-'
          <<setw(2)
          <<setfill('0')
         <<(int)one.stEnd.uMonth
        <<'-'
        <<setw(2)
        <<setfill('0')
        <<(int)one.stEnd.uDay
        <<' '
        <<setw(2)
        <<setfill('0')
        <<(int)one.stEnd.uHour
        <<':'
        <<setw(2)
        <<setfill('0')
        <<(int)one.stEnd.uMinute
        <<':'
        <<setw(2)
        <<setfill('0')
        <<(int)one.stEnd.uSecond
        <<'.'<<setw(3)
        <<setfill('0')
        <<one.stEnd.uMSecond;
        cout<<" Duration:"<<setprecision(6)<<one.dDurationTime<<endl;
    }


    /// 根据SGP4生成六根数
    Satellite::CSGP4 tmpSGP4(tle.stTLE.sLine1,tle.stTLE.sLine2);

    Math::CVector vKepler = tmpSGP4.ClassicalElements();
    Satellite_Element tmpElement;
    tmpElement.elemType = SAT_TWOBODY;
    tmpElement.stKepler.dA = vKepler(0);
    tmpElement.stKepler.dE = vKepler(1);
    tmpElement.stKepler.dI = vKepler(2);
    tmpElement.stKepler.dRAAN = vKepler(3);
    tmpElement.stKepler.dW = vKepler(4);
    tmpElement.stKepler.dMA = vKepler(5);
    vector<Satellite_Element> vOut = CreateConstellatory(tmpElement,5,6);

    CDate dataBJ(tmpSGP4.GetTLEEpoch(),BJ);

    int nYear,nMonth,nDay,nHour,nMinute;
    double dSeconds;
    dataBJ.GetDate(nYear,nMonth,nDay,nHour,nMinute,dSeconds);
    cout<<dataBJ<<endl;
    for(Satellite_Element one : vOut)
    {
        cout<<"---------------------------\ndA:\t"<<one.stKepler.dA<<'\t'
            <<"dE:\t"<<one.stKepler.dE<<'\t'
            <<"dI:\t"<<one.stKepler.dI*DR2D<<'\n'
            <<"dRAAN:\t"<<one.stKepler.dRAAN*DR2D<<'\t'
            <<"dW:\t"<<one.stKepler.dW*DR2D<<'\t'
            <<"dMA:\t"<<one.stKepler.dMA*DR2D<<endl;
    }

     ifstream inputFile("gps-ops.txt");

     char sBuffer[75];

     vector<Satellite_Element> vAllSatellite;
     Satellite_Element tmpSatellite;
     tmpSatellite.elemType = SAT_TLE;

     while(inputFile)
     {
         inputFile.getline(sBuffer,sizeof(sBuffer));
         inputFile.getline(sBuffer,sizeof(sBuffer));
         tmpSatellite.stTLE.sLine1 = sBuffer;
         inputFile.getline(sBuffer,sizeof (sBuffer));
         tmpSatellite.stTLE.sLine2 = sBuffer;
         if(!tmpSatellite.stTLE.sLine1.empty() && !tmpSatellite.stTLE.sLine2.empty())
         {
             vAllSatellite.push_back(tmpSatellite);
         }
     }

     vector<double> vResultPDOP = CalPDOP(vAllSatellite,startTime,endTime,60.0,goundPos);

     for(auto one : vResultPDOP)
     {
         cout<<"PDOP:"<<one<<endl;
     }

    /// 释放资源
    CloseSatelliteToolKit();
    return 0;
}
