#include <iostream>
#include <iomanip>
#include "SatelliteToolKit.h"
#include "sofam.h"
#include "Vector.h"
#include "SGP4.h"

using namespace std;

int main()
{
    string sErrInfo;
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
                                             45*DD2R,45*DD2R,Rectangle);
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

    for(Satellite_Element one : vOut)
    {
        cout<<"---------------------------\ndA:\t"<<one.stKepler.dA<<'\t'
            <<"dE:\t"<<one.stKepler.dE<<'\t'
            <<"dI:\t"<<one.stKepler.dI*DR2D<<'\n'
            <<"dRAAN:\t"<<one.stKepler.dRAAN*DR2D<<'\t'
            <<"dW:\t"<<one.stKepler.dW*DR2D<<'\t'
            <<"dMA:\t"<<one.stKepler.dMA*DR2D<<endl;
    }

    //    if(!SGP4(startTime,endTime,60,sLine1,sLine2,stPos))
    //    {
    //        cout<<"Cal SGP4 faild"<<endl;
    //    }

    //    double dL,dB,dHeight;
    //    Aerospace::CCoorSys::XYZ2LBH(-2941580.75, 5963327.50, 1695900.63,
    //                                 dL,dB,dHeight);
    //    cout<<dL*DR2D<<'\t'<<dB*DR2D<<'\t'<<dHeight<<endl;
    //    PV satPv={-2941580.75, 5963327.50, 1695900.63,
    //             75.5630264, 2130.42725, -7320.25488};

    //    Pos geoPos = {116.3,39.9,0.};
    //    Pos tgtPos = {-2582857.7500000000,5622775.0000000000,1542029.2500000000};

    //    PV satPv={-3668081.0000000000,-1923860.1250000000,6218450.5000000000 ,
    //              -4326.4951171875000,5839.6411132812500,-747.25616455078125};
    //    Pos tgtPos = {-3361188.2500000000,4893817.0000000000,2323193.2500000000};

    //    Pos att = {0.00000000000000000,0.00000000000000000,0.00000000000000000 };
    //    for(int i=0; i<100;++i)
    //    {
    //        cout<<RectangleVisible(satPv,tgtPos,att,0.017453292519943295,0.87266462599716477,Rota_PYR)<<'\t';
    //    }

    return 0;
}
