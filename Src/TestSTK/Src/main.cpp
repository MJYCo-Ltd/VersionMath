#include <cmath>
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
#include <Satellite/JPLEphemeris.h>
#include <Satellite/TimeSys.h>
#include <Satellite/CoorSys.h>
#include <Satellite/J2.h>

using namespace std;
int main()
{
    string sErrInfo;

    /// 初始化运算库
    if(!InitSatelliteToolKit(".",sErrInfo))
    {
        cout<<sErrInfo<<endl;
    }
    GisMath::InitGis(GisMath::WGS_84);

//    for(int k=0;k<100;++k)
//    {
//        CVector vTest(6);
//        vector<CVector> vPushData;
//        int nNum(10000);
//        /// 通过遍历
//        for(int i=0;i<nNum;++i)
//        {
//            for(int j=nNum-i;j>0;--j)
//            {
//                vPushData.push_back(vTest);
//            }
//            vPushData.clear();
//        }
//        cout<<"K="<<k<<endl;
//    }
//    return(0);

//    Math::CVector vPos(384263.875,7061764.50,2416784.75);
//    cout<<"Distance: "<<vPos.Length()<<" Earth R:"<<R_Earth2<<endl;
//    double dSinA = asin(R_Earth2 / vPos.Length());
//    cout<<"Sin Angle: "<<dSinA*DR2D<<endl;
//    /// 向内0.15
//    dSinA -= 0.15* DD2R;
//    Math::CVector vZPlane(0,0,1);

//    /// 计算一个旋转轴
//    Math::CVector vRotateAix = CVecMat::Cross(vPos,vZPlane);
//    Math::CQuaternion vRotateQua(vRotateAix,dSinA);
//    Math::CVector vRotatedPos = vRotateQua * -vPos;

//    /// 判断两个向量的夹角是否在Sin角度内
//    double dDot = CVecMat::Dot(vRotatedPos,-vPos);
//    double dCosB = acos(dDot/vRotatedPos.Length()/vPos.Length());
//    cout<<"Cos Angle: "<<dCosB*DR2D<<endl;

    /// 计算卫星最近碰撞距离
    if(0)
    {
        Math::CVector vPos(378617,7.06624e+06,2.40458e+06),vRotatedPos(-0.886557,-0.374096,-0.306347),vInsertPos;
        double dDot = CVecMat::Dot(vRotatedPos,-vPos);
        double dCosB = acos(dDot/vRotatedPos.Length()/vPos.Length());
        cout<<"Cos Angle: "<<dCosB*DR2D<<endl;
        if(!GisMath::CalLineInterEarth(vPos,vRotatedPos,vInsertPos))
        {
            cout<<"not Insert"<<endl;
        }
        else
        {
            cout<<setprecision(12)<<vInsertPos<<endl;
        }
        cout<<"call end"<<endl;

        Satellite_Element tle;
        tle.elemType = SAT_TLE;
        tle.stTLE.sLine1 = "1 90101U          20061.16666667  .00000000   00000- +00000-0 0   872";
        tle.stTLE.sLine2 = "2 90101  98.5033 141.9267 0001830 317.3362 227.2961 14.38038875-1795754448Z";

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

        cout<<"A:\t"<<tmpElement.stKepler.dA<<endl;
        cout<<"E:\t"<<tmpElement.stKepler.dE<<endl;
        cout<<"I:\t"<<tmpElement.stKepler.dI<<endl;
        cout<<"RAAN:\t"<<tmpElement.stKepler.dRAAN<<endl;
        cout<<"W:\t"<<tmpElement.stKepler.dW<<endl;
        cout<<"MA:\t"<<tmpElement.stKepler.dMA<<endl;

        BJTime nowBJTime;
        nowBJTime.nYear = 2021;
        nowBJTime.uMonth = 9;
        nowBJTime.uDay = 1;
        nowBJTime.uHour = 17;
        nowBJTime.uMinute = 43;
        nowBJTime.uSecond = 0;
        nowBJTime.uMSecond = 0;
        tmpElement.stKepler.stEpoch = nowBJTime;

        Kepler oribit = InsertOribit(tmpElement,90*DD2R,nowBJTime);
        cout<<"A:\t"<<oribit.dA<<endl;
        cout<<"E:\t"<<oribit.dE<<endl;
        cout<<"I:\t"<<oribit.dI<<endl;
        cout<<"RAAN:\t"<<oribit.dRAAN<<endl;
        cout<<"W:\t"<<oribit.dW<<endl;
        cout<<"MA:\t"<<oribit.dMA<<endl;

//        tmpElement.stKepler.dA = R_Earth;
//        tmpElement.stKepler.dE = 0;
//        tmpElement.stKepler.dI = 0;
//        tmpElement.stKepler.dRAAN = 0;
//        tmpElement.stKepler.dW = 0;
//        tmpElement.stKepler.dMA = 0;
    }

    /// 计算面积
    if(0)
    {
        Math::CVector vCalArea(109.006,38.3493,113.45,39.703,114.27,34.9285);
        cout<<GisMath::CalclutaGeoArea(vCalArea,GisMath::Msqu)<<endl;
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

    if(0)
    {
        Math::CVector vKepler;
        Satellite_Element tle;
        tle.elemType = SAT_TLE;
        tle.stTLE.sLine1 = "1 90101U          20061.16666667  .00000000   00000- +00000-0 0   872";
        tle.stTLE.sLine2 = "2 90101  98.5033 141.9267 0001830 317.3362 227.2961 14.38038875-1795754448Z";
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

        //    cout<<tmpSGP4<<endl;
        //    for(int i=0;i<100;++i)
        //    {
        //        cout<<setprecision(6)<<tmpSGP4.CalPV(58909.0000000000000000)<<endl;
        //    }

        vKepler = tmpSGP4.ClassicalElements();
        Satellite_Element tmpElement;
        tmpElement.elemType = SAT_TWOBODY;
        tmpElement.stKepler.dA = vKepler(0);
        tmpElement.stKepler.dE = vKepler(1);
        tmpElement.stKepler.dI = vKepler(2);
        tmpElement.stKepler.dRAAN = vKepler(3);
        tmpElement.stKepler.dW = vKepler(4);
        tmpElement.stKepler.dMA = vKepler(5);
        vector<Satellite_Element> vOut = CreateConstellatory(tmpElement,5,6,1,D2PI);

        CDate dataBJ(tmpSGP4.GetTLEEpoch(),BJ);

        int nYear,nMonth,nDay,nHour,nMinute;
        double dSeconds;
        dataBJ.GetDate(nYear,nMonth,nDay,nHour,nMinute,dSeconds);
        cout<<dataBJ<<endl;
        for(int nIndex=0,nSize=vOut.size();nIndex<nSize;++nIndex)
        {
            const Satellite_Element& one = vOut[nIndex];
            cout<<"-------------"<<nIndex<<"-------------\ndA:\t"<<one.stKepler.dA<<'\t'
               <<"dE:\t"<<one.stKepler.dE<<'\t'
              <<"dI:\t"<<one.stKepler.dI*DR2D<<'\n'
             <<"dRAAN:\t"<<one.stKepler.dRAAN*DR2D<<'\t'
            <<"dW:\t"<<one.stKepler.dW*DR2D<<'\t'
            <<"dMA:\t"<<one.stKepler.dMA*DR2D<<endl;
        }


        cout<<"=============================\ndA:\t"<<vKepler(0)<<'\t'
           <<"dE:\t"<<vKepler(1)<<'\t'
          <<"dI:\t"<<vKepler(2)*DR2D<<'\n'
         <<"dRAAN:\t"<<vKepler(3)*DR2D<<'\t'
        <<"dW:\t"<<vKepler(4)*DR2D<<'\t'
        <<"dMA:\t"<<vKepler(5)*DR2D<<endl;
        cout<<"DMJD:\t"<<tmpSGP4.GetTLEEpoch()<<endl;

        char tleOut[2][73];
        CDate dataTemp(2020,3,1,0,0,0,UTC);
        vKepler(0)=6814585.9394080322;
        vKepler(1)=0.00029999999999999997;
        vKepler(2)=0.69868322484136203;
        vKepler(3)=5.1470178520333256;
        vKepler(4)=0.88560449707145072;
        vKepler(5)=4.2943616805810239;

        Satellite::CKepler::Classical2TLE(vKepler,dataTemp.GetJD(),90103,tleOut);
        cout<<tleOut[0]<<tleOut[1]<<endl;

        tmpSGP4.SetTLE(tleOut[0],tleOut[1]);
        cout<<tmpSGP4<<endl;

        cout<<setprecision(6)<<tmpSGP4.CalPV(58909.0000000000000000)<<endl;

        vKepler = tmpSGP4.ClassicalElements();
        cout<<"=============================\ndA:\t"<<vKepler(0)<<'\t'
           <<"dE:\t"<<vKepler(1)<<'\t'
          <<"dI:\t"<<vKepler(2)*DR2D<<'\n'
         <<"dRAAN:\t"<<vKepler(3)*DR2D<<'\t'
        <<"dW:\t"<<vKepler(4)*DR2D<<'\t'
        <<"dMA:\t"<<vKepler(5)*DR2D<<endl;
        cout<<"DMJD:\t"<<tmpSGP4.GetTLEEpoch()<<endl;
    }

    /// 测试
    if(1)
    {
        Satellite_Element tmpElement;
        tmpElement.elemType = SAT_TWOBODY;
        tmpElement.stKepler.dA = 42167e3;
        tmpElement.stKepler.dE = 0;
        tmpElement.stKepler.dI = 0;
        tmpElement.stKepler.dRAAN = 327*DD2R;
        tmpElement.stKepler.dW = 0;
        tmpElement.stKepler.dMA = 0;
        tmpElement.stKepler.stEpoch.nYear = 2018;
        tmpElement.stKepler.stEpoch.uMonth = 3;
        tmpElement.stKepler.stEpoch.uDay = 1;
        tmpElement.stKepler.stEpoch.uHour = 8;
        tmpElement.stKepler.stEpoch.uMinute=00;
        tmpElement.stKepler.stEpoch.uSecond=0;
        tmpElement.stKepler.stEpoch.uMSecond=0;

        BJTime endTime;
        endTime.nYear = 2019;
        endTime.uMonth = 3;
        endTime.uDay = 1;
        endTime.uHour = 8;
        endTime.uMinute=0;
        endTime.uSecond=0;
        endTime.uMSecond=0;

//        Pos stPos,stPos1;
//        stPos.dX = 120.257291;
//        stPos.dY = 24.080290;
//        stPos.dZ = 0;
//        stPos1.dX=stPos1.dY=stPos1.dZ=0;


//        vector<Period> periods = VisiblePeriod(beginTime,endTime,tmpElement,stPos,Rota_PRY,stPos1,30*DD2R,30*DD2R,eEllipse);
//        cout<<"periods Size:"<<periods.size()<<endl;

        SatellitePos pos;
        clock_t begine = clock();
//        for(int i=0;i<1000;++i)
        {
            TwoBody(tmpElement.stKepler.stEpoch,endTime,60,tmpElement.stKepler,pos);

        }
        clock_t end = clock();
        std::cout<<double (end-begine)/CLOCKS_PER_SEC<<std::endl;
        std::cout<<pos.vECF.size()<<endl;
        cout<<setprecision(16)<<setw(25)<<pos.vJ2000[0].stP.dX<<','<<pos.vJ2000[0].stP.dY<<','<<pos.vJ2000[0].stP.dZ<<endl;
        cout<<setprecision(16)<<setw(25)<<pos.vJ2000[pos.vJ2000.size()-1].stP.dX<<','
            <<pos.vJ2000[pos.vJ2000.size()-1].stP.dY<<','<<pos.vJ2000[pos.vJ2000.size()-1].stP.dZ<<endl;

        cout<<setprecision(16)<<setw(25)<<pos.vECF[0].stP.dX<<','<<pos.vECF[0].stP.dY<<','<<pos.vECF[0].stP.dZ<<endl;
        cout<<setprecision(16)<<setw(25)<<pos.vECF[pos.vECF.size()-1].stP.dX<<','
            <<pos.vECF[pos.vECF.size()-1].stP.dY<<','<<pos.vECF[pos.vECF.size()-1].stP.dZ<<endl;
//        for(auto one : pos.vJ2000)
//        {
//            cout<<one.stP.dX<<','<<one.stP.dY<<','<<one.stP.dZ<<endl;
//        }
//        for(auto one : pos.vECF)
//        {
//            cout<<one.stP.dX<<','<<one.stP.dY<<','<<one.stP.dZ<<endl;
//        }
//        for(auto one : pos.vLLA)
//        {
//            cout<<one.dX<<','<<one.dY<<','<<one.dZ<<endl;
//        }
        double dX(-4.1283223155120395E7),dY(8587564.931477092),dZ(71438.5903900467);
        double dL,dB,dLat;
        GisMath::XYZ2LBH(dX,dY,dZ,dL,dB,dLat);
        std::cout<<dL*DR2D<<','<<dB*DR2D<<','<<dLat<<std::endl;

        Satellite::CSGP4 time("1 41886U 16078C   21102.11088259  .00000442  00000-0  24494-4 0  9996","2 41886  34.9600 357.9085 0015283 204.3567 155.6423 15.15463304239455");
        double dDouble = time.GetTLEEpoch();
        CDate data(dDouble,BJ);
        cout<<data<<endl;

//        CDate data(2021,3,2,8,0,0);
        CTimeSys timeSys(data);

        std::cout<<setw(25)<<Aerospace::CJPLEphemeris::GetInstance()->GetSunPos(timeSys.GetTT());
        std::cout<<setw(25)<<Aerospace::CJPLEphemeris::GetInstance()->GetMoonPos(timeSys.GetTT());

        CQuaternion tmpQuat(CVector::X_AXIS,45*DD2R);
        CVector v(1,0,0);
        int nTime(0);
        begine = clock();
        CQuaternion tmpQuat2;
        for(int i=0;i<nTime;++i)
        {
            tmpQuat2.Rebuild(CVector::X_AXIS,0);
            tmpQuat *= tmpQuat2;
            tmpQuat2.Rebuild(CVector::X_AXIS,0);
            tmpQuat *= tmpQuat2;
            tmpQuat2.Rebuild(CVector::X_AXIS,0);
            tmpQuat *= tmpQuat2;
            tmpQuat2.Rebuild(CVector::X_AXIS,0);
            tmpQuat *= tmpQuat2;
            tmpQuat2.Rebuild(CVector::X_AXIS,0);
            tmpQuat *= tmpQuat2;
            tmpQuat2.Rebuild(CVector::X_AXIS,0);
            tmpQuat *= tmpQuat2;
            tmpQuat2.Rebuild(CVector::X_AXIS,0);
            tmpQuat *= tmpQuat2;
            tmpQuat2.Rebuild(CVector::X_AXIS,0);
            tmpQuat *= tmpQuat2;

            tmpQuat.Translate(v,v);
        }
        std::cout<<double (clock()-begine) / CLOCKS_PER_SEC<<std::endl;
        std::cout<<setprecision(11)<<setw(20)<<CCoorSys::J20002ECF(dDouble)<<std::endl;

        CVector vJ2000(pos.vJ2000[0].stP.dX,pos.vJ2000[0].stP.dY,pos.vJ2000[0].stP.dZ);
        CDate testData(pos.stStart.nYear,pos.stStart.uMonth,pos.stStart.uDay,
                       pos.stStart.uHour,pos.stStart.uMinute,pos.stStart.uSecond+pos.stStart.uMSecond/1000.,BJ);
        CVector vICRF = CCoorSys::J20002MOD(testData.GetMJD()) * vJ2000;
        std::cout<<setprecision(11)<<setw(25)<<vICRF<<std::endl;
        vICRF = CCoorSys::J20002TOD(testData.GetMJD()) * vJ2000;
        std::cout<<setprecision(11)<<setw(25)<<vICRF<<std::endl;
        vICRF = CCoorSys::J20002TEME(testData.GetMJD()) * vJ2000;
        std::cout<<setprecision(11)<<setw(25)<<vICRF<<std::endl;

////        CDate testData(2009,1,1,4,0,0,UTC);

        CVector vTest(6678140,0,28.5*DD2R,0,0,0);
        Satellite::CJ2 j2(6678140,0,28.5*DD2R,0,0,0);
//        for(int i=0;i<86400;++i)
        {
            std::cout<<setprecision(11)<<setw(25)<<j2.CalPV(0)<<endl;
            std::cout<<setprecision(11)<<setw(25)<<j2.CalPV(86400)<<endl;
        }
//        vJ2000.Set(-4953.733255,-3936.230829,-2136.381869);
//        vICRF = CCoorSys::TEME2ECF(testData.GetMJD()) * vJ2000;
//        std::cout<<setprecision(11)<<setw(20)<<vICRF<<std::endl;
    }

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
//    GisMath::InitGis(GisMath::WGS_84);

//    Math::CVector vEye(3);
//    Math::CVector vDir(3);
//    Math::CVector vOut(3);

//    vDir(0) = tan(1*DD2R);
//    vDir(1) = vDir(0);
//    vDir(2) = 1;
//    if(GisMath::CalLineInterEarth(vEye,vDir,vOut))
//    {
//        cout<<setw(20)<<vOut<<endl;
//    }
//    double dEarth_a(6378137);
//    double dEarth_b(6356752.314245179);
//    vDir(0) = 1.;
//    vDir(1) = 0.;
//    vDir(2) = 0.;
//    cout<<dEarth_a<<endl;
//    if(GisMath::CalLineInterEarth(vEye,vDir,vOut))
//    {
//        cout<<setw(20)<<vOut<<endl;
//    }



    /// 释放资源
    CloseSatelliteToolKit();
    return 0;
}
