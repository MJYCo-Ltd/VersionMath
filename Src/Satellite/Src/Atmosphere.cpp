#include <Math/Vector.h>
#include <Math/MathCommonAlgorithm.h>
#include <GisMath/GisMath.h>
#include <Satellite/Atmosphere.h>
#include <Satellite/CoorSys.h>
#include "CommonAlgorithm.h"
#include "NrlMsise_00/nrlmsise-00.h"
#include "sofam.h"

using namespace Aerospace;
using namespace Math;
using namespace Physical;
CAtmosphere::CAtmosphere(AtmosModel type):m_typeAtmos(type),m_bInit(true)
{
    SetModelName();
}

CAtmosphere::~CAtmosphere()
{
}

void CAtmosphere::SetModelName()
{
    switch(m_typeAtmos)
    {
    case HarrisPriester:
        m_strModelName = "Harris-Priester";
        break;
    case NRLMSISE2000:
        m_strModelName = "NRLMSISE 2000";
        break;
    case STANDARD1976:
        m_strModelName = "1976 Standard";
        break;
    default:
        m_strModelName = "UnKnow";
        break;
    }
}

double CAtmosphere::GetDensity(const double &dMJD, const CVector &rSat, const CVector &rSun)
{
    return(GetDensity(dMJD,rSat,rSun,m_typeAtmos));
}

/// 获取大气的密度信息
double CAtmosphere::GetDensity(const double &dMJD, const CVector &rSat, const CVector &rSun,
                                   AtmosModel typeAtmos)
{
    switch(typeAtmos)
    {
    case HarrisPriester:
        /// STK 相应的密度模型有效范围为 0 ~ 1000 km
        /// 而本模型 有效范围为 100~1000km
        return(Density_HP(rSun,rSat));
    case NRLMSISE2000:
    {
        static nrlmsise_flags flag={{1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1}};
        static nrlmsise_input input;
        static nrlmsise_output output;

        double dL,dB,dH;
        GisMath::XYZ2LBH(rSat(0),rSat(1),rSat(2),dL,dB,dH);

        ///  STK对应的大气密度模型的有效范围为 0~1000km
        ///  而本模型有效范围为 >0km，且当低于 0km时，按照0km计算
        ///  为了与STK比对，此处进行处理，当dH高于1000km时，直接返回

        /// 因nrlmsise_input需要的参数 分别 为 [km] [deg] [deg]
        /// 因此需要进行转换
        input.alt = dH*1.e-3;
        input.g_lat = dB*DR2D;
        input.g_long = dL*DR2D;
        input.f107A = 150.;
        input.f107  = 150.;
        input.ap    = 4.;

        /// 以2000年1月1日0时0分0秒为基准根据年份减少天数
        double dTempLoal = dMJD - DJM00 + 0.5;

        /// DJY 是365.25 认为每年都是这么多天，而我们平时采用的是闰年操作
        /// 因此离闰年越远，需要减去的天数越大。
        int nJudge = int(dTempLoal/DJY) % 4;

        if(dTempLoal < 0)
        {
            switch(nJudge)
            {
            case 0:
                dTempLoal -= 0.25;
                break;
            case -1:
                dTempLoal -= 0.5;
                break;
            case -2:
                dTempLoal -= 0.75;
                break;
            }
        }
        else
        {
            switch(nJudge)
            {
            case 1:
                dTempLoal -= 0.75;
                break;
            case 2:
                dTempLoal -= 0.5;
                break;
            case 3:
                dTempLoal -= 0.25;
                break;
            }
        }

        /// 设置当前日期是这一年的第几天
        input.doy = Modulo(dTempLoal,DJY);

        /// 因约简儒略日开始时间是正午，因此要将这半天加上
        /// 设置当前时间，在这一天的秒数
        input.sec = Frac(dMJD-DJM00+0.5)*DAYSEC;

        /// 将UTC时间转换成当地时间
        input.lst = input.sec/3600.0 + input.g_long/15.0;

        gtd7d(&input, &flag, &output);
        return(output.d[5]*1.e-3);
    }
        break;
    case STANDARD1976:
        /// STK 相应的密度模型有效范围为 86 ~ 1000 km
        /// 而本模型的有效范围为 0~1000km，且当物体低于0时按照0km计算
        /// 将ECF转换成地理坐标
        double dL,dB,dH;
        GisMath::XYZ2LBH(rSat(0),rSat(1),rSat(2),dL,dB,dH);
        return(SA76(dH));
    default:
        return(0.);
    }
    return(0.);
}
