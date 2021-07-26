#include <Satellite/IRESInfo.h>
#include <Satellite/JPLEphemeris.h>
#include <SatelliteToolKit/SatelliteToolKit.h>

using namespace Aerospace;

bool InitSatelliteToolKit(const std::string& sPath,std::string& sErrorInfo)
{
    /// 初始化地球极移数据
    CIRESInfo* pIERS = CIRESInfo::GetInstance();

    if(nullptr == pIERS)
    {
        sErrorInfo = "Can't Create IERS!";
        return(false);
    }

    if(!pIERS->IsInit())
    {
        string filePath = sPath + "/Data/STKData/finals2000A.data";
        if(!pIERS->Init(filePath))
        {
            sErrorInfo = filePath + " does not exist or is malformed!";
            return false;
        }
    }

    /// 初始化JPL星历数据
    CJPLEphemeris* pJPL = CJPLEphemeris::GetInstance();
    if(nullptr == pJPL)
    {
        sErrorInfo = "Can't Create JPLEphemeris!";
        return(false);
    }

    if(!pJPL->IsInit())
    {
        string filePath = sPath + "/Data/STKData/ephem";
        if(!pJPL->Init(filePath))
        {
            sErrorInfo = filePath + " does not exist or is malformed!";
            return false;
        }
    }

    return(true);
}

void CloseSatelliteToolKit()
{
    /// 初始化JPL星历数据
    CJPLEphemeris* pJPL = CJPLEphemeris::GetInstance();
    if(pJPL->IsInit())
    {
        CJPLEphemeris::Realse();
    }
}
