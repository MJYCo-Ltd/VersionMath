#include "IRESInfo.h"
#include "JPLEphemeris.h"
#include "SatelliteToolKit.h"

using namespace Aerospace;

bool InitSatelliteToolKit(std::string& sErrorInfo)
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
        char filePath[] = "STKData/finals2000A.data";
        if(!pIERS->Init(filePath))
        {
            sErrorInfo = "STKData/finals2000A.data does not exist or is malformed!";
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
        char filePath[] = "STKData/ephem";
        if(!pJPL->Init(filePath))
        {
            sErrorInfo = "STKData/ephem does not exist or is malformed!";
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
