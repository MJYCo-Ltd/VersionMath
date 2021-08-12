# 开启utf-8 编码方式支持
win32{
    QMAKE_CXXFLAGS += -utf-8
    QMAKE_CXXFLAGS += /wd"4100"
    RC_FILE = Src/$${TARGET}_Version.rc
    DLLDESTDIR = $${SDK_PATH}/../Bin
    DESTDIR = $${SDK_PATH}/Lib
}



### Linux 或 Mac 环境
unix{
    DESTDIR = $${SDK_PATH}/../Bin/stklib
    VERSION = 2.0.0
    QMAKE_LFLAGS += -Wl,-rpath,.
    QMAKE_CXXFLAGS += -fvisibility=hidden
}

INCLUDEPATH *= $${SDK_PATH}/Inc/
LIBS *= -L$$DESTDIR

contains(SDK_CONFIG,MATH){
    CONFIG(debug, debug|release) {
      LIBS *= -lMathd
    }else{
      LIBS *= -lMath
    }
}

contains(SDK_CONFIG,SATELLITE){
    CONFIG(debug, debug|release) {
      LIBS *= -lSatellited
    }else{
      LIBS *= -lSatellite
    }
}

contains(SDK_CONFIG,GIS_MATH){
    CONFIG(debug, debug|release) {
      LIBS *= -lGisMathd
    }else{
      LIBS *= -lGisMath
    }
}

contains(SDK_CONFIG,SATELLITE_TOOL_KIT){
    CONFIG(debug, debug|release) {
      LIBS *= -lSatelliteToolKitd
    }else{
      LIBS *= -lSatelliteToolKit
    }
}

CONFIG(debug, debug|release){
  TARGET = $$join(TARGET,,,d)
}


contains(TEMPLATE, "app"){
    win32:DESTDIR = $${SDK_PATH}/../Bin
}

