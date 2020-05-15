# 开启utf-8 编码方式支持
win32-msvc:QMAKE_CXXFLAGS += -utf-8
win32-msvc:QMAKE_CXXFLAGS += /wd"4100"
win32-msvc:RC_FILE = Src/$${TARGET}_Version.rc
win32-msvc:DLLDESTDIR = $${SDK_PATH}/../Bin
win32-msvc:DESTDIR = $${SDK_PATH}/Lib

### Linux 或 Mac 环境
unix{
    DESTDIR = $${SDK_PATH}/Bin
    VERSION = 2.0.0
    QMAKE_LFLAGS += -Wl,-rpath,.
}

LIBS *= -L$$DESTDIR

contains(SDK_CONFIG,MATH){
    INCLUDEPATH *= $${SDK_PATH}/Inc/Math
    CONFIG(debug, debug|release) {
      LIBS *= -lMathD
    }else{
      LIBS *= -lMath
    }
}

contains(SDK_CONFIG,SOFA){
    INCLUDEPATH *= $${SDK_PATH}/Inc/Sofa
    CONFIG(debug, debug|release) {
      LIBS *= -lSofaD
    }else{
      LIBS *= -lSofa
    }
}

contains(SDK_CONFIG,SATELLITE){
    INCLUDEPATH *= $${SDK_PATH}/Inc/Satellite
    CONFIG(debug, debug|release) {
      LIBS *= -lSatelliteD
    }else{
      LIBS *= -lSatellite
    }
}

contains(SDK_CONFIG,GISMATH){
    INCLUDEPATH *= $${SDK_PATH}/Inc/GisMath
    CONFIG(debug, debug|release) {
      LIBS *= -lGisMathD
    }else{
      LIBS *= -lGisMath
    }
}

CONFIG(debug, debug|release){
  TARGET = $$join(TARGET,,,D)
}


