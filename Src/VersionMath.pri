win32{
    INCLUDEPATH *= $$PWD/../Inc
    LIBS *= -L$$PWD/../Lib
    QMAKE_CXXFLAGS += -utf-8
    QMAKE_CXXFLAGS += /wd"4100"
    RC_FILE = Src/$${TARGET}_Version.rc

    contains(TEMPLATE,"lib"){
        DLLDESTDIR = $$PWD/../../Bin
        DESTDIR = $$PWD/../Lib
    }else{
        DESTDIR = $$PWD/../../Bin
    }
}

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

CONFIG(debug,debug|release){
    TARGET = $$join(TARGET,,,d)
}
