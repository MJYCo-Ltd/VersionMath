#-------------------------------------------------
#
# Project created by QtCreator 2015-11-22T13:21:00
#
#-------------------------------------------------

CONFIG -= qt

TEMPLATE = lib

DEFINES += GIS_LIBRARY

SDK_PATH=$$PWD/../..
SDK_CONFIG *= SOFA MATH
include($$PWD/../VersionMath.pri)

INCLUDEPATH *= $$SDK_PATH/Inc/GisMath Inc

HEADERS += \
    Inc/GisMath_Common.h

SOURCES += \
    Src/GisMath_Common.cpp \
    Src/GisMath_Geographic.cpp \
    Src/GisMath_Transform.cpp

LIBS *= -L$$PWD/Lib
win32-msvc2015{
    CONFIG(debug, debug|release) {
      LIBS *= -lproj_5_0_d
    }else{
      LIBS *= -lproj_5_0
    }
}else{
    CONFIG(debug, debug|release) {
      LIBS *= -lprojd
    }else{
      LIBS *= -lproj
    }
}
