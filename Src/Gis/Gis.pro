#-------------------------------------------------
#
# Project created by QtCreator 2015-11-22T13:21:00
#
#-------------------------------------------------

QT       -= gui

TEMPLATE = lib

DEFINES += GIS_LIBRARY
STK_PATH=$${PWD}/../..
STK_CONFIG += Math
INCLUDEPATH *= $${STK_PATH}/Public/Include/Gis
INCLUDEPATH *= $${STK_PATH}/Public/Include/Sofa
include($$STK_PATH/stkconfig.pri)
SOURCES += Src/GisMath.cpp

CONFIG (debug, debug|release){
    TARGET = GisD
}
else{
    TARGET = Gis
}

### Linux 或 Mac 环境
unix{
    DESTDIR = $${PWD}/../../Bin
    VERSION = 2.0.0
    QMAKE_LFLAGS += -Wl,-rpath,.
}
