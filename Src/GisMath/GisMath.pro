#-------------------------------------------------
#
# Project created by QtCreator 2015-11-22T13:21:00
#
#-------------------------------------------------

CONFIG -= qt

TEMPLATE = lib

DEFINES += GIS_LIBRARY

SDK_PATH=$$PWD/../..
SDK_CONFIG *= MATH
include($$PWD/../VersionMath.pri)

INCLUDEPATH *= $$SDK_PATH/Inc/GisMath Inc

HEADERS += \
    Inc/GisMath_Common.h

SOURCES += \
    ../Satellite/Src/Sofa/anp.cpp \
    ../Satellite/Src/Sofa/eform.cpp \
    ../Satellite/Src/Sofa/gc2gd.cpp \
    ../Satellite/Src/Sofa/gc2gde.cpp \
    ../Satellite/Src/Sofa/gd2gc.cpp \
    ../Satellite/Src/Sofa/gd2gce.cpp \
    ../Satellite/Src/Sofa/zp.cpp \
    Src/GisMath_Common.cpp \
    Src/GisMath_Geographic.cpp \
    Src/GisMath_Transform.cpp \
    Src/geodesic.c
