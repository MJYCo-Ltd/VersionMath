#-------------------------------------------------
#
# Project created by QtCreator 2015-11-22T13:21:00
#
#-------------------------------------------------

CONFIG -= qt

TEMPLATE = lib

DEFINES += GIS_LIBRARY

SDK_CONFIG *= MATH
include($$PWD/../VersionMath.pri)

INCLUDEPATH *= Inc $$PWD/../Satellite/Inc/Sofa

HEADERS += \
    Inc/GisMath_Common.h

SOURCES += \
    Src/GisMath_Common.cpp \
    Src/GisMath_Geographic.cpp \
    Src/GisMath_Transform.cpp \
    Src/geodesic.c \
    ../Satellite/Src/Sofa/anp.c \
    ../Satellite/Src/Sofa/eform.c \
    ../Satellite/Src/Sofa/gc2gd.c \
    ../Satellite/Src/Sofa/gc2gde.c \
    ../Satellite/Src/Sofa/gd2gc.c \
    ../Satellite/Src/Sofa/gd2gce.c \
    ../Satellite/Src/Sofa/zp.c
