CONFIG -= qt

TEMPLATE = lib
DEFINES += SATELLITETOOLKIT_LIBRARY

CONFIG += c++11

# The following define makes your compiler emit warnings if you use
# any Qt feature that has been marked deprecated (the exact warnings
# depend on your compiler). Please consult the documentation of the
# deprecated API in order to know how to port your code away from it.
DEFINES += QT_DEPRECATED_WARNINGS

# You can also make your code fail to compile if it uses deprecated APIs.
# In order to do so, uncomment the following line.
# You can also select to disable deprecated APIs only up to a certain version of Qt.
#DEFINES += QT_DISABLE_DEPRECATED_BEFORE=0x060000    # disables all the APIs deprecated before Qt 6.0.0

SDK_PATH=$$PWD/../..
SDK_CONFIG *= MATH GIS_MATH SATELLITE
include($$PWD/../VersionMath.pri)

INCLUDEPATH *= $$SDK_PATH/Inc/SatelliteToolKit

SOURCES += \
    ../Satellite/Src/Sofa/anp.cpp \
    Src/SatelliteToolKitCommon.cpp \
    Src/SatelliteToolKit_Angle.cpp \
    Src/SatelliteToolKit_Init.cpp \
    Src/SatelliteToolKit_Intersect.cpp \
    Src/SatelliteToolKit_Oribit.cpp \
    Src/SatelliteToolKit_Period.cpp \
    Src/SatelliteToolKit_Positioning.cpp \
    Src/SatelliteToolKit_Rotate.cpp \
    Src/SatelliteToolKit_Visible.cpp


HEADERS += \
    Inc/SatelliteToolKitCommon.h
