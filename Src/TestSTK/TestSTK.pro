TEMPLATE = app
CONFIG += console c++11
CONFIG -= app_bundle
CONFIG -= qt

SDK_PATH=$$PWD/../..
SDK_CONFIG *= MATH GIS_MATH SATELLITE_TOOL_KIT SATELLITE
include($$PWD/../VersionMath.pri)
# 配置依赖的算法、及头文件

SOURCES += \
        Src/main.cpp
