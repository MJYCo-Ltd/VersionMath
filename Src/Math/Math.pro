#-------------------------------------------------
#
# Project created by QtCreator 2015-11-22T13:21:00
#
#-------------------------------------------------

CONFIG -= qt

TEMPLATE = lib

DEFINES *= MATH_LIBRARY

SDK_PATH=$$PWD/../..
SDK_CONFIG *= SOFA
include($$PWD/../VersionMath.pri)

INCLUDEPATH *= $$SDK_PATH/Inc/Math

# 目标文件的输出路  end

SOURCES += \
    Src/DE.cpp \
    Src/EulerAngle.cpp \
    Src/GJ4.cpp \
    Src/Intpol.cpp \
    Src/Matrix.cpp \
    Src/Quaternion.cpp \
    Src/RK4.cpp \
    Src/RKF78.cpp \
    Src/Transit.cpp \
    Src/VecMat.cpp \
    Src/Vector.cpp \
    Src/YPRAngle.cpp \
    Src/MathCommonAlgorithm.cpp
