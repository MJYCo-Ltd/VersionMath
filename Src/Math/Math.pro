#-------------------------------------------------
#
# Project created by QtCreator 2015-11-22T13:21:00
#
#-------------------------------------------------

QT       -= gui

TEMPLATE = lib

DEFINES *= MATH_LIBRARY

INCLUDEPATH *= ../../Public/Include/Sofa/
INCLUDEPATH *= ../../Public/Include/Math/

# 目标文件的输出路 begin
### win32 环境
win32{
    contains(TEMPLATE, "lib") {
        win32-msvc2010{
            RC_FILE = Src/Math_Version.rc
            DESTDIR = $${PWD}/../../Lib
            DLLDESTDIR = $${PWD}/../../Bin/VS2010
        }
        win32-msvc{
            RC_FILE = Src/Math_Version.rc
            DESTDIR = $${PWD}/../../Lib
            DLLDESTDIR = $${PWD}/../../Bin
        }
        win32-msvc2015{
            RC_FILE = Src/Math_Version.rc
            DESTDIR = $${PWD}/../../Lib/VS2015
            DLLDESTDIR = $${PWD}/../../Bin/VS2015
        }
        win32-g++{
            DESTDIR = $${PWD}/../../Bin/g++
        }
    } else {
      DESTDIR = $${PWD}/../../Bin
    }
}

### Linux  Mac 环境
unix{
    DESTDIR = $${PWD}/../../Bin
    VERSION = 2.0.0
    QMAKE_LFLAGS += -Wl,-rpath,.
}

CONFIG(debug, debug|release) {
  TARGET = $$join(TARGET,,,D)
}else{
}

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
    Src/VecMat.cpp \
    Src/Vector.cpp \
    Src/YPRAngle.cpp \
    Src/MathCommonAlgorithm.cpp
