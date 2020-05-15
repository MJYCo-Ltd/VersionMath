
TEMPLATE = lib
CONFIG -= lib_bundle
INCLUDEPATH *= Inc

# Sofa
INCLUDEPATH *= ../../Public/Include/Sofa/

# Math
INCLUDEPATH *= ../../Public/Include/Math/

# SAT
INCLUDEPATH *= ../../Public/Include/Satellite/

#
DEFINES *= ALGORITHM_LIBRARY


# 目标文件的输出路径 begin
### win32 环境
win32{
    contains(TEMPLATE, "lib") {
        win32-msvc2010{
            RC_FILE = Src/Satellite_Version.rc
            DESTDIR = $${PWD}/../../Lib
            DLLDESTDIR = $${PWD}/../../Bin/VS2010
        }
        win32-msvc{
            RC_FILE = Src/Satellite_Version.rc
            DESTDIR = $${PWD}/../../Lib
            DLLDESTDIR = $${PWD}/../../Bin
        }
        win32-msvc2015{
            RC_FILE = Src/Satellite_Version.rc
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

### Linux 或 Mac 环境
unix{
    DESTDIR = $${PWD}/../../Bin
    VERSION = 2.0.0
    QMAKE_LFLAGS += -Wl,-rpath,.
}

CONFIG(debug, debug|release) {
  TARGET = $$join(TARGET,,,D)
  LIBS *= -L$${DESTDIR} -lSofaD -lMathD
}else{
  LIBS *= -L$${DESTDIR} -lSofa -lMath
}

# 目标文件的输出路径  end

SOURCES += \
    Src/jpl_eph/jpleph.cpp \
    Src/SGP4/sgp4ext.cpp \
    Src/SGP4/sgp4io.cpp \
    Src/SGP4/sgp4unit.cpp \
    Src/NrlMsise_00/nrlmsise-00.c \
    Src/NrlMsise_00/nrlmsise-00_data.c \
    Src/CommonAlgorithm.cpp \
    Src/Atmosphere.cpp \
    Src/CoorSys.cpp \
    Src/Date.cpp \
    Src/Force.cpp \
    Src/IRESInfo.cpp \
    Src/JPLEphemeris.cpp \
    Src/Kepler.cpp \
    Src/OpticalAlg.cpp \
    Src/SGP4.cpp \
    Src/STKGraveModel.cpp \
    Src/TimeSys.cpp \
    Src/Transit.cpp


HEADERS += \
    Inc/SGP4/sgp4ext.h \
    Inc/SGP4/sgp4io.h \
    Inc/SGP4/sgp4unit.h \
    Inc/jpl_eph/jpl_int.h \
    Inc/jpl_eph/jpleph.h \
    Inc/NrlMsise_00/nrlmsise-00.h \
    Inc/CommonAlgorithm.h
