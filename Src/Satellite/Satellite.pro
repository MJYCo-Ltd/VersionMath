TEMPLATE = lib
CONFIG -= qt
INCLUDEPATH *= Inc

DEFINES *= ALGORITHM_LIBRARY

SDK_PATH=$$PWD/../..
SDK_CONFIG *= SOFA MATH GIS_MATH
include($$PWD/../VersionMath.pri)
INCLUDEPATH *= $$SDK_PATH/Inc/Satellite

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
    Src/TimeSys.cpp


HEADERS += \
    Inc/SGP4/sgp4ext.h \
    Inc/SGP4/sgp4io.h \
    Inc/SGP4/sgp4unit.h \
    Inc/jpl_eph/jpl_int.h \
    Inc/jpl_eph/jpleph.h \
    Inc/NrlMsise_00/nrlmsise-00.h \
    Inc/CommonAlgorithm.h
