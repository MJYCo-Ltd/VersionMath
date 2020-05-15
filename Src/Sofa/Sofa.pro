#-------------------------------------------------
#
# Project created by QtCreator 2015-02-27T14:26:31
#
#-------------------------------------------------

CONFIG -= qt

TEMPLATE = lib

# 预定义宏
DEFINES += SOFA_LIBRARY

SDK_PATH = $$PWD/../..

INCLUDEPATH *= $$SDK_PATH/Inc/Sofa
win32-msvc:QMAKE_CXXFLAGS += -utf-8
win32-msvc:QMAKE_CXXFLAGS += /wd"4100"
win32-msvc:RC_FILE = Src/$${TARGET}_Version.rc
win32-msvc:DLLDESTDIR = $${SDK_PATH}/../Bin
win32-msvc:DESTDIR = $${SDK_PATH}/Lib

### Linux 或 Mac 环境
unix{
    DESTDIR = $${SDK_PATH}/Bin
    VERSION = 2.0.0
    QMAKE_LFLAGS += -Wl,-rpath,.
}

CONFIG(debug, debug|release){
  TARGET = $$join(TARGET,,,D)
}
# 目标文件的输出路径  end

# 源代码
SOURCES += \
    Src/zr.cpp \
    Src/zpv.cpp \
    Src/zp.cpp \
    Src/xys06a.cpp \
    Src/xys00b.cpp \
    Src/xys00a.cpp \
    Src/xy06.cpp \
    Src/utcut1.cpp \
    Src/utctai.cpp \
    Src/ut1utc.cpp \
    Src/ut1tt.cpp \
    Src/ut1tai.cpp \
    Src/ttut1.cpp \
    Src/tttdb.cpp \
    Src/tttcg.cpp \
    Src/tttai.cpp \
    Src/trxpv.cpp \
    Src/trxp.cpp \
    Src/tr.cpp \
    Src/tf2d.cpp \
    Src/tf2a.cpp \
    Src/tdbtt.cpp \
    Src/tdbtcb.cpp \
    Src/tcgtt.cpp \
    Src/tcbtdb.cpp \
    Src/taiutc.cpp \
    Src/taiut1.cpp \
    Src/taitt.cpp \
    Src/sxpv.cpp \
    Src/sxp.cpp \
    Src/starpv.cpp \
    Src/starpm.cpp \
    Src/sp00.cpp \
    Src/seps.cpp \
    Src/sepp.cpp \
    Src/s06a.cpp \
    Src/s06.cpp \
    Src/s2xpv.cpp \
    Src/s2pv.cpp \
    Src/s2p.cpp \
    Src/s2c.cpp \
    Src/s00b.cpp \
    Src/s00a.cpp \
    Src/s00.cpp \
    Src/rz.cpp \
    Src/ry.cpp \
    Src/rxr.cpp \
    Src/rxpv.cpp \
    Src/rxp.cpp \
    Src/rx.cpp \
    Src/rv2m.cpp \
    Src/rm2v.cpp \
    Src/refco.cpp \
    Src/pxp.cpp \
    Src/pvxpv.cpp \
    Src/pvup.cpp \
    Src/pvu.cpp \
    Src/pvtob.cpp \
    Src/pvstar.cpp \
    Src/pvppv.cpp \
    Src/pvmpv.cpp \
    Src/pvm.cpp \
    Src/pvdpv.cpp \
    Src/pv2s.cpp \
    Src/pv2p.cpp \
    Src/prec76.cpp \
    Src/pr00.cpp \
    Src/ppsp.cpp \
    Src/ppp.cpp \
    Src/pom00.cpp \
    Src/pnm80.cpp \
    Src/pnm06a.cpp \
    Src/pnm00b.cpp \
    Src/pnm00a.cpp \
    Src/pn06a.cpp \
    Src/pn06.cpp \
    Src/pn00b.cpp \
    Src/pn00a.cpp \
    Src/pn00.cpp \
    Src/pn.cpp \
    Src/pmsafe.cpp \
    Src/pmpx.cpp \
    Src/pmp.cpp \
    Src/pmat76.cpp \
    Src/pmat06.cpp \
    Src/pmat00.cpp \
    Src/pm.cpp \
    Src/plan94.cpp \
    Src/pfw06.cpp \
    Src/pdp.cpp \
    Src/pb06.cpp \
    Src/pas.cpp \
    Src/pap.cpp \
    Src/p06e.cpp \
    Src/p2s.cpp \
    Src/p2pv.cpp \
    Src/obl80.cpp \
    Src/obl06.cpp \
    Src/nutm80.cpp \
    Src/nut80.cpp \
    Src/nut06a.cpp \
    Src/nut00b.cpp \
    Src/nut00a.cpp \
    Src/numat.cpp \
    Src/num06a.cpp \
    Src/num00b.cpp \
    Src/num00a.cpp \
    Src/ldsun.cpp \
    Src/ldn.cpp \
    Src/ld.cpp \
    Src/jdcalf.cpp \
    Src/jd2cal.cpp \
    Src/ir.cpp \
    Src/hfk5z.cpp \
    Src/h2fk5.cpp \
    Src/gst94.cpp \
    Src/gst06a.cpp \
    Src/gst06.cpp \
    Src/gst00b.cpp \
    Src/gst00a.cpp \
    Src/gmst82.cpp \
    Src/gmst06.cpp \
    Src/gmst00.cpp \
    Src/gd2gce.cpp \
    Src/gd2gc.cpp \
    Src/gc2gde.cpp \
    Src/gc2gd.cpp \
    Src/fw2xy.cpp \
    Src/fw2m.cpp \
    Src/fk52h.cpp \
    Src/fk5hz.cpp \
    Src/fk5hip.cpp \
    Src/fave03.cpp \
    Src/faur03.cpp \
    Src/fasa03.cpp \
    Src/fapa03.cpp \
    Src/faom03.cpp \
    Src/fane03.cpp \
    Src/fame03.cpp \
    Src/fama03.cpp \
    Src/falp03.cpp \
    Src/fal03.cpp \
    Src/faju03.cpp \
    Src/faf03.cpp \
    Src/fae03.cpp \
    Src/fad03.cpp \
    Src/era00.cpp \
    Src/eqeq94.cpp \
    Src/epv00.cpp \
    Src/epj2jd.cpp \
    Src/epj.cpp \
    Src/epb2jd.cpp \
    Src/epb.cpp \
    Src/eors.cpp \
    Src/eo06a.cpp \
    Src/eform.cpp \
    Src/eect00.cpp \
    Src/ee06a.cpp \
    Src/ee00b.cpp \
    Src/ee00a.cpp \
    Src/ee00.cpp \
    Src/dtf2d.cpp \
    Src/dtdb.cpp \
    Src/dat.cpp \
    Src/d2tf.cpp \
    Src/d2dtf.cpp \
    Src/cr.cpp \
    Src/cpv.cpp \
    Src/cp.cpp \
    Src/cal2jd.cpp \
    Src/c2txy.cpp \
    Src/c2tpe.cpp \
    Src/c2teqx.cpp \
    Src/c2tcio.cpp \
    Src/c2t06a.cpp \
    Src/c2t00b.cpp \
    Src/c2t00a.cpp \
    Src/c2s.cpp \
    Src/c2ixys.cpp \
    Src/c2ixy.cpp \
    Src/c2ibpn.cpp \
    Src/c2i06a.cpp \
    Src/c2i00b.cpp \
    Src/c2i00a.cpp \
    Src/bpn2xy.cpp \
    Src/bp06.cpp \
    Src/bp00.cpp \
    Src/bi00.cpp \
    Src/atoiq.cpp \
    Src/atoi13.cpp \
    Src/atoc13.cpp \
    Src/atioq.cpp \
    Src/atio13.cpp \
    Src/aticqn.cpp \
    Src/aticq.cpp \
    Src/atic13.cpp \
    Src/atco13.cpp \
    Src/atciqz.cpp \
    Src/atciqn.cpp \
    Src/atciq.cpp \
    Src/atci13.cpp \
    Src/apio13.cpp \
    Src/apio.cpp \
    Src/aper13.cpp \
    Src/aper.cpp \
    Src/apcs13.cpp \
    Src/apcs.cpp \
    Src/apco13.cpp \
    Src/apco.cpp \
    Src/apci13.cpp \
    Src/apci.cpp \
    Src/apcg13.cpp \
    Src/apcg.cpp \
    Src/anpm.cpp \
    Src/anp.cpp \
    Src/af2a.cpp \
    Src/ab.cpp \
    Src/a2tf.cpp \
    Src/a2af.cpp \
    Src/g2icrs.cpp \
    Src/icrs2g.cpp \
    Src/eceq06.cpp \
    Src/ecm06.cpp \
    Src/eqec06.cpp \
    Src/lteceq.cpp \
    Src/ltecm.cpp \
    Src/lteqec.cpp \
    Src/ltp.cpp \
    Src/ltpb.cpp \
    Src/ltpecl.cpp \
    Src/ltpequ.cpp
