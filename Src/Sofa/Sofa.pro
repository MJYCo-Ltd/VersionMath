#-------------------------------------------------
#
# Project created by QtCreator 2015-02-27T14:26:31
#
#-------------------------------------------------

QT       -= gui

TEMPLATE = lib

# 头文件路径
INCLUDEPATH += ../../Public/Include/Sofa/

# 预定义宏
DEFINES += SOFA_LIBRARY

# 目标文件的输出路径 begin
### win32 环境

win32{
    contains(TEMPLATE, "lib"){
        win32-msvc2010{
            RC_FILE = Source/Sofa_Version.rc
            DESTDIR = $${PWD}/../../Lib
            DLLDESTDIR = $${PWD}/../../Bin/VS2010
        }
        win32-msvc{
            RC_FILE = Source/Sofa_Version.rc
            DESTDIR = $${PWD}/../../Lib
            DLLDESTDIR = $${PWD}/../../Bin
        }
        win32-msvc2015{
            RC_FILE = Source/Sofa_Version.rc
            DESTDIR = $${PWD}/../../Lib/VS2015
            DLLDESTDIR = $${PWD}/../../Bin/VS2015
        }
        win32-g++{
            DESTDIR = $${PWD}/../../Bin/g++
        }
    }else{
      DESTDIR = $${PWD}/../../Bin
    }
}

### Linux 或 Mac 环境

macx{
    DESTDIR = $${PWD}/../../Bin
    VERSION = 2.0.0
}

CONFIG(debug, debug|release){
  TARGET = $$join(TARGET,,,D)
}

# 目标文件的输出路径  end

# 源代码
SOURCES += \
    Source/zr.cpp \
    Source/zpv.cpp \
    Source/zp.cpp \
    Source/xys06a.cpp \
    Source/xys00b.cpp \
    Source/xys00a.cpp \
    Source/xy06.cpp \
    Source/utcut1.cpp \
    Source/utctai.cpp \
    Source/ut1utc.cpp \
    Source/ut1tt.cpp \
    Source/ut1tai.cpp \
    Source/ttut1.cpp \
    Source/tttdb.cpp \
    Source/tttcg.cpp \
    Source/tttai.cpp \
    Source/trxpv.cpp \
    Source/trxp.cpp \
    Source/tr.cpp \
    Source/tf2d.cpp \
    Source/tf2a.cpp \
    Source/tdbtt.cpp \
    Source/tdbtcb.cpp \
    Source/tcgtt.cpp \
    Source/tcbtdb.cpp \
    Source/taiutc.cpp \
    Source/taiut1.cpp \
    Source/taitt.cpp \
    Source/sxpv.cpp \
    Source/sxp.cpp \
    Source/starpv.cpp \
    Source/starpm.cpp \
    Source/sp00.cpp \
    Source/seps.cpp \
    Source/sepp.cpp \
    Source/s06a.cpp \
    Source/s06.cpp \
    Source/s2xpv.cpp \
    Source/s2pv.cpp \
    Source/s2p.cpp \
    Source/s2c.cpp \
    Source/s00b.cpp \
    Source/s00a.cpp \
    Source/s00.cpp \
    Source/rz.cpp \
    Source/ry.cpp \
    Source/rxr.cpp \
    Source/rxpv.cpp \
    Source/rxp.cpp \
    Source/rx.cpp \
    Source/rv2m.cpp \
    Source/rm2v.cpp \
    Source/refco.cpp \
    Source/pxp.cpp \
    Source/pvxpv.cpp \
    Source/pvup.cpp \
    Source/pvu.cpp \
    Source/pvtob.cpp \
    Source/pvstar.cpp \
    Source/pvppv.cpp \
    Source/pvmpv.cpp \
    Source/pvm.cpp \
    Source/pvdpv.cpp \
    Source/pv2s.cpp \
    Source/pv2p.cpp \
    Source/prec76.cpp \
    Source/pr00.cpp \
    Source/ppsp.cpp \
    Source/ppp.cpp \
    Source/pom00.cpp \
    Source/pnm80.cpp \
    Source/pnm06a.cpp \
    Source/pnm00b.cpp \
    Source/pnm00a.cpp \
    Source/pn06a.cpp \
    Source/pn06.cpp \
    Source/pn00b.cpp \
    Source/pn00a.cpp \
    Source/pn00.cpp \
    Source/pn.cpp \
    Source/pmsafe.cpp \
    Source/pmpx.cpp \
    Source/pmp.cpp \
    Source/pmat76.cpp \
    Source/pmat06.cpp \
    Source/pmat00.cpp \
    Source/pm.cpp \
    Source/plan94.cpp \
    Source/pfw06.cpp \
    Source/pdp.cpp \
    Source/pb06.cpp \
    Source/pas.cpp \
    Source/pap.cpp \
    Source/p06e.cpp \
    Source/p2s.cpp \
    Source/p2pv.cpp \
    Source/obl80.cpp \
    Source/obl06.cpp \
    Source/nutm80.cpp \
    Source/nut80.cpp \
    Source/nut06a.cpp \
    Source/nut00b.cpp \
    Source/nut00a.cpp \
    Source/numat.cpp \
    Source/num06a.cpp \
    Source/num00b.cpp \
    Source/num00a.cpp \
    Source/ldsun.cpp \
    Source/ldn.cpp \
    Source/ld.cpp \
    Source/jdcalf.cpp \
    Source/jd2cal.cpp \
    Source/ir.cpp \
    Source/hfk5z.cpp \
    Source/h2fk5.cpp \
    Source/gst94.cpp \
    Source/gst06a.cpp \
    Source/gst06.cpp \
    Source/gst00b.cpp \
    Source/gst00a.cpp \
    Source/gmst82.cpp \
    Source/gmst06.cpp \
    Source/gmst00.cpp \
    Source/gd2gce.cpp \
    Source/gd2gc.cpp \
    Source/gc2gde.cpp \
    Source/gc2gd.cpp \
    Source/fw2xy.cpp \
    Source/fw2m.cpp \
    Source/fk52h.cpp \
    Source/fk5hz.cpp \
    Source/fk5hip.cpp \
    Source/fave03.cpp \
    Source/faur03.cpp \
    Source/fasa03.cpp \
    Source/fapa03.cpp \
    Source/faom03.cpp \
    Source/fane03.cpp \
    Source/fame03.cpp \
    Source/fama03.cpp \
    Source/falp03.cpp \
    Source/fal03.cpp \
    Source/faju03.cpp \
    Source/faf03.cpp \
    Source/fae03.cpp \
    Source/fad03.cpp \
    Source/era00.cpp \
    Source/eqeq94.cpp \
    Source/epv00.cpp \
    Source/epj2jd.cpp \
    Source/epj.cpp \
    Source/epb2jd.cpp \
    Source/epb.cpp \
    Source/eors.cpp \
    Source/eo06a.cpp \
    Source/eform.cpp \
    Source/eect00.cpp \
    Source/ee06a.cpp \
    Source/ee00b.cpp \
    Source/ee00a.cpp \
    Source/ee00.cpp \
    Source/dtf2d.cpp \
    Source/dtdb.cpp \
    Source/dat.cpp \
    Source/d2tf.cpp \
    Source/d2dtf.cpp \
    Source/cr.cpp \
    Source/cpv.cpp \
    Source/cp.cpp \
    Source/cal2jd.cpp \
    Source/c2txy.cpp \
    Source/c2tpe.cpp \
    Source/c2teqx.cpp \
    Source/c2tcio.cpp \
    Source/c2t06a.cpp \
    Source/c2t00b.cpp \
    Source/c2t00a.cpp \
    Source/c2s.cpp \
    Source/c2ixys.cpp \
    Source/c2ixy.cpp \
    Source/c2ibpn.cpp \
    Source/c2i06a.cpp \
    Source/c2i00b.cpp \
    Source/c2i00a.cpp \
    Source/bpn2xy.cpp \
    Source/bp06.cpp \
    Source/bp00.cpp \
    Source/bi00.cpp \
    Source/atoiq.cpp \
    Source/atoi13.cpp \
    Source/atoc13.cpp \
    Source/atioq.cpp \
    Source/atio13.cpp \
    Source/aticqn.cpp \
    Source/aticq.cpp \
    Source/atic13.cpp \
    Source/atco13.cpp \
    Source/atciqz.cpp \
    Source/atciqn.cpp \
    Source/atciq.cpp \
    Source/atci13.cpp \
    Source/apio13.cpp \
    Source/apio.cpp \
    Source/aper13.cpp \
    Source/aper.cpp \
    Source/apcs13.cpp \
    Source/apcs.cpp \
    Source/apco13.cpp \
    Source/apco.cpp \
    Source/apci13.cpp \
    Source/apci.cpp \
    Source/apcg13.cpp \
    Source/apcg.cpp \
    Source/anpm.cpp \
    Source/anp.cpp \
    Source/af2a.cpp \
    Source/ab.cpp \
    Source/a2tf.cpp \
    Source/a2af.cpp \
    Source/g2icrs.cpp \
    Source/icrs2g.cpp \
    Source/eceq06.cpp \
    Source/ecm06.cpp \
    Source/eqec06.cpp \
    Source/lteceq.cpp \
    Source/ltecm.cpp \
    Source/lteqec.cpp \
    Source/ltp.cpp \
    Source/ltpb.cpp \
    Source/ltpecl.cpp \
    Source/ltpequ.cpp
