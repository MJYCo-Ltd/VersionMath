#ifndef SOFAHDEF
#define SOFAHDEF

/*
**  - - - - - - -
**   s o f a . h
**  - - - - - - -
**
**  Prototype function declarations for SOFA library.
**
**  This file is part of the International Astronomical Union's
**  SOFA (Standards Of Fundamental Astronomy) software collection.
**
**  This revision:   2015 January 28
**
**  SOFA release 2015-02-09
**
**  Copyright (C) 2015 IAU SOFA Board.  See notes at end.
**
** 注意：如果本文件没有说明，向量对应的 笛卡尔坐标系下的向量
*/
#include "Sofa_global.h"
#include "sofam.h"
#include "math.h"

#ifdef __cplusplus
extern "C" {
#endif

/* Astronomy/历法部分 */

/**
 * @brief 将时分秒转换成天
 * @param s      正负号
 * @param ihour  时 [0~23]
 * @param imin   分 [0~59]
 * @param sec    秒 [0~59.99999]
 * @param days   对应的天数
 * @return 0表示没有错误
 */
SOFASHARED_EXPORT int  iauTf2d(char s, int ihour, int imin, double sec, double *days);
SOFASHARED_EXPORT void iauD2tf(int ndp, double days, char *sign, int ihmsf[4]);

/** 通过年月日计算儒略日
 * @param iy 年
 * @param im 月
 * @param id 日
 * @param djm0 约简儒略日起始天数
 * @param djm  从约简儒略日开始的天数
 */
SOFASHARED_EXPORT int iauCal2jd(int iy, int im, int id, double *djm0, double *djm);
SOFASHARED_EXPORT int iauJd2cal(double dj1, double dj2,
                                int *iy, int *im, int *id, double *fd);
SOFASHARED_EXPORT int iauJdcalf(int ndp, double dj1, double dj2, int iymdf[4]);

/// 精确的儒略日计算方法(包含闰秒，或者说包含跳秒)
SOFASHARED_EXPORT int iauD2dtf(const char *scale, int ndp, double d1, double d2,
                               int *iy, int *im, int *id, int ihmsf[4]);
/// 精确的将儒略日转换成 日期
SOFASHARED_EXPORT int iauDtf2d(const char *scale, int iy, int im, int id,
                               int ihr, int imn, double sec, double *d1, double *d2);


/// 将儒略日转换成 贝塞尔历元（1984年之前使用）
SOFASHARED_EXPORT double iauEpb(double dj1, double dj2);
SOFASHARED_EXPORT void iauEpb2jd(double epb, double *djm0, double *djm);

/// 将儒略日转换成 儒略日历元（现在的标准历元）
SOFASHARED_EXPORT double iauEpj(double dj1, double dj2);
SOFASHARED_EXPORT void iauEpj2jd(double epj, double *djm0, double *djm);

/* Astronomy/时间转换——时间尺度 */
/// 近似计算 TDB-TT 的秒数
/// @param ut    当天的不要天数的UT1儒略日时间
/// @param elong 经度向东为正                [rad]
/// @param u     距离自转轴的距离             [km]
/// @param v     距离赤道到的距离 北为正       [km]
/// @return TDB-TT秒数                      [s]
SOFASHARED_EXPORT double iauDtdb(double date1, double date2,
                                 double ut, double elong, double u, double v);

/// TAI和TT时间转换
/// @brief TT-TAI = 32.184秒 见 sofam.h 中 TTMTAI
SOFASHARED_EXPORT int iauTaitt(double tai1, double tai2, double *tt1, double *tt2);
SOFASHARED_EXPORT int iauTttai(double tt1, double tt2, double *tai1, double *tai2);

/// TAI和UT1时间转换
/// @param dta UT1-TAI 的秒数
SOFASHARED_EXPORT int iauTaiut1(double tai1, double tai2, double dta,
                                double *ut11, double *ut12);
SOFASHARED_EXPORT int iauUt1tai(double ut11, double ut12, double dta,
                                double *tai1, double *tai2);

/// 获取 TAI-UTC 的秒数
/// @brief TAI和UTC相差整数秒，IERS 公报C进行公布
/// @return deltat TAI-UTC
SOFASHARED_EXPORT int iauDat(int iy, int im, int id, double fd, double *deltat);

/// TAI和UTC时间转换
SOFASHARED_EXPORT int iauTaiutc(double tai1, double tai2, double *utc1, double *utc2);
SOFASHARED_EXPORT int iauUtctai(double utc1, double utc2, double *tai1, double *tai2);

/// TCB和TDB时间转换
SOFASHARED_EXPORT int iauTcbtdb(double tcb1, double tcb2, double *tdb1, double *tdb2);
SOFASHARED_EXPORT int iauTdbtcb(double tdb1, double tdb2, double *tcb1, double *tcb2);

/// TCG和TT时间转换
SOFASHARED_EXPORT int iauTcgtt(double tcg1, double tcg2, double *tt1, double *tt2);
SOFASHARED_EXPORT int iauTttcg(double tt1, double tt2, double *tcg1, double *tcg2);

/// TT和TDB时间转换
/// @param dtr TDB-TT 的秒数
SOFASHARED_EXPORT int iauTttdb(double tt1, double tt2, double dtr,
                               double *tdb1, double *tdb2);
SOFASHARED_EXPORT int iauTdbtt(double tdb1, double tdb2, double dtr,
                               double *tt1, double *tt2);

/// TT和UT1时间转换
/// @param dt TT-UT1 的秒数
SOFASHARED_EXPORT int iauTtut1(double tt1, double tt2, double dt,
                               double *ut11, double *ut12);
SOFASHARED_EXPORT int iauUt1tt(double ut11, double ut12, double dt,
                               double *tt1, double *tt2);

/// UT1和UTC时间转换
/// @param dut1 UT1-UTC 的秒数
/// @brief 数据可以从IERS公报B得到
SOFASHARED_EXPORT int iauUt1utc(double ut11, double ut12, double dut1,
                                double *utc1, double *utc2);
SOFASHARED_EXPORT int iauUtcut1(double utc1, double utc2, double dut1,
                                double *ut11, double *ut12);
/* Astronomy/地球自转和恒星时 */
//////////////////////瞬时春分点与平春分点的计算/////////////
/// GAST = GMST + EE + CT(可以不要)
/// GST、GAST 两个名词等同
/// GAST 是ECI转ECF时需要用的 时角

/// @brief 瞬时真春分点与平春分点的差值
/// @param date1 date2 date1+date2=JD(TT)
/// @param epsa 黄赤交角              [rad]
/// @param dpsi 黄经章动              [rad]
/// @return 瞬时平春分点到平春分点的角度  [rad]
SOFASHARED_EXPORT double iauEe00(double date1, double date2, double epsa, double dpsi);
SOFASHARED_EXPORT double iauEe00a(double date1, double date2);
SOFASHARED_EXPORT double iauEe00b(double date1, double date2);
SOFASHARED_EXPORT double iauEe06a(double date1, double date2);

/// @brief 瞬时真春分点与平春分点的差值 使用的补充项
/// @brief GAST = GMST + CT +　EE
/// @param date1 date2 date1+date2=JD(TT)
/// @return 瞬时平春分点到平春分点的角度补充项 [rad]
SOFASHARED_EXPORT double iauEect00(double date1, double date2);
/// @param date1 date2 date1+date2=JD(TDB)
/// @return 瞬时平春分点到平春分点的角度      [rad]
SOFASHARED_EXPORT double iauEqeq94(double date1, double date2);
//////////////////////瞬时春分点与平春分点的计算 end/////////////

/// @brief 地球自转角
/// @param dj1 dj2 dj1+dj2=JD(UT1)
/// @return 地球自转角               [rad]
SOFASHARED_EXPORT double iauEra00(double dj1, double dj2);

/// @brief 格林尼治平恒星时(GMST)：真春分点和格林尼治子午线之间的夹角
/// @param uta utb uta+utb=JD(UT1)
/// @param tta ttb tta+ttb=JD(TT)
/// @return GMST                   [rad]
SOFASHARED_EXPORT double iauGmst00(double uta, double utb, double tta, double ttb);
SOFASHARED_EXPORT double iauGmst06(double uta, double utb, double tta, double ttb);
/// @param dj1 dj2 dj1+dj2=JD(UT1)
SOFASHARED_EXPORT double iauGmst82(double dj1, double dj2);

/// @brief 格林尼治真恒星时(GAST):测量的是真春分点的时角
/// @param uta utb uta+utb=JD(UT1)
/// @param tta ttb tta+ttb=JD(TT)
/// @return GAST                   [rad]
SOFASHARED_EXPORT double iauGst00a(double uta, double utb, double tta, double ttb);
SOFASHARED_EXPORT double iauGst00b(double uta, double utb);
SOFASHARED_EXPORT double iauGst06(double uta, double utb, double tta, double ttb,
                                  double rnpb[3][3]);
SOFASHARED_EXPORT double iauGst06a(double uta, double utb, double tta, double ttb);
SOFASHARED_EXPORT double iauGst94(double uta, double utb);

/* Astronomy/EclipticCoordinates */
SOFASHARED_EXPORT void iauEceq06(double date1, double date2, double dl, double db,
               double *dr, double *dd);
SOFASHARED_EXPORT void iauEcm06(double date1, double date2, double rm[3][3]);
SOFASHARED_EXPORT void iauEqec06(double date1, double date2, double dr, double dd,
               double *dl, double *db);
SOFASHARED_EXPORT void iauLteceq(double epj, double dl, double db, double *dr, double *dd);
SOFASHARED_EXPORT void iauLtecm(double epj, double rm[3][3]);
SOFASHARED_EXPORT void iauLteqec(double epj, double dr, double dd, double *dl, double *db);

/* Astronomy/银河系坐标 */
SOFASHARED_EXPORT void iauG2icrs ( double dl, double db, double *dr, double *dd );
SOFASHARED_EXPORT void iauIcrs2g ( double dr, double dd, double *dl, double *db );

/* Astronomy/地球地理 */
/// @brief 获取地球椭球系数
/// @param n 1(WGS84)2(GRS80)3(WGS72)
/// @param a 长半轴          [m]
/// @param f 扁率
SOFASHARED_EXPORT int iauEform(int n, double *a, double *f);

/// @brief 地心坐标系转换成地理坐标系
/// @param n 1(WGS84)2(GRS80)3(WGS72)
/// @param xyz 单位        [m]
/// @param elong  经度     [rad] (向东为正)
/// @param phi    纬度     [rad] (向北为正)
/// @param height 高度     [m]
SOFASHARED_EXPORT int iauGc2gd(int n, double xyz[3],
                               double *elong, double *phi, double *height);
/// @param a 长半轴
/// @param f 扁率
/// @attention a xyz height 必须是相同的单位 米 千米 等
SOFASHARED_EXPORT int iauGc2gde(double a, double f, double xyz[3],
                                double *elong, double *phi, double *height);

/// @brief 地理坐标系转换成地心坐标系
SOFASHARED_EXPORT int iauGd2gc(int n, double elong, double phi, double height,
                               double xyz[3]);
SOFASHARED_EXPORT int iauGd2gce(double a, double f,
                                double elong, double phi, double height, double xyz[3]);
/// @param elong  WGS84 下的经度 [rad]
/// @param phi    WGS84 下的纬度 [rad]
/// @param height 高度          [m]
/// @param xp,yp 极移数据        [rad]
/// @param sp    TIO定位角       [rad]    由章动可以计算见 iauS06a
/// @param theta 地球自转角       [rad]
/// @param pv    CIRS坐标系下 观察者的位置速度 单位 m m/s CIRS(天球中间参考系，对GCRS进行岁差章转动后的地心参考系)
SOFASHARED_EXPORT void iauPvtob(double elong, double phi, double height, double xp,
                                double yp, double sp, double theta, double pv[2][3]);
/* Astronomy/岁差 章动 极移 */
/// GCRS 到 ITRS 坐标系转换矩阵 (包括岁差、章动、极移)
/// @param tta tta+ttb = dJD(TT)
/// @param uta uta+utb = dJD(UT1)
/// @param xp,yp  极移数据     [rad]
SOFASHARED_EXPORT void iauC2t00a(double tta, double ttb, double uta, double utb,
                                 double xp, double yp, double rc2t[3][3]);
SOFASHARED_EXPORT void iauC2t00b(double tta, double ttb, double uta, double utb,
                                 double xp, double yp, double rc2t[3][3]);
SOFASHARED_EXPORT void iauC2t06a(double tta, double ttb, double uta, double utb,
                                 double xp, double yp, double rc2t[3][3]);

/// @brief IAU 2000A model
/// @param date1 TT 儒略日时间
/// @param date2 TT 儒略日时间 date1+date2 = dJD(TT)
/// 精度为 0.2 mas
SOFASHARED_EXPORT void iauNum00a(double date1, double date2, double rmatn[3][3]);
/// @attention 相较于 num00a 00b速度快，精度为1mas（1毫弧秒或者1毫角秒）
SOFASHARED_EXPORT void iauNum00b(double date1, double date2, double rmatn[3][3]);
SOFASHARED_EXPORT void iauNum06a(double date1, double date2, double rmatn[3][3]);
SOFASHARED_EXPORT void iauNutm80(double date1, double date2, double rmatn[3][3]);

/// @brief IAU 2000A model
/// @param date1 date2 date1+date2=JD(TT)
SOFASHARED_EXPORT void iauNut00a(double date1, double date2, double *dpsi, double *deps);
SOFASHARED_EXPORT void iauNut00b(double date1, double date2, double *dpsi, double *deps);
SOFASHARED_EXPORT void iauNut06a(double date1, double date2, double *dpsi, double *deps);
SOFASHARED_EXPORT void iauNut80(double date1, double date2, double *dpsi, double *deps);


SOFASHARED_EXPORT void iauBi00(double *dpsibi, double *depsbi, double *dra);
SOFASHARED_EXPORT void iauBp00(double date1, double date2,
                               double rb[3][3], double rp[3][3], double rbp[3][3]);
SOFASHARED_EXPORT void iauBp06(double date1, double date2,
                               double rb[3][3], double rp[3][3], double rbp[3][3]);
SOFASHARED_EXPORT void iauBpn2xy(double rbpn[3][3], double *x, double *y);
SOFASHARED_EXPORT void iauC2i00a(double date1, double date2, double rc2i[3][3]);
SOFASHARED_EXPORT void iauC2i00b(double date1, double date2, double rc2i[3][3]);
SOFASHARED_EXPORT void iauC2i06a(double date1, double date2, double rc2i[3][3]);
SOFASHARED_EXPORT void iauC2ibpn(double date1, double date2, double rbpn[3][3],
                                 double rc2i[3][3]);
SOFASHARED_EXPORT void iauC2ixy(double date1, double date2, double x, double y,
                                double rc2i[3][3]);
SOFASHARED_EXPORT void iauC2ixys(double x, double y, double s, double rc2i[3][3]);
SOFASHARED_EXPORT void iauC2tcio(double rc2i[3][3], double era, double rpom[3][3],
                                 double rc2t[3][3]);
SOFASHARED_EXPORT void iauC2teqx(double rbpn[3][3], double gst, double rpom[3][3],
                                 double rc2t[3][3]);
SOFASHARED_EXPORT void iauC2tpe(double tta, double ttb, double uta, double utb,
                                double dpsi, double deps, double xp, double yp,
                                double rc2t[3][3]);
SOFASHARED_EXPORT void iauC2txy(double tta, double ttb, double uta, double utb,
                                double x, double y, double xp, double yp,
                                double rc2t[3][3]);
SOFASHARED_EXPORT double iauEo06a(double date1, double date2);
SOFASHARED_EXPORT double iauEors(double rnpb[3][3], double s);
SOFASHARED_EXPORT void iauFw2m(double gamb, double phib, double psi, double eps,
                               double r[3][3]);
SOFASHARED_EXPORT void iauFw2xy(double gamb, double phib, double psi, double eps,
                                double *x, double *y);
SOFASHARED_EXPORT void iauLtp(double epj, double rp[3][3]);
SOFASHARED_EXPORT void iauLtpb(double epj, double rpb[3][3]);
SOFASHARED_EXPORT void iauLtpecl(double epj, double vec[3]);
SOFASHARED_EXPORT void iauLtpequ(double epj, double veq[3]);
SOFASHARED_EXPORT void iauNumat(double epsa, double dpsi, double deps, double rmatn[3][3]);
SOFASHARED_EXPORT double iauObl06(double date1, double date2);
SOFASHARED_EXPORT double iauObl80(double date1, double date2);
SOFASHARED_EXPORT void iauP06e(double date1, double date2,
                               double *eps0, double *psia, double *oma, double *bpa,
                               double *bqa, double *pia, double *bpia,
                               double *epsa, double *chia, double *za, double *zetaa,
                               double *thetaa, double *pa,
                               double *gam, double *phi, double *psi);
SOFASHARED_EXPORT void iauPb06(double date1, double date2,
                               double *bzeta, double *bz, double *btheta);
SOFASHARED_EXPORT void iauPfw06(double date1, double date2,
                                double *gamb, double *phib, double *psib, double *epsa);
SOFASHARED_EXPORT void iauPmat00(double date1, double date2, double rbp[3][3]);
SOFASHARED_EXPORT void iauPmat06(double date1, double date2, double rbp[3][3]);
SOFASHARED_EXPORT void iauPmat76(double date1, double date2, double rmatp[3][3]);
SOFASHARED_EXPORT void iauPn00(double date1, double date2, double dpsi, double deps,
                               double *epsa,
                               double rb[3][3], double rp[3][3], double rbp[3][3],
                               double rn[3][3], double rbpn[3][3]);
SOFASHARED_EXPORT void iauPn00a(double date1, double date2,
                                double *dpsi, double *deps, double *epsa,
                                double rb[3][3], double rp[3][3], double rbp[3][3],
                                double rn[3][3], double rbpn[3][3]);
SOFASHARED_EXPORT void iauPn00b(double date1, double date2,
                                double *dpsi, double *deps, double *epsa,
                                double rb[3][3], double rp[3][3], double rbp[3][3],
                                double rn[3][3], double rbpn[3][3]);
SOFASHARED_EXPORT void iauPn06(double date1, double date2, double dpsi, double deps,
                               double *epsa,
                               double rb[3][3], double rp[3][3], double rbp[3][3],
                               double rn[3][3], double rbpn[3][3]);
SOFASHARED_EXPORT void iauPn06a(double date1, double date2,
                                double *dpsi, double *deps, double *epsa,
                                double rb[3][3], double rp[3][3], double rbp[3][3],
                                double rn[3][3], double rbpn[3][3]);
SOFASHARED_EXPORT void iauPnm00a(double date1, double date2, double rbpn[3][3]);
SOFASHARED_EXPORT void iauPnm00b(double date1, double date2, double rbpn[3][3]);
SOFASHARED_EXPORT void iauPnm06a(double date1, double date2, double rnpb[3][3]);
SOFASHARED_EXPORT void iauPnm80(double date1, double date2, double rmatpn[3][3]);
SOFASHARED_EXPORT void iauPom00(double xp, double yp, double sp, double rpom[3][3]);
SOFASHARED_EXPORT void iauPr00(double date1, double date2, double *dpsipr, double *depspr);
SOFASHARED_EXPORT void iauPrec76(double ep01, double ep02, double ep11, double ep12,
               double *zeta, double *z, double *theta);
SOFASHARED_EXPORT double iauS00(double date1, double date2, double x, double y);
SOFASHARED_EXPORT double iauS00a(double date1, double date2);
SOFASHARED_EXPORT double iauS00b(double date1, double date2);
SOFASHARED_EXPORT double iauS06(double date1, double date2, double x, double y);
SOFASHARED_EXPORT double iauS06a(double date1, double date2);
SOFASHARED_EXPORT double iauSp00(double date1, double date2);
SOFASHARED_EXPORT void iauXy06(double date1, double date2, double *x, double *y);
SOFASHARED_EXPORT void iauXys00a(double date1, double date2,
               double *x, double *y, double *s);
SOFASHARED_EXPORT void iauXys00b(double date1, double date2,
               double *x, double *y, double *s);
SOFASHARED_EXPORT void iauXys06a(double date1, double date2,
               double *x, double *y, double *s);

/* Astronomy/天体测量 */
SOFASHARED_EXPORT void iauAb(double pnat[3], double v[3], double s, double bm1,
                             double ppr[3]);
SOFASHARED_EXPORT void iauApcg(double date1, double date2,
                               double ebpv[2][3], double ehp[3],
                               iauASTROM *astrom);
SOFASHARED_EXPORT void iauApcg13(double date1, double date2, iauASTROM *astrom);
SOFASHARED_EXPORT void iauApci(double date1, double date2,
                               double ebpv[2][3], double ehp[3],
                               double x, double y, double s,
                               iauASTROM *astrom);
SOFASHARED_EXPORT void iauApci13(double date1, double date2,
                                 iauASTROM *astrom, double *eo);
SOFASHARED_EXPORT void iauApco(double date1, double date2,
                               double ebpv[2][3], double ehp[3],
                               double x, double y, double s, double theta,
                               double elong, double phi, double hm,
                               double xp, double yp, double sp,
                               double refa, double refb,
                               iauASTROM *astrom);
SOFASHARED_EXPORT int iauApco13(double utc1, double utc2, double dut1,
                                double elong, double phi, double hm, double xp, double yp,
                                double phpa, double tk, double rh, double wl,
                                iauASTROM *astrom, double *eo);
SOFASHARED_EXPORT void iauApcs(double date1, double date2, double pv[2][3],
                               double ebpv[2][3], double ehp[3],
                               iauASTROM *astrom);
SOFASHARED_EXPORT void iauApcs13(double date1, double date2, double pv[2][3],
                                 iauASTROM *astrom);
SOFASHARED_EXPORT void iauAper(double theta, iauASTROM *astrom);
SOFASHARED_EXPORT void iauAper13(double ut11, double ut12, iauASTROM *astrom);
SOFASHARED_EXPORT void iauApio(double sp, double theta,
                               double elong, double phi, double hm, double xp, double yp,
                               double refa, double refb,
                               iauASTROM *astrom);
SOFASHARED_EXPORT int iauApio13(double utc1, double utc2, double dut1,
                                double elong, double phi, double hm, double xp, double yp,
                                double phpa, double tk, double rh, double wl,
                                iauASTROM *astrom);
SOFASHARED_EXPORT void iauAtci13(double rc, double dc,
                                 double pr, double pd, double px, double rv,
                                 double date1, double date2,
                                 double *ri, double *di, double *eo);
SOFASHARED_EXPORT void iauAtciq(double rc, double dc, double pr, double pd,
                                double px, double rv, iauASTROM *astrom,
                                double *ri, double *di);
SOFASHARED_EXPORT void iauAtciqn(double rc, double dc, double pr, double pd,
                                 double px, double rv, iauASTROM *astrom,
                                 int n, iauLDBODY b[], double *ri, double *di);
SOFASHARED_EXPORT void iauAtciqz(double rc, double dc, iauASTROM *astrom,
                                 double *ri, double *di);
SOFASHARED_EXPORT int iauAtco13(double rc, double dc,
                                double pr, double pd, double px, double rv,
                                double utc1, double utc2, double dut1,
                                double elong, double phi, double hm, double xp, double yp,
                                double phpa, double tk, double rh, double wl,
                                double *aob, double *zob, double *hob,
                                double *dob, double *rob, double *eo);
SOFASHARED_EXPORT void iauAtic13(double ri, double di,
                                 double date1, double date2,
                                 double *rc, double *dc, double *eo);
SOFASHARED_EXPORT void iauAticq(double ri, double di, iauASTROM *astrom,
                                double *rc, double *dc);
SOFASHARED_EXPORT void iauAticqn(double ri, double di, iauASTROM *astrom,
                                 int n, iauLDBODY b[], double *rc, double *dc);
SOFASHARED_EXPORT int iauAtio13(double ri, double di,
                                double utc1, double utc2, double dut1,
                                double elong, double phi, double hm, double xp, double yp,
                                double phpa, double tk, double rh, double wl,
                                double *aob, double *zob, double *hob,
                                double *dob, double *rob);
SOFASHARED_EXPORT void iauAtioq(double ri, double di, iauASTROM *astrom,
                                double *aob, double *zob,
                                double *hob, double *dob, double *rob);
SOFASHARED_EXPORT int iauAtoc13(const char *type, double ob1, double ob2,
                                double utc1, double utc2, double dut1,
                                double elong, double phi, double hm, double xp, double yp,
                                double phpa, double tk, double rh, double wl,
                                double *rc, double *dc);
SOFASHARED_EXPORT int iauAtoi13(const char *type, double ob1, double ob2,
                                double utc1, double utc2, double dut1,
                                double elong, double phi, double hm, double xp, double yp,
                                double phpa, double tk, double rh, double wl,
                                double *ri, double *di);
SOFASHARED_EXPORT void iauAtoiq(const char *type,
                                double ob1, double ob2, iauASTROM *astrom,
                                double *ri, double *di);
SOFASHARED_EXPORT void iauLd(double bm, double p[3], double q[3], double e[3],
                             double em, double dlim, double p1[3]);
SOFASHARED_EXPORT void iauLdn(int n, iauLDBODY b[], double ob[3], double sc[3],
                              double sn[3]);
SOFASHARED_EXPORT void iauLdsun(double p[3], double e[3], double em, double p1[3]);
SOFASHARED_EXPORT void iauPmpx(double rc, double dc, double pr, double pd,
                               double px, double rv, double pmt, double vob[3],
                               double pco[3]);
SOFASHARED_EXPORT int iauPmsafe(double ra1, double dec1, double pmr1, double pmd1,
                                double px1, double rv1,
                                double ep1a, double ep1b, double ep2a, double ep2b,
                                double *ra2, double *dec2, double *pmr2, double *pmd2,
                                double *px2, double *rv2);
SOFASHARED_EXPORT void iauRefco(double phpa, double tk, double rh, double wl,
                                double *refa, double *refb);

/* Astronomy/星历表 */
/// @brief 地球在太阳系中的位置和速度
/// @param date1 + date2 = dJD(TDB)
/// @param pvh 以太阳为中心的坐标系的 位置、速度 [Au、Au/d]
/// @param pvb 以太阳质心为中心的坐标系的位置、速度 [Au、Au/d]
/// @attention 只能计算 1900-2100 年
SOFASHARED_EXPORT int iauEpv00(double date1, double date2,
             double pvh[2][3], double pvb[2][3]);

/// @brief 各行星在太阳系的位置和速度
/// @param date1 + date2 = dJD(TDB)
/// @param np [1-8] [Mercury,Venus,EMB,Mars,Jupiter,Saturn,Uranus,Neptune]
/// @param pv 以太阳为中心的坐标系的 位置、速度 [Au、Au/d]
/// @attention 只能计算 1000-3000 年
SOFASHARED_EXPORT int iauPlan94(double date1, double date2, int np, double pv[2][3]);

/* Astronomy/用于章动等公式中的基本幅角 */
SOFASHARED_EXPORT double iauFad03(double t);
SOFASHARED_EXPORT double iauFae03(double t);
SOFASHARED_EXPORT double iauFaf03(double t);
SOFASHARED_EXPORT double iauFaju03(double t);
SOFASHARED_EXPORT double iauFal03(double t);
SOFASHARED_EXPORT double iauFalp03(double t);
SOFASHARED_EXPORT double iauFama03(double t);
SOFASHARED_EXPORT double iauFame03(double t);
SOFASHARED_EXPORT double iauFane03(double t);
SOFASHARED_EXPORT double iauFaom03(double t);
SOFASHARED_EXPORT double iauFapa03(double t);
SOFASHARED_EXPORT double iauFasa03(double t);
SOFASHARED_EXPORT double iauFaur03(double t);
SOFASHARED_EXPORT double iauFave03(double t);

/* Astronomy/SpaceMotion */
SOFASHARED_EXPORT int iauPmsafe(double ra1, double dec1, double pmr1, double pmd1,
              double px1, double rv1, double ep1a, double ep1b,
              double ep2a, double ep2b, double *ra2, double *dec2,
              double *pmr2, double *pmd2, double *px2, double *rv2);
SOFASHARED_EXPORT int iauPvstar(double pv[2][3], double *ra, double *dec,
              double *pmr, double *pmd, double *px, double *rv);
SOFASHARED_EXPORT int iauStarpv(double ra, double dec,
              double pmr, double pmd, double px, double rv,
              double pv[2][3]);

/* Astronomy/星表转换 */
SOFASHARED_EXPORT void iauFk52h(double r5, double d5,
              double dr5, double dd5, double px5, double rv5,
              double *rh, double *dh,
              double *drh, double *ddh, double *pxh, double *rvh);
SOFASHARED_EXPORT void iauFk5hip(double r5h[3][3], double s5h[3]);
SOFASHARED_EXPORT void iauFk5hz(double r5, double d5, double date1, double date2,
              double *rh, double *dh);
SOFASHARED_EXPORT void iauH2fk5(double rh, double dh,
              double drh, double ddh, double pxh, double rvh,
              double *r5, double *d5,
              double *dr5, double *dd5, double *px5, double *rv5);
SOFASHARED_EXPORT void iauHfk5z(double rh, double dh, double date1, double date2,
              double *r5, double *d5, double *dr5, double *dd5);
SOFASHARED_EXPORT int iauStarpm(double ra1, double dec1,
              double pmr1, double pmd1, double px1, double rv1,
              double ep1a, double ep1b, double ep2a, double ep2b,
              double *ra2, double *dec2,
              double *pmr2, double *pmd2, double *px2, double *rv2);

/* VectorMatrix/AngleOps */

/**
 * @brief 将度分秒转成弧度
 * @param s      正负号
 * @param ideg   度   [0~359]
 * @param iamin  分   [0~59]
 * @param asec   秒   [0~59.9999]
 * @param rad    计算好的弧度值
 * @return 0     表示计算完成
 */
SOFASHARED_EXPORT int  iauAf2a(char s, int ideg, int iamin, double asec, double *rad);
SOFASHARED_EXPORT void iauA2af(int ndp, double angle, char *sign, int idmsf[4]);

/**
 * @brief 将时分秒转成地球旋转角度
 * @param s      正负号
 * @param ihour  时 [0~23]
 * @param imin   分 [0~59]
 * @param sec    秒 [0~59.99999]
 * @param rad    地球旋转的弧度值 [rad]
 * @return 0表示转换成功
 */
SOFASHARED_EXPORT int  iauTf2a(char s, int ihour, int imin, double sec, double *rad);
SOFASHARED_EXPORT void iauA2tf(int ndp, double angle, char *sign, int ihmsf[4]);

/**
 * @brief iauAnp  将弧度限定到 [0~2π)
 * @brief iauAnpm 将弧度限定到 [-π~π)
 * @param a 任意弧度值
 * @return
 */
SOFASHARED_EXPORT double iauAnp(double a);
SOFASHARED_EXPORT double iauAnpm(double a);


/* VectorMatrix/BuildRotations */
SOFASHARED_EXPORT void iauRx(double phi, double r[3][3]);
SOFASHARED_EXPORT void iauRy(double theta, double r[3][3]);
SOFASHARED_EXPORT void iauRz(double psi, double r[3][3]);

/* VectorMatrix/CopyExtendExtract */
SOFASHARED_EXPORT void iauCp(double p[3], double c[3]);
SOFASHARED_EXPORT void iauCpv(double pv[2][3], double c[2][3]);
SOFASHARED_EXPORT void iauCr(double r[3][3], double c[3][3]);
SOFASHARED_EXPORT void iauP2pv(double p[3], double pv[2][3]);
SOFASHARED_EXPORT void iauPv2p(double pv[2][3], double p[3]);

/* VectorMatrix/Initialization */
SOFASHARED_EXPORT void iauIr(double r[3][3]);
SOFASHARED_EXPORT void iauZp(double p[3]);
SOFASHARED_EXPORT void iauZpv(double pv[2][3]);
SOFASHARED_EXPORT void iauZr(double r[3][3]);

/* VectorMatrix/MatrixOps */
SOFASHARED_EXPORT void iauRxr(double a[3][3], double b[3][3], double atb[3][3]);
SOFASHARED_EXPORT void iauTr(double r[3][3], double rt[3][3]);

/* VectorMatrix/MatrixVectorProducts */
SOFASHARED_EXPORT void iauRxp(double r[3][3], double p[3], double rp[3]);
SOFASHARED_EXPORT void iauRxpv(double r[3][3], double pv[2][3], double rpv[2][3]);
SOFASHARED_EXPORT void iauTrxp(double r[3][3], double p[3], double trp[3]);
SOFASHARED_EXPORT void iauTrxpv(double r[3][3], double pv[2][3], double trpv[2][3]);

/* VectorMatrix/四元数与矩阵转换 */
/// @brief 将矩阵转换成四元数
/// @param r 旋转矩阵
/// @return w w的模是旋转角度 [rad] 方向为旋转轴的方向
SOFASHARED_EXPORT void iauRm2v(double r[3][3], double w[3]);

/// @brief 将四元数转换成旋转矩阵
/// @return r r为旋转矩阵
SOFASHARED_EXPORT void iauRv2m(double w[3], double r[3][3]);

/* VectorMatrix/方位角与角距 */
/// @brief b点相对于a点的方位角
/// @param a 向量 a
/// @param b 向量 b
/// @return 方位角 [rad]
SOFASHARED_EXPORT double iauPap(double a[3], double b[3]);

/// @brief b 点相对于 a点的方位角 球面坐标
/// @param al a 点的经度 [rad]
/// @param ap a 点的纬度 [rad]
/// @param bl b 点的经度 [rad]
/// @param bp b 点的纬度 [rad]
/// @return 方位角       [rad]
SOFASHARED_EXPORT double iauPas(double al, double ap, double bl, double bp);

/// @brief 向量a 与 向量b 的角距(夹角)
/// @param a 向量 a
/// @param b 向量 b
/// @return 角距        [rad]
/// @attention 经过处理在 0 和 Pi 时精度依然很高
SOFASHARED_EXPORT double iauSepp(double a[3], double b[3]);

/// @brief 向量a 与 向量b 的角距(夹角) 球面坐标
/// @param al a 点的经度 [rad]
/// @param ap a 点的纬度 [rad]
/// @param bl b 点的经度 [rad]
/// @param bp b 点的纬度 [rad]
/// @return 角距        [rad]
SOFASHARED_EXPORT double iauSeps(double al, double ap, double bl, double bp);

/* VectorMatrix/球坐标与笛卡尔坐标的转换 */
/// @brief 求笛卡尔坐标对应的球坐标
/// @param p 笛卡尔坐标系位置
/// @return theta 经度        [rad]
/// @return phi   纬度        [rad]
SOFASHARED_EXPORT void iauC2s(double p[3], double *theta, double *phi);

/// @brief 求笛卡尔坐标对应的球坐标
/// 参数同上
/// @return r 距离     单位与传入的分量相同
SOFASHARED_EXPORT void iauP2s(double p[3], double *theta, double *phi, double *r);

/// @brief 求笛卡尔下的位置/速度 对应的球坐标及其变化率
/// @param pv 笛卡尔位置/速度
/// @return theta 等同上
/// @return td 经度的变化率  [rad/t] t与速度的单位相同
/// @return pd 纬度的变化率  [rad/t] t与速度的单位相同
/// @return rd 距离的变化率  [r/t] r与位置的分量单位相同 t与速度的单位相同
SOFASHARED_EXPORT void iauPv2s(double pv[2][3],
             double *theta, double *phi, double *r,
             double *td, double *pd, double *rd);

//// 这几个是上面三个的逆变换
SOFASHARED_EXPORT void iauS2c(double theta, double phi, double c[3]);
SOFASHARED_EXPORT void iauS2p(double theta, double phi, double r, double p[3]);
SOFASHARED_EXPORT void iauS2pv(double theta, double phi, double r,
             double td, double pd, double rd,
             double pv[2][3]);

/* VectorMatrix/向量操作 */
/// @brief 求向量a与向量b的点乘(内积)
/// @param 向量 a
/// @param 向量 b
/// @return a·b 向量a与向量b的点乘
SOFASHARED_EXPORT double iauPdp(double a[3], double b[3]);
/// @return adb[0] = a[0] . b[0]
/// @return adb[1] = a[0] . b[1] + a[1] . b[0]
SOFASHARED_EXPORT void iauPvdpv(double a[2][3], double b[2][3], double adb[2]);

/// @brief 求向量p的模
/// @param 向量 p
/// @return 向量p的模
SOFASHARED_EXPORT double iauPm(double p[3]);
/// @return r 向量pv[0]的模
/// @return s 向量pv[1]的模
SOFASHARED_EXPORT void iauPvm(double pv[2][3], double *r, double *s);

/// @brief 求向量a减向量b
/// @param a 向量a
/// @param b 向量b
/// @return amb = a - b
SOFASHARED_EXPORT void iauPmp(double a[3], double b[3], double amb[3]);
/// @return amb[0] = a[0] - b[0]
/// @return amb[1] = a[1] - b[1]
SOFASHARED_EXPORT void iauPvmpv(double a[2][3], double b[2][3], double amb[2][3]);

/// @brief 根据向量求 向量的模 以及归一化的向量
/// @param p  向量 p
/// @return r 向量的模
/// @return u 归一化的向量
SOFASHARED_EXPORT void iauPn(double p[3], double *r, double u[3]);

/// @brief 求向量的和
/// @param a 向量 a
/// @param b 向量 b
/// @return apb = a + b
SOFASHARED_EXPORT void iauPpp(double a[3], double b[3], double apb[3]);
/// @return apb[0] = a[0] + b[0]
/// @return apb[1] = a[1] + b[1]
SOFASHARED_EXPORT void iauPvppv(double a[2][3], double b[2][3], double apb[2][3]);

/// @brief 求向量的和
/// @param a 向量 a
/// @param s 缩放系数
/// @param b 向量 b
/// @return apsb = a + s*b
SOFASHARED_EXPORT void iauPpsp(double a[3], double s, double b[3], double apsb[3]);

/// @brief 根据 位置/速度 时间 更新 位置/速度
/// @param dt  往后推算的时间  时间单位与pv中v的时间尺度相同
/// @param pv  位置/速度      p[0] 位置 p[1] 速度
/// @return upv 更新后的位置/速度 upv[0] = p[0] + p[1]*dt;upv[1] = p[1];
SOFASHARED_EXPORT void iauPvu(double dt, double pv[2][3], double upv[2][3]);
SOFASHARED_EXPORT void iauPvup(double dt, double pv[2][3], double p[3]);

/// @brief 求向量a与向量b的叉乘(外积)
/// @param 向量 a
/// @param 向量 b
/// @return axb = a X b 向量a与向量b的叉乘
SOFASHARED_EXPORT void iauPxp(double a[3], double b[3], double axb[3]);
/// @return axb[0] = a[0] X b[0]
/// @return axb[1] = a[0] X b[1] + a[1] X b[0]
SOFASHARED_EXPORT void iauPvxpv(double a[2][3], double b[2][3], double axb[2][3]);

/// @brief 对向量进行缩放
/// @param s  缩放系数
/// @param p  向量
/// @return sp = s * p 缩放后的向量
SOFASHARED_EXPORT void iauSxp(double s, double p[3], double sp[3]);
/// @return spv = s * pv
SOFASHARED_EXPORT void iauSxpv(double s, double pv[2][3], double spv[2][3]);
/// @return spv[0] = s1 * pv[0]
/// @return spv[1] = s2 * pv[0]
SOFASHARED_EXPORT void iauS2xpv(double s1, double s2, double pv[2][3], double spv[2][3]);

#ifdef __cplusplus
}
#endif

#endif

/*----------------------------------------------------------------------
**
**  Copyright (C) 2015
**  Standards Of Fundamental Astronomy Board
**  of the International Astronomical Union.
**
**  =====================
**  SOFA Software License
**  =====================
**
**  NOTICE TO USER:
**
**  BY USING THIS SOFTWARE YOU ACCEPT THE FOLLOWING SIX TERMS AND
**  CONDITIONS WHICH APPLY TO ITS USE.
**
**  1. The Software is owned by the IAU SOFA Board ("SOFA").
**
**  2. Permission is granted to anyone to use the SOFA software for any
**     purpose, including commercial applications, free of charge and
**     without payment of royalties, subject to the conditions and
**     restrictions listed below.
**
**  3. You (the user) may copy and distribute SOFA source code to others,
**     and use and adapt its code and algorithms in your own software,
**     on a world-wide, royalty-free basis.  That portion of your
**     distribution that does not consist of intact and unchanged copies
**     of SOFA source code files is a "derived work" that must comply
**     with the following requirements:
**
**     a) Your work shall be marked or carry a statement that it
**        (i) uses routines and computations derived by you from
**        software provided by SOFA under license to you; and
**        (ii) does not itself constitute software provided by and/or
**        endorsed by SOFA.
**
**     b) The source code of your derived work must contain descriptions
**        of how the derived work is based upon, contains and/or differs
**        from the original SOFA software.
**
**     c) The names of all routines in your derived work shall not
**        include the prefix "iau" or "sofa" or trivial modifications
**        thereof such as changes of case.
**
**     d) The origin of the SOFA components of your derived work must
**        not be misrepresented;  you must not claim that you wrote the
**        original software, nor file a patent application for SOFA
**        software or algorithms embedded in the SOFA software.
**
**     e) These requirements must be reproduced intact in any source
**        distribution and shall apply to anyone to whom you have
**        granted a further right to modify the source code of your
**        derived work.
**
**     Note that, as originally distributed, the SOFA software is
**     intended to be a definitive implementation of the IAU standards,
**     and consequently third-party modifications are discouraged.  All
**     variations, no matter how minor, must be explicitly marked as
**     such, as explained above.
**
**  4. You shall not cause the SOFA software to be brought into
**     disrepute, either by misuse, or use for inappropriate tasks, or
**     by inappropriate modification.
**
**  5. The SOFA software is provided "as is" and SOFA makes no warranty
**     as to its use or performance.   SOFA does not and cannot warrant
**     the performance or results which the user may obtain by using the
**     SOFA software.  SOFA makes no warranties, express or implied, as
**     to non-infringement of third party rights, merchantability, or
**     fitness for any particular purpose.  In no event will SOFA be
**     liable to the user for any consequential, incidental, or special
**     damages, including any lost profits or lost savings, even if a
**     SOFA representative has been advised of such damages, or for any
**     claim by any third party.
**
**  6. The provision of any version of the SOFA software under the terms
**     and conditions specified herein does not imply that future
**     versions will also be made available under the same terms and
**     conditions.
*
**  In any published work or commercial product which uses the SOFA
**  software directly, acknowledgement (see www.iausofa.org) is
**  appreciated.
**
**  Correspondence concerning SOFA software should be addressed as
**  follows:
**
**      By email:  sofa@ukho.gov.uk
**      By post:   IAU SOFA Center
**                 HM Nautical Almanac Office
**                 UK Hydrographic Office
**                 Admiralty Way, Taunton
**                 Somerset, TA1 2DN
**                 United Kingdom
**
**--------------------------------------------------------------------*/
