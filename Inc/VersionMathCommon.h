#ifndef _INCLUDE_VERSION_MATH_COMMON_H
#define _INCLUDE_VERSION_MATH_COMMON_H

/* T0 */
const double T0(-273.20);

/* Pi */
const double DPI(3.141592653589793238462643);

/* 2Pi */
const double D2PI(6.283185307179586476925287);

/* Radians to degrees */
const double DR2D (57.29577951308232087679815);

/* Degrees to radians */
const double DD2R(1.745329251994329576923691e-2);

/* Radians to arcseconds */
/**
 * @brief 弧度转秒
 */
const double DR2AS(206264.8062470963551564734);

/* Arcseconds to radians */
/**
 * @brief 秒到弧度的转换
 */
const double DAS2R(4.848136811095359935899141e-6);

/* Seconds of time to radians */
/**
 * @brief 每秒地球的旋转弧度
 */
const double DS2R(7.272205216643039903848712e-5);

/* Arcseconds in a full circle */
const double TURNAS(1296000.0);

/* Milliarcseconds to radians */
const double DMAS2R(DAS2R / 1e3);

/* Length of tropical year B1900 (days) */
const double DTY(365.242198781);

/* Seconds per day. */
const double DAYSEC(86400.0);
const double SECDAY(1.0/DAYSEC);
/* Days per Julian year */
const double DJY(365.25);

/* Days per Julian century */
const double DJC(36525.0);

/* Days per Julian millennium */
const double DJM(365250.0);

/* Reference epoch (J2000.0), Julian Date */
const double DJ00(2451545.0);

/* Julian Date of Modified Julian Date zero */
const double DJM0(2400000.5);

/* Reference epoch (J2000.0), Modified Julian Date */
const double DJM00(51544.5);

/* 1977 Jan 1.0 as MJD */
const double DJM77(43144.0);

/* TT minus TAI (s) */
const double TTMTAI(32.184);

/* Astronomical unit (m, IAU 2012)*/
const double DAU(149597870.7e3);

/* Speed of light (m/s) */
const double CMPS(299792458.0);

/* Light time for 1 au (s) */
const double AULT(DAU/CMPS);

/* Speed of light (AU per day) */
const double DC(DAYSEC / AULT);

/* L_G = 1 - d(TT)/d(TCG) */
const double ELG(6.969290134e-10);

/* L_B = 1 - d(TDB)/d(TCB), and TDB (s) at TAI 1977/1/1.0 */
const double ELB(1.550519768e-8);
const double TDB0(-6.55e-5);

/* Schwarzschild radius of the Sun (au) */
/* = 2 * 1.32712440041e20 / (2.99792458e8)^2 / 1.49597870700e11 */
const double SRS(1.97412574336e-8);

/* Speed of light  [m/s]; IAU 1976 */
const double DLIGHT   = 299792458.0;

// Gravitational coefficient
const double GM_Earth    = 398600.4415e+9;    // [m^3/s^2]; JGM3
//const double GM_Sun      = 1.32712438e+20;    // [m^3/s^2]; IAU 1976
//const double GM_Moon     = GM_Earth/81.300587;// [m^3/s^2]; DE200
const double GM_Sun      = 1.327122e+20;    // [m^3/s^2]; STK
const double GM_Moon     = 4.902801076e+12;// [m^3/s^2]; STK

const double DF_ZERO = 1e-12;
const double R_Earth     = 6378245.0000000000/*6378.137e3*/;// Radius Earth [m]; WGS-84
const double R_Earth2    = 6356755.2881600000;

const double R_Sun = 695990000.;// Radius Sun [m];
const double R_Moon = 1738000.; // Radius Moon [m];
// Solar radiation pressure at 1 AU
const double P_Sol       = 4.560E-6;          // [N/m^2] (~1367 W/m^2); IERS 96

/* Reference ellipsoids */
const int WGS84(1);
const int GRS80(2);
const int WGS72(3);
const int BJ54(4);
const int CGCS2000(5);
#endif
