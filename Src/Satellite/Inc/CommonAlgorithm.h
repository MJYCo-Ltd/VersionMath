#ifndef YTY_COMMONALGORITHM_H
#define YTY_COMMONALGORITHM_H

/*****************************************
  作用：各个需要的基础算法，通用算法
  注意：此类是非线程安全的，在多线程中使用需要对变量
       进行加锁保护，否则运算结果不能保证
       对外没有接口，让外部使用，只能内部使用
 *****************************************/

namespace Math{
class CVector;
class CMatrix;
}

/// 用于存放 tle数据的结构体
typedef struct
{
  double epoch, xndt2o, xndd6o, bstar;
  double xincl, xnodeo, eo, omegao, xmo, xno;
  int norad_number, bulletin_number, revolution_number;
  char classification;    /* "U" = unclassified;  only type I've seen */
  char ephemeris_type;
  char intl_desig[9];
} tle_t;

////////////天文学基础///////////


/**
 * @brief 计算岁差的各项值
 * @param epoch_from
 * @param epoch_to
 * @param eta
 * @param pi
 * @param p
 */
void compute_ecliptic_precession_angles( const double epoch_from,
       const double epoch_to, double *eta, double *pi, double *p);

/**
 * @brief J2000参考系下岁差算法
 * @param Mjd_1 历元时间        [MJD UTT]
 * @param Mjd_2 要计算岁差的时间 [MJD UTT]
 * @return 岁差旋转矩阵
 */
void PrecMatrix (double Mjd_2, double dTemp[3][3]);

/**
 * @brief J2000参考系下章动矩阵
 * @param Mjd_TT 要计算章动的时间 [MJD UTT]
 * @return 章动旋转矩阵
 */
void NutMatrix (double Mjd_TT, double &dpsi, double &eps, double &deps, double &dom, double dTemp[3][3]);

/**
 * @brief EquationOfEquinoxes
 * @param dpsi
 * @param eps
 * @param om
 * @return
 */
double EquationOfEquinoxes(double dpsi, double deps, double eps, double om);

/**
 * @brief 格林尼治真恒星时
 * @param Mjd_UT1 要计算的时间         [MJD UT1]
 * @return 从平春分点到格林尼治经线的角度 [rad]
 */
double GAST (double Mjd_UT1, double dpsi, double deps, double eps, double om);

/**
 * @brief J2000参考系下地球旋转矩阵
 * @param Mjd_UT1 要计算的时间  [MJD UT1]
 * @return 地球旋转矩阵
 */
void GHAMatrix (double Mjd_UT1, double dpsi, double deps, double eps, double om, double dTemp[3][3]);

/**
 * @brief 极移矩阵
 * @param dX [rad]
 * @param dY [rad]
 * @return 极移矩阵
 */
void PoleMatrix(const double& dX, const double& dY, double dTemp[3][3]);
////////////天文学基础 end///////////

////////////轨道处理////////////////

/**
 * @brief 解广义开普勒方程
 * @param ksi   e*cos(w)
 * @param eta   e*sin(w)
 * @param lamda w + M   [rad]
 * @return E+w          [rad]
 * @attention e偏心率 w近地点幅角 E偏近点角 M平近点角
 */
double KeplerFunc2(double ksi,double eta,double lamda);

/**
 * @brief 计算短周期项
 * @param MElem 平均根数
 * @param ZS    短周期项(a_s,i_s,Omega_s,ksi_s,eta_s,lamda_s)
 * @param AS2   a的二阶长期项
 */
void Short2(const Math::CVector & MElem,double ZS[6],double & AS2);

/**
 * @brief 被FindEta调用
 * @return
 */
double F(double eta, double m, double l);

/**
 * @brief FindEta
 * @param r_a 在a时刻的位置 (x,y,z)[m]
 * @param r_b 在并时刻的位置 (x,y,z)[m]
 * @param tau
 * @return
 */
double FindEta(const Math::CVector& r_a, const Math::CVector& r_b, double tau);

/**
 * @brief 被 CKepler::ClassicalElements 调用
 */
double remaining_terms( const double ival);

/**
 * @brief 将Kepler的六根数转换成TLE的
 * @param epoch_from 历元时间 2000. [JD/365.25]
 * @param epoch_to   要计算的历元 [JD/365.25]
 * @param incl       轨道倾角   [rad]
 * @param asc_node   升交点赤经 [rad]
 * @param arg_per    近地点幅角 [rad]
 */
void convert_elements( const double epoch_from, const double epoch_to,
                       double *incl, double *asc_node, double *arg_per);

/**
 * @brief 给TLE添加校验值
 * @param buff
 * @return 失败返回 0 正确返回 1
 * @attention 失败的话，则认为buff对应的不是TLE
 */
int add_tle_checksum_data( char *buff);

/**
 * @brief 将值放到 TLE字符串中
 */
void put_sci( char *obuff, double ival);

/**
 * @brief 将TLE数据写入到字符串中
 */
void write_elements_in_tle_format( char buff[2][73], const tle_t *tle);
////////////轨道处理end////////////////


////////////大气////////////////////////

/**
 * @brief Harris-Priester大气密度模型
 * @param vecRSun  太阳的ECI位置 (x,y,z) [m]
 * @param r_tod    卫星ECF位置  (x,y,z) [m]
 * @return 大气密度 [kg/m^3]
 */
double Density_HP(const Math::CVector &vecRSun, const Math::CVector& r_tod );

/**
 * @brief 美国1976标准大气模型
 * @param 距地表高度   [m]
 * @return 大气密度 [kg/m^3]
 */
double SA76(double dH);
////////////大气 end ///////////////////
/**
 * @brief CalP
 * @param dT
 * @return
 * @attention CTimeSys::GetTCB 中调用
 */
double CalP(double dT);
#endif // YTY_COMMONALGORITHM_H

