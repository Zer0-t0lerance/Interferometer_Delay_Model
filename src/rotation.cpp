#include "functions.h"
#include <cmath>

// Предполагается, что твои константы лежат здесь
#include "constants.h" 

namespace ariadna {

// Фундаментальные константы
constexpr double PI = 3.14159265358979323846;
constexpr double SEC_TO_RAD = PI / (180.0 * 3600.0);
constexpr double OMEGA_EARTH = 7.292115146706979e-5; // Угловая скорость Земли, рад/с

// ============================================================================
// 1. ПРЕЦЕССИЯ (P)
// ============================================================================
Eigen::Matrix3d calc_precession(double TDB) {
    double T2 = TDB * TDB;
    double T3 = TDB * T2;

    // Углы прецессии (из градусов в радианы)
    double zeta  = (0.011180860 * TDB + 1.464e-6 * T2 + 8.7e-8 * T3) * PI / 180.0;
    double z     = (0.011180860 * TDB + 5.307e-6 * T2 + 8.8e-8 * T3) * PI / 180.0;
    double theta = (0.009717173 * TDB - 2.068e-6 * T2 - 2.0e-8 * T3) * PI / 180.0;

    // Rz(-z) * Ry(theta) * Rz(-zeta)
    return Eigen::AngleAxisd(-z, Eigen::Vector3d::UnitZ()) *
           Eigen::AngleAxisd(theta, Eigen::Vector3d::UnitY()) *
           Eigen::AngleAxisd(-zeta, Eigen::Vector3d::UnitZ()).toRotationMatrix();
}

// ============================================================================
// 2. НУТАЦИЯ (N)
// ============================================================================
Eigen::Matrix3d calc_nutation(double TDB, double& dpsi, double& deps, double& eps_true) {
    double T2 = TDB * TDB;
    double T3 = TDB * T2;

    // Фундаментальные аргументы (в радианах)
    double l  = (134.96298139 + 477198.8673981 * TDB + 0.0086972 * T2 + 1.78e-5 * T3) * PI / 180.0;
    double lp = (357.52772333 +  35999.05034   * TDB - 0.0001603 * T2 - 3.3e-6  * T3) * PI / 180.0;
    double F  = ( 93.27191028 + 483202.0175381 * TDB - 0.0036825 * T2 + 3.1e-6  * T3) * PI / 180.0;
    double D  = (297.85036306 + 445267.11148   * TDB - 0.0019142 * T2 + 5.3e-6  * T3) * PI / 180.0;
    double Om = (125.04452222 -   1934.1362611 * TDB + 0.0020708 * T2 + 2.2e-6  * T3) * PI / 180.0;

    dpsi = 0.0;
    deps = 0.0;

    // Модульный цикл по 106 коэффициентам из твоего constants.cpp
    for (int i = 0; i < 106; ++i) {
        double arg = cnst::coef_a1[i] * l + cnst::coef_a2[i] * lp + cnst::coef_a3[i] * F + cnst::coef_a4[i] * D + cnst::coef_a5[i] * Om;
        dpsi += (cnst::coef_S[i] + cnst::coef_dS[i] * TDB) * std::sin(arg);
        deps += (cnst::coef_C[i] + cnst::coef_dC[i] * TDB) * std::cos(arg);
    }

    // Перевод в радианы
    dpsi *= 1e-4 * SEC_TO_RAD;
    deps *= 1e-4 * SEC_TO_RAD;

    // Истинный наклон эклиптики
    double eps_mean = (23.43929111 - 0.013004167 * TDB - 1.6389e-7 * T2 + 5.036e-7 * T3) * PI / 180.0;
    eps_true = eps_mean + deps;

    // Rx(-eps_true) * Rz(-dpsi) * Rx(eps_mean)
    return Eigen::AngleAxisd(-eps_true, Eigen::Vector3d::UnitX()) *
           Eigen::AngleAxisd(-dpsi, Eigen::Vector3d::UnitZ()) *
           Eigen::AngleAxisd(eps_mean, Eigen::Vector3d::UnitX()).toRotationMatrix();
}

// ============================================================================
// 3. ВРАЩЕНИЕ ЗЕМЛИ (S) И ЕГО ПРОИЗВОДНАЯ
// ============================================================================
void calc_earth_rotation(double JD_UT1, double dpsi, double eps_true, 
                         Eigen::Matrix3d& S, Eigen::Matrix3d& dS_dt) {
    double dUT1 = JD_UT1 - 2451545.0;
    double T_UT1 = dUT1 / 36525.0;

    // Гринвичское среднее звездное время (GMST) в секундах
    double gmst_sec = 24110.54841 + 8640184.812866 * T_UT1 + 0.093104 * T_UT1 * T_UT1 - 6.2e-6 * std::pow(T_UT1, 3);
    
    // Перевод в радианы. Добавляем полный оборот Земли (dUT1 * 2 * PI), чтобы Земля вращалась каждый день!
    double gmst_rad = std::fmod(gmst_sec * 15.0 * SEC_TO_RAD + dUT1 * 2.0 * PI, 2.0 * PI);
    
    // Уравнение равноденствий -> Истинное время (GAST)
    double gast_rad = gmst_rad + dpsi * std::cos(eps_true);
    if (gast_rad < 0) gast_rad += 2.0 * PI;

    // Матрица вращения
    S = Eigen::AngleAxisd(-gast_rad, Eigen::Vector3d::UnitZ()).toRotationMatrix();

    // Производная матрицы вращения по времени (матрица угловой скорости)
    // dRz/dt = Rz(-GAST) * [0, w, 0; -w, 0, 0; 0, 0, 0]
    Eigen::Matrix3d Omega;
    Omega <<  0.0,         OMEGA_EARTH, 0.0,
             -OMEGA_EARTH, 0.0,         0.0,
              0.0,         0.0,         0.0;
    
    dS_dt = S * Omega;
}

// ============================================================================
// 4. ДВИЖЕНИЕ ПОЛЮСА (M)
// ============================================================================
Eigen::Matrix3d calc_polar_motion(double xp, double yp) {
    // Вращение Ry(-xp) * Rx(-yp)
    return Eigen::AngleAxisd(-xp, Eigen::Vector3d::UnitY()) *
           Eigen::AngleAxisd(-yp, Eigen::Vector3d::UnitX()).toRotationMatrix();
}

// ============================================================================
// ИТОГОВАЯ ФУНКЦИЯ: ITRF -> J2000 (Координаты и Скорости)
// ============================================================================
void get_r2000_matrices(double JD_TDB, double JD_UT1, double xp, double yp, 
                        Eigen::Matrix3d& R2000, Eigen::Matrix3d& dR2000_dt) {
    
    double TDB = (JD_TDB - 2451545.0) / 36525.0;

    Eigen::Matrix3d P = calc_precession(TDB);
    
    double dpsi, deps, eps_true;
    Eigen::Matrix3d N = calc_nutation(TDB, dpsi, deps, eps_true);

    Eigen::Matrix3d S, dS_dt;
    calc_earth_rotation(JD_UT1, dpsi, eps_true, S, dS_dt);

    Eigen::Matrix3d M = calc_polar_motion(xp, yp);

    // Матрица поворота для координат: X_j2000 = R2000 * X_itrf
    R2000 = P * N * S * M;

    // Матрица поворота для скоростей: V_j2000 = dR2000_dt * X_itrf
    dR2000_dt = P * N * dS_dt * M;
}

} // namespace ariadna