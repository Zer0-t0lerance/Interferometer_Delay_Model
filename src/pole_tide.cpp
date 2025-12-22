// POLE_TIDE.cpp

#include "functions.h"

namespace ariadna {

void POLE_TIDE(double cent, double lat_geod, double lon_geod,
               double xp, double yp, double xp_rate, double yp_rate,
               const Eigen::Matrix3d& vw_i, const Eigen::MatrixXd& r2000,
               Eigen::Vector3d& dx_poltide, Eigen::Vector3d& dv_poltide) {
    
    using namespace cnst;
    using Eigen::Vector3d;
    using Eigen::Matrix3d;
    using std::sin;
    using std::cos;
    using std::pow;

    // Константы
    const double G3 = G3_LOVE_NUMBER;
    const double C_rad = CARCRAD; // arcsec to rad

    // 1. Расчет разницы полюсов: (фактический - средний)
    
    // Эпоха в годах с 2000.0
    double t_year = 2000.0 + cent * JUL_CENT / 365.25; 
    double dt = t_year - T0_YEAR_2000; // t - t0 [years]

    // Среднее положение полюса по линейному тренду [arcsec]
    double xm_2 = X_MEAN_T0_ARCSEC + dt * X_DOT_MEAN_T0_ARCSEC_PER_YEAR;
    double ym_2 = Y_MEAN_T0_ARCSEC + dt * Y_DOT_MEAN_T0_ARCSEC_PER_YEAR;

    // Вариация полюса относительно среднего тренда [rad]
    double dm_2_x = (xp - xm_2 * C_rad);
    double dm_2_y = (yp - ym_2 * C_rad);
    
    // Скорость вариации полюса относительно среднего тренда [rad/day]
    // Примечание: d(xm/ym)/dt = 0, так как тренд - константа по времени, но 
    // Fortran-код использует dxdt/dydt (скорость полюса), поэтому dm_2_v - это просто dxdt/dydt.
    // Переводим xp_rate и yp_rate из [rad/day] в [rad/s] для dv
    double dm_2_vx = xp_rate;
    double dm_2_vy = yp_rate; 

    // 2. Расчет топоцентрических смещений (drse) и скоростей (drse_v)
    
    // Формулы (IERS Conventions 2000, Eq 6.1, 6.2)
    // S_r = -32 * G3 * cos(2*lat) * (xm*cos(lon) + ym*sin(lon)) [mm] -> [m]
    // S_phi = 32 * G3 * sin(2*lat) * (xm*cos(lon) + ym*sin(lon)) [mm] -> [m]
    // S_lambda = 32 * G3 * sin(lat) * (xm*sin(lon) - ym*cos(lon)) [mm] -> [m]
    
    // 32.0 * G3 * 1e-3 = 32.0 * 0.0036 * 0.001 = 1.152e-4
    const double FACTOR_D = 32.0 * G3 * 0.001; 
    
    double cos_lat = cos(lat_geod);
    double sin_lat = sin(lat_geod);
    double cos_2lat = cos(2.0 * lat_geod);
    double sin_2lat = sin(2.0 * lat_geod);
    double cos_lon = cos(lon_geod);
    double sin_lon = sin(lon_geod);

    // Вспомогательные термины для смещения
    double M_d = dm_2_x * cos_lon + dm_2_y * sin_lon; // (x-xm)cos(lon) + (y-ym)sin(lon)
    double N_d = dm_2_x * sin_lon - dm_2_y * cos_lon; // (x-xm)sin(lon) - (y-ym)cos(lon)

    // Топоцентрические смещения (Up, North, East) [m]
    Vector3d drse;
    drse(0) = -FACTOR_D * cos_2lat * M_d;  // Up
    drse(1) = FACTOR_D * sin_2lat * M_d;   // North
    drse(2) = FACTOR_D * sin_lat * N_d;    // East
    
    // Вспомогательные термины для скорости
    // Скорости полюса dX/dt, dY/dt
    double M_dv = dm_2_vx * cos_lon + dm_2_vy * sin_lon;
    double N_dv = dm_2_vx * sin_lon - dm_2_vy * cos_lon;
    
    // Топоцентрические скорости (Up, North, East) [m/s]
    Vector3d drse_v;
    drse_v(0) = -FACTOR_D * cos_2lat * M_dv / SECDAY; // /SECDAY для перевода [rad/day] -> [rad/s]
    drse_v(1) = FACTOR_D * sin_2lat * M_dv / SECDAY;
    drse_v(2) = FACTOR_D * sin_lat * N_dv / SECDAY;

    // 3. Преобразование из топоцентрической (Up, North, East) в J2000.0

    // VW_i - матрица перехода из VEN (Up, North, East) в ITRF (краст-фиксированную X,Y,Z).
    // dr_itrf = VW_i * drse
    Vector3d dx_itrf = vw_i * drse;
    Vector3d dv_itrf = vw_i * drse_v;

    // R2000 - матрица перехода из ITRF в J2000.0
    // dx_j2000 = R * dx_itrf
    const Matrix3d R = r2000.block<3, 3>(0, 0); // R2000(3,3,0)
    dx_poltide = R * dx_itrf;

    // dv_j2000 = R_dot * dx_itrf + R * dv_itrf
    const Matrix3d R_dot = r2000.block<3, 3>(0, 3); // R2000(3,3,1)
    dv_poltide = R_dot * dx_itrf + R * dv_itrf;
}

} // namespace ariadna