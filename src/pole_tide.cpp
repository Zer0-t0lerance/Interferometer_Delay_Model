#include "functions.h"
#include <cmath>

namespace ariadna {

void POLE_TIDE(double cent, double lat_geod, double lon_geod,
               double xp, double yp, double xp_rate, double yp_rate,
               const Eigen::Matrix3d& vw_i, const Eigen::Matrix3d& r2000,
               const Eigen::Matrix3d& dr2000_dt,
               Eigen::Vector3d& dx_poltide, Eigen::Vector3d& dv_poltide) 
{
    // IERS 2000 Constants for Pole Tide (строго по Фортрану)
    const double factor_r = -0.032; // Radial factor [m/arcsec] (-32.D0 * 1.D-3)
    const double factor_t = 0.009;  // Transverse factor [m/arcsec] (9.D0 * 1.D-3)

    // Вековой дрейф полюса (t_0 = 2000)
    const double x_mean = 0.054;
    const double y_mean = 0.357;
    const double xdot_mean = 0.00083; // arcsec/year
    const double ydot_mean = 0.00395; // arcsec/year

    double dt_years = cent * 100.0;
    
    // Среднее положение полюса
    double x_meant = x_mean + xdot_mean * dt_years;
    double y_meant = y_mean + ydot_mean * dt_years;

    // Смещение полюса (arcsec)
    double m1 = xp - x_meant;
    double m2 = -(yp - y_meant); 

    // Скорости изменения смещения полюса (arcsec/sec)
    double dx_mean_dt = xdot_mean * 100.0 / cnst::JUL_CENT / cnst::SECDAY;
    double dy_mean_dt = ydot_mean * 100.0 / cnst::JUL_CENT / cnst::SECDAY;
    double dm1_dt = xp_rate - dx_mean_dt;
    double dm2_dt = -(yp_rate - dy_mean_dt);

    // Геометрия станции
    double colat = cnst::HALFPI - lat_geod; // Коширота
    double cos_lon = std::cos(lon_geod);
    double sin_lon = std::sin(lon_geod);
    double cos_colat = std::cos(colat);
    double cos_2colat = std::cos(2.0 * colat);
    double sin_2colat = std::sin(2.0 * colat);

    // --- Топоцентрические смещения (Up, East, North) ---
    // В Фортране drse(3) считается для Юга, а затем инвертируется. Мы считаем сразу Север.
    double drse_up    = factor_r * sin_2colat * (m1 * cos_lon + m2 * sin_lon);
    double drse_east  = factor_t * cos_colat  * (m1 * sin_lon - m2 * cos_lon);
    double drse_south = -factor_t * cos_2colat * (m1 * cos_lon + m2 * sin_lon);
    
    Eigen::Vector3d dr_ven(drse_up, drse_east, -drse_south);

    // --- Топоцентрические скорости (Up, East, North) ---
    double drsev_up    = factor_r * sin_2colat * (dm1_dt * cos_lon + dm2_dt * sin_lon);
    double drsev_east  = factor_t * cos_colat  * (dm1_dt * sin_lon - dm2_dt * cos_lon);
    double drsev_south = -factor_t * cos_2colat * (dm1_dt * cos_lon + dm2_dt * sin_lon);

    Eigen::Vector3d dv_ven(drsev_up, drsev_east, -drsev_south);

    // --- Преобразование VEN -> ITRF ---
    Eigen::Vector3d dx_itrf = vw_i * dr_ven;
    Eigen::Vector3d dv_itrf = vw_i * dv_ven;

    // --- Преобразование ITRF -> J2000 ---
    dx_poltide = r2000 * dx_itrf;
    dv_poltide = r2000 * dv_itrf + dr2000_dt * dx_itrf;
}

} // namespace ariadna