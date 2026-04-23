#include "functions.h"
#include <cmath>

namespace ariadna {

void POLE_TIDE(double cent, double lat_geod, double lon_geod,
               double xp, double yp, double xp_rate, double yp_rate,
               const Eigen::Matrix3d& vw_i, const Eigen::MatrixXd& r2000,
               Eigen::Vector3d& dx_poltide, Eigen::Vector3d& dv_poltide) {
    
    // IERS 2010 Constants for Pole Tide
    const double factor_r = -0.033; // Radial factor [m/arcsec]
    const double factor_t = 0.009;  // Transverse factor [m/arcsec]

    // Mean pole path (IERS 2010 linear trend approximation)
    double t_year = 2000.0 + cent * cnst::JUL_CENT / 365.25; 
    double dt = t_year - 2000.0;
    
    double xm = 0.054 + 0.00083 * dt; // arcsec
    double ym = 0.357 + 0.00395 * dt; // arcsec

    double m1 = xp - xm;
    double m2 = -(yp - ym); 

    double cos_lat = std::cos(lat_geod);
    double sin_lat = std::sin(lat_geod);
    double cos_lon = std::cos(lon_geod);
    double sin_lon = std::sin(lon_geod);
    double cos_2lat = std::cos(2.0 * lat_geod);

    // Displacements in Local Topocentric (VEN: Up, East, North)
    double u_r = factor_r * sin_lat * cos_lat * (m1 * cos_lon + m2 * sin_lon);
    double u_e = factor_t * sin_lat * (m1 * sin_lon - m2 * cos_lon);
    double u_n = -factor_t * cos_2lat * (m1 * cos_lon + m2 * sin_lon);

    Eigen::Vector3d dr_ven(u_r, u_e, u_n);

    // Velocities
    double dm1_dt = xp_rate;
    double dm2_dt = -yp_rate;

    double v_r = factor_r * sin_lat * cos_lat * (dm1_dt * cos_lon + dm2_dt * sin_lon);
    double v_e = factor_t * sin_lat * (dm1_dt * sin_lon - dm2_dt * cos_lon);
    double v_n = -factor_t * cos_2lat * (dm1_dt * cos_lon + dm2_dt * sin_lon);

    Eigen::Vector3d dv_ven(v_r, v_e, v_n);

    // Convert VEN to ITRF using vw_i
    Eigen::Vector3d dx_itrf = vw_i * dr_ven;
    Eigen::Vector3d dv_itrf = vw_i * dv_ven;

    // Convert ITRF to J2000 using r2000
    dx_poltide = r2000.block<3, 3>(0, 0) * dx_itrf;
    dv_poltide = r2000.block<3, 3>(0, 3) * dx_itrf + r2000.block<3, 3>(0, 0) * dv_itrf;
}
}