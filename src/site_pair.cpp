// site_pair.cpp
//
// Посуточная сборка координат ПАРЫ станций в системе J2000.0 (роль SUBROUTINE SITE
// в основном цикле ARIADNA): для каждой из двух станций вычисляет геофизические
// поправки (твёрдые приливы, океаническая нагрузка, прилив полюса, атмосферная
// нагрузка, термодеформация антенны) и через SITE_INST переводит в J2000 с
// производными.
//
// Конвенции подтверждены сверкой с дампами:
//   * Солнце/Луна ДЛЯ ПРИЛИВОВ — ГЕОЦЕНТРИЧЕСКИЕ (jpleph), 3x2 (поз/скор) [м, м/с];
//   * r2000 (ITRF->J2000) и его производные dr2000, d2r2000.

#include "functions.h"

namespace ariadna {

static void site_one(const SitePrep& s,
                     int mjd, double utc, double jd, double ut1_sec, double cent,
                     const Eigen::VectorXd& f, const Eigen::VectorXd& fd,
                     const Eigen::Vector2d& gast,
                     const Eigen::Matrix<double, 3, 2>& sun_geo,
                     const Eigen::Matrix<double, 3, 2>& moon_geo,
                     double xp, double yp, double xp_rate, double yp_rate,
                     const Eigen::Matrix3d& r2000, const Eigen::Matrix3d& dr2000,
                     const Eigen::Matrix3d& d2r2000,
                     Eigen::Vector3d& x_j2000, Eigen::Vector3d& v_j2000, Eigen::Vector3d& a_j2000) {
    // Твёрдые земные приливы: r2000 передаётся блоком 9x3 (r0=R, r1=dR).
    Eigen::MatrixXd r2000_blk = Eigen::MatrixXd::Zero(9, 3);
    r2000_blk.block<3, 3>(0, 0) = r2000;
    r2000_blk.block<3, 3>(3, 0) = dr2000;
    Eigen::Vector3d dx_tide, dv_tide;
    SITE_TIDE_SOLID(s.xsta_itrf, s.lat_gcen, s.lon_gcen, sun_geo, moon_geo,
                    f, fd, s.vw_i, gast, r2000_blk, dx_tide, dv_tide);

    // Океаническая нагрузка.
    Eigen::Vector3d dx_oc, dv_oc;
    SITE_TIDE_OC(jd, ut1_sec, s.tide_data, s.vw_i, r2000, dr2000, dx_oc, dv_oc);

    // Прилив полюса.
    Eigen::Vector3d dx_pol, dv_pol;
    POLE_TIDE(cent, s.lat_geod, s.lon_gcen, xp, yp, xp_rate, yp_rate,
              s.vw_i, r2000, dr2000, dx_pol, dv_pol);

    // Атмосферная нагрузка.
    Eigen::Vector3d dx_atm, dv_atm;
    SITE_ATM40(mjd, utc, s.pres, s.dPdt, s.atm_load, s.vw_i, r2000, dr2000, dx_atm, dv_atm);

    // Термодеформация антенны (Nothnagel 2009): расширение фундамента+колонны.
    Eigen::Vector3d dx_temp, dv_temp;
    THERM_DEF40(s.tC, s.dTdt, s.def_par, s.vw_i, r2000, dr2000, dx_temp, dv_temp);

    // Суммирование и перевод в J2000 с производными.
    SITE_INST(s.xsta_itrf, r2000, dr2000, d2r2000,
              dx_tide, dv_tide, dx_oc, dv_oc, dx_pol, dv_pol, dx_atm, dv_atm, dx_temp, dv_temp,
              x_j2000, v_j2000, a_j2000);
}

void site_pair(const SitePrep& s1, const SitePrep& s2,
               int mjd, double utc, double jd, double ut1_sec, double cent,
               const Eigen::VectorXd& f, const Eigen::VectorXd& fd,
               const Eigen::Vector2d& gast,
               const Eigen::Matrix<double, 3, 2>& sun_geo,
               const Eigen::Matrix<double, 3, 2>& moon_geo,
               double xp, double yp, double xp_rate, double yp_rate,
               const Eigen::Matrix3d& r2000, const Eigen::Matrix3d& dr2000,
               const Eigen::Matrix3d& d2r2000,
               std::vector<Eigen::Vector3d>& xsta_j2000,
               std::vector<Eigen::Vector3d>& vsta_j2000,
               std::vector<Eigen::Vector3d>& asta_j2000) {
    xsta_j2000.resize(2); vsta_j2000.resize(2); asta_j2000.resize(2);
    site_one(s1, mjd, utc, jd, ut1_sec, cent, f, fd, gast, sun_geo, moon_geo,
             xp, yp, xp_rate, yp_rate, r2000, dr2000, d2r2000,
             xsta_j2000[0], vsta_j2000[0], asta_j2000[0]);
    site_one(s2, mjd, utc, jd, ut1_sec, cent, f, fd, gast, sun_geo, moon_geo,
             xp, yp, xp_rate, yp_rate, r2000, dr2000, d2r2000,
             xsta_j2000[1], vsta_j2000[1], asta_j2000[1]);
}

} // namespace ariadna
