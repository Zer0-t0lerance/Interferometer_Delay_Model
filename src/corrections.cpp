#include "functions.h"
#include "structures.h"
#include "constants.h"

// Заглушки для функций коррекций

void site_tide_solid(const Eigen::Vector3d& xsta, double lat_gcen, double lon_gcen,
                     const Eigen::Vector3d& sun, const Eigen::Vector3d& moon,
                     const Eigen::VectorXd& f, const Eigen::VectorXd& fd,
                     const Eigen::Matrix3d& vw_i, double gast,
                     const Eigen::MatrixXd& r2000,
                     Eigen::Vector3d& dxtide, Eigen::Vector3d& dvtide) {
    // Реализация на основе SITE_TIDE_SOLID.for.txt
    Eigen::Vector3d xsun = sun;
    Eigen::Vector3d xmon = moon;
    Eigen::Vector3d vsun = Eigen::Vector3d::Zero(); // Заглушка для скорости
    Eigen::Vector3d vmon = Eigen::Vector3d::Zero();

    Eigen::Vector3d xsta_j2000 = r2000.block<3,3>(0,0) * xsta;
    Eigen::Vector3d vsta_j2000 = r2000.block<3,3>(0,3) * xsta;

    double rsta = xsta_j2000.norm();
    double rsun = xsun.norm();
    double rmon = xmon.norm();

    Eigen::Vector3d xsta_n = xsta_j2000 / rsta;
    Eigen::Vector3d xsun_n = xsun / rsun;
    Eigen::Vector3d xmon_n = xmon / rmon;

    double SCSUN = xsta_n.dot(xsun_n);
    double SCMON = xsta_n.dot(xmon_n);

    double w1 = 1.5 * sin(lat_gcen) * sin(lat_gcen) - 0.5;
    double H2 = cnst::H02 + cnst::H22 * w1;
    double L2 = cnst::L02 + cnst::L22 * w1;

    double H2_2 = H2 / 2.0;
    double P2SUN = 3.0 * (H2_2 - L2) * SCSUN * SCSUN - H2_2;
    double P2MON = 3.0 * (H2_2 - L2) * SCMON * SCMON - H2_2;

    double P3SUN = 2.5 * (cnst::H3 - 3.0 * cnst::L3) * SCSUN * SCSUN * SCSUN + 1.5 * (cnst::L3 - cnst::H3) * SCSUN;
    double P3MON = 2.5 * (cnst::H3 - 3.0 * cnst::L3) * SCMON * SCMON * SCMON + 1.5 * (cnst::L3 - cnst::H3) * SCMON;

    double X2SUN = 3.0 * L2 * SCSUN;
    double X2MON = 3.0 * L2 * SCMON;
    double X3SUN = 1.5 * cnst::L3 * (5.0 * SCSUN * SCSUN - 1.0);
    double X3MON = 1.5 * cnst::L3 * (5.0 * SCMON * SCMON - 1.0);

    double RE = cnst::AE;
    double FGM11 = RE * pow(RE / rsun, 3) * cnst::GSUN / cnst::GEARTH;
    double FGM12 = RE * pow(RE / rmon, 3) * cnst::GMOON / cnst::GEARTH;
    double FGM21 = RE / rsun * FGM11;
    double FGM22 = RE / rmon * FGM12;

    dxtide(0) = FGM11 * (X2SUN * xsun_n(0) + P2SUN * xsta_n(0)) +
                FGM12 * (X2MON * xmon_n(0) + P2MON * xsta_n(0)) +
                FGM21 * (X3SUN * xsun_n(0) + P3SUN * xsta_n(0)) +
                FGM22 * (X3MON * xmon_n(0) + P3MON * xsta_n(0));
    dxtide(1) = FGM11 * (X2SUN * xsun_n(1) + P2SUN * xsta_n(1)) +
                FGM12 * (X2MON * xmon_n(1) + P2MON * xsta_n(1)) +
                FGM21 * (X3SUN * xsun_n(1) + P3SUN * xsta_n(1)) +
                FGM22 * (X3MON * xmon_n(1) + P3MON * xsta_n(1));
    dxtide(2) = FGM11 * (X2SUN * xsun_n(2) + P2SUN * xsta_n(2)) +
                FGM12 * (X2MON * xmon_n(2) + P2MON * xsta_n(2)) +
                FGM21 * (X3SUN * xsun_n(2) + P3SUN * xsta_n(2)) +
                FGM22 * (X3MON * xmon_n(2) + P3MON * xsta_n(2));

    // Аналогично для dvtide, но упрощено
    dvtide.setZero();
}

void site_tide_oc(int mjd_start, double ut1_sec, const ariadna::OceanTideData& tide_data,
                  const Eigen::Matrix3d& vw_i, const Eigen::MatrixXd& r2000,
                  Eigen::Vector3d& dx_octide, Eigen::Vector3d& dv_octide) {
    dx_octide.setZero();
    dv_octide.setZero();
    // Реализовать океанические приливы
}

void pole_tide(double cent, double lat_geod, double lon_geod,
               double xp, double yp, double xp_rate, double yp_rate,
               const Eigen::Matrix3d& vw_i, const Eigen::MatrixXd& r2000,
               Eigen::Vector3d& dx_poltide, Eigen::Vector3d& dv_poltide) {
    dx_poltide.setZero();
    dv_poltide.setZero();
    // Реализовать прилив полюса
}

void site_atm40(double dPdt, const Eigen::Matrix3d& vw_i, const Eigen::MatrixXd& r2000,
                Eigen::Vector3d& dx_atm, Eigen::Vector3d& dv_atm) {
    dx_atm.setZero();
    dv_atm.setZero();
    // Реализовать атмосферную нагрузку
}

void dmeteo1_dt(const std::vector<Observation>& observations, const std::vector<Station>& stations, double t_mean, std::vector<Eigen::Vector3d>& site_meteo, int& ndeg, std::vector<Eigen::VectorXd>& t_coef, std::vector<Eigen::VectorXd>& p_coef, std::vector<Eigen::VectorXd>& hum_coef) {
    // Заглушка: вычисление коэффициентов полинома для метео параметров
    ndeg = 2; // Степень полинома
    int n_stations = stations.size();
    t_coef.resize(n_stations, Eigen::VectorXd::Zero(ndeg + 1));
    p_coef.resize(n_stations, Eigen::VectorXd::Zero(ndeg + 1));
    hum_coef.resize(n_stations, Eigen::VectorXd::Zero(ndeg + 1));
    site_meteo.resize(n_stations, Eigen::Vector3d::Zero());
    // В реальности нужно аппроксимировать данные из observations
}

void dmeteo2_dt(const Station& station, int n_stations, int ndeg, const Observation& obs, double t_mean, const std::vector<Eigen::VectorXd>& t_coef, const std::vector<Eigen::VectorXd>& p_coef, const std::vector<Eigen::VectorXd>& hum_coef, double& dtdt, double& dpdt, double& dhumdt) {
    // Заглушка: вычисление производных метео параметров
    dtdt = 0.0;
    dpdt = 0.0;
    dhumdt = 0.0;
    // В реальности вычислить из полинома
}

void site_inst(const Eigen::Vector3d& xsta_itrf, const Eigen::Vector3d& vsta_itrf,
               const Eigen::MatrixXd& r2000,
               const Eigen::Vector3d& dx_tide, const Eigen::Vector3d& dv_tide,
               const Eigen::Vector3d& dx_octide, const Eigen::Vector3d& dv_octide,
               const Eigen::Vector3d& dx_poltide, const Eigen::Vector3d& dv_poltide,
               const Eigen::Vector3d& dx_atm, const Eigen::Vector3d& dv_atm,
               const Eigen::Vector3d& dx_temp, const Eigen::Vector3d& dv_temp,
               Eigen::Vector3d& xsta_j2000t, Eigen::Vector3d& vsta_j2000t,
               Eigen::Vector3d& asta_j2000) {
    // Простая сумма
    xsta_j2000t = xsta_itrf + dx_tide + dx_octide + dx_poltide + dx_atm + dx_temp;
    vsta_j2000t = vsta_itrf + dv_tide + dv_octide + dv_poltide + dv_atm + dv_temp;
    asta_j2000.setZero();
}
