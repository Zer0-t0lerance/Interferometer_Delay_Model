// process_obs.cpp
//
// compute_delay_obs -- полная теоретическая задержка ОДНОГО наблюдения: связывает
// уже реализованные и сверенные с ARIADNA блоки:
//   site_pair -> baseline -> aber_source -> trop_delay -> mount_tel -> theor_delay.
//
// Все "средовые" величины (время, эфемериды, EOP, r2000, gast) подготавливает
// вызывающий код (process_ariadna) и передаёт сюда. Конвенции (проверены сверкой):
//   * Солнце/Луна для theor_delay -- SSB; для приливов (в site_pair) -- геоцентрические;
//   * Datmc (trop_delay) и dtau_off (mount_tel) ТРАНСПОНИРУЮТСЯ перед theor_delay
//     (у них [строка=станция], у theor_delay [строка=задержка/скорость]);
//   * термодеформация dt_temp пока не учитывается (нет THERM_DEF40).

#include "functions.h"

namespace ariadna {

void compute_delay_obs(const SitePrep& s1, const SitePrep& s2,
                       const Eigen::Vector3d& K_s, const Observation& obs,
                       int mjd, double utc, double jd0, double ct, double cent, double ut1_sec,
                       const Eigen::VectorXd& f, const Eigen::VectorXd& fd,
                       const Eigen::Vector2d& gast,
                       const Eigen::Matrix3d& Earth, const Eigen::Matrix3d& Sun, const Eigen::Matrix3d& Moon,
                       const Eigen::Matrix<double, 3, 2>& sun_geo, const Eigen::Matrix<double, 3, 2>& moon_geo,
                       double xp, double yp, double xp_rate, double yp_rate,
                       const Eigen::Matrix3d& r2000, const Eigen::Matrix3d& dr2000, const Eigen::Matrix3d& d2r2000,
                       double& tau, double& dtau) {
    // 1. Координаты станций в J2000 со всеми поправками.
    std::vector<Eigen::Vector3d> xs, vs, as_;
    site_pair(s1, s2, mjd, utc, jd0, ut1_sec, cent, f, fd, gast, sun_geo, moon_geo,
              xp, yp, xp_rate, yp_rate, r2000, dr2000, d2r2000, xs, vs, as_);

    // 2. Вектор базы.
    Eigen::MatrixXd Xm(3, 2), Vm(3, 2), Am(3, 2);
    Xm.col(0) = xs[0]; Xm.col(1) = xs[1];
    Vm.col(0) = vs[0]; Vm.col(1) = vs[1];
    Am.col(0) = as_[0]; Am.col(1) = as_[1];
    Eigen::MatrixXd blm(3, 2); Eigen::Vector3d bcfs;
    baseline(r2000, Xm, Vm, Am, blm, bcfs);

    // 3. Аберрация -> элевации/азимуты (earth: только скорость SSB).
    std::vector<Eigen::Matrix3d> r2000_v = {r2000, dr2000};
    std::vector<Eigen::Matrix3d> vw = {s1.vw_i, s2.vw_i};
    Eigen::Matrix3d Earth_vel = Eigen::Matrix3d::Zero(); Earth_vel.col(1) = Earth.col(1);
    Eigen::Matrix2d e2, az2;
    aber_source(obs, r2000_v, K_s, Earth_vel, vs, vw, e2, az2);

    // 4. Тропосфера.
    Station st1, st2;
    st1.lat_geod = s1.lat_geod; st1.h_geod = s1.h_geod;
    st2.lat_geod = s2.lat_geod; st2.h_geod = s2.h_geod;
    Eigen::MatrixXd eM = e2, azM = az2;
    Eigen::MatrixXd dd, dw, dhmf, dwmf, dgn, dge, zd, zw;
    double jd_full = jd0 + ct; // JD в шкале TT/TDB для тропосферы (day-of-year)
    trop_delay(obs, jd_full, 0.0, st1, st2, eM, azM, dd, dw, dhmf, dwmf, dgn, dge, zd, zw);

    // 5. Инструментальная задержка оси монтировки.
    Eigen::MatrixXd r2000_9(9, 3); r2000_9.setZero();
    r2000_9.block<3, 3>(0, 0) = r2000.transpose();
    r2000_9.block<3, 3>(3, 0) = dr2000.transpose();
    std::vector<Station> mst(2);
    mst[0].axsty = s1.axsty; mst[0].offs = s1.offs; mst[0].lat_geod = s1.lat_geod;
    mst[1].axsty = s2.axsty; mst[1].offs = s2.offs; mst[1].lat_geod = s2.lat_geod;
    std::vector<Eigen::Vector3d> ks{K_s};
    Eigen::MatrixXd doff_dl, d_dax, dtau_mt;
    mount_tel(obs, r2000_9, mst, ks, vw, eM, azM, doff_dl, d_dax, dtau_mt);

    // 6. Теоретическая задержка (Datmc и dtau_off транспонируем; Солнце/Луна SSB; dt_temp=0).
    Eigen::Matrix<double, 3, 2> bl; bl.col(0) = blm.col(0); bl.col(1) = blm.col(1);
    Eigen::Matrix2d Dd = Eigen::Matrix2d(dd).transpose();
    Eigen::Matrix2d Dw = Eigen::Matrix2d(dw).transpose();
    Eigen::Matrix2d Dtau = Eigen::Matrix2d(dtau_mt).transpose();
    Eigen::Matrix2d Dtmp = Eigen::Matrix2d::Zero(); // THERM_DEF40 не реализован
    theor_delay(bl, xs, vs, as_, K_s, Earth, Sun, Moon, Dd, Dw, Dtau, Dtmp, tau, dtau);
}

} // namespace ariadna
