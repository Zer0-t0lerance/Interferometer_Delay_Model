// theor_delay.cpp
//
// Порт SUBROUTINE THEOR_DELAYcorr (версия, которую реально вызывает главная
// программа ARIADNA) — теоретическая ВАКУУМНАЯ задержка, скорректированная за
// релятивистские эффекты (гравитационная задержка Шапиро от Солнца/Луны/Земли +
// геометрия с поправками на движение станций), плюс инструментальная задержка
// смещения оси монтировки (dtau_off).
//
// Отличия от прежнего порта THEOR_DELAY4_10:
//   * добавлены члены term2e/term2f («New terms for Radioastron»);
//   * тропосфера ВНУТРИ не добавляется (в corr Datmc обнулена) — её прибавляет
//     вызывающий код (trop_delay) отдельно; поэтому Datmc_d/Datmc_w убраны из сигнатуры;
//   * термодеформация (dt_temp) убрана — в corr она учитывается в координатах станций
//     (через THERM_DEF), а не в задержке;
//   * финальная сборка суммирует dtau_off ОБЕИХ станций: dtau_off(1,1)+dtau_off(2,1).
//
// Планеты (гравитация Шапиро от планет) пока не считаются (эффект ~пс).

#include "functions.h"
#include <cmath>

namespace ariadna {

void theor_delay(const Eigen::Matrix<double, 3, 2>& base_line,
                 const std::vector<Eigen::Vector3d>& xsta,
                 const std::vector<Eigen::Vector3d>& vsta,
                 const std::vector<Eigen::Vector3d>& asta,
                 const Eigen::Vector3d& K_s,
                 const Eigen::Matrix3d& Earth,
                 const Eigen::Matrix3d& Sun,
                 const Eigen::Matrix3d& Moon,
                 const Eigen::Matrix2d& dtau_off,
                 double& t2_t1, double& dt2_t1) {

    double C = cnst::C;
    double C2 = C * C;
    double C3 = C2 * C;
    double Gamma = 1.0;

    Eigen::Vector3d b = base_line.col(0);
    Eigen::Vector3d db = base_line.col(1);

    double K_starB = K_s.dot(b);
    double dotK_starB = K_s.dot(db);

    Eigen::Vector3d X_1 = Earth.col(0) + xsta[0];
    Eigen::Vector3d X_2 = Earth.col(0) + xsta[1];
    Eigen::Vector3d dotX_1 = Earth.col(1) + vsta[0];
    Eigen::Vector3d dotX_2 = Earth.col(1) + vsta[1];

    Eigen::Vector3d R_geocen = Earth.col(0) - Sun.col(0);
    Eigen::Vector3d V_geocen = Earth.col(1) - Sun.col(1);
    double R_geo = R_geocen.norm();

    double U_Sun = cnst::GSUN / C2;
    double U = U_Sun / R_geo;
    double dR_geo = R_geocen.dot(V_geocen) / R_geo;
    double dU = -(U_Sun / (R_geo * R_geo)) * dR_geo;

    // ===================== SHAPIRO SUN =====================
    double C_Sun = (1.0 + Gamma) * cnst::GSUN / C3;
    double t1_Sun = K_s.dot(Sun.col(0) - X_1) / C;
    if (t1_Sun < 0.0) t1_Sun = 0.0;

    Eigen::Vector3d X1_Sun_t1 = Sun.col(0) - Sun.col(1) * t1_Sun;
    Eigen::Vector3d R1_Sun_t1 = X_1 - X1_Sun_t1;
    Eigen::Vector3d R2_Sun_t1 = X_2 - Earth.col(1) * (K_starB / C) - X1_Sun_t1;

    double w1_s = R1_Sun_t1.norm();
    double numer1 = w1_s + K_s.dot(R1_Sun_t1);
    double w2_s = R2_Sun_t1.norm();
    double denom1 = w2_s + K_s.dot(R2_Sun_t1);
    double delta_t_grav_Sun = C_Sun * std::log(numer1 / denom1);

    double C_Sun1 = C_Sun * cnst::GSUN / C2;
    Eigen::Vector3d N_hat = R1_Sun_t1 / w1_s;
    Eigen::Vector3d NplusK = N_hat + K_s;
    double V1 = b.dot(NplusK);
    double V2 = numer1 * numer1;
    double add_grav_Sun1 = C_Sun1 * V1 / V2;

    // ===================== SHAPIRO MOON =====================
    double C_Moon = (1.0 + Gamma) * cnst::GMOON / C3;
    double t1_Moon = K_s.dot(Moon.col(0) - X_1) / C;
    if (t1_Moon < 0.0) t1_Moon = 0.0;

    Eigen::Vector3d X1_Moon_t1 = Moon.col(0) - Moon.col(1) * t1_Moon;
    Eigen::Vector3d R1_Moon_t1 = X_1 - X1_Moon_t1;
    Eigen::Vector3d R2_Moon_t1 = X_2 - Earth.col(1) * (K_starB / C) - X1_Moon_t1;

    double w1_m = R1_Moon_t1.norm();
    double numer2 = w1_m + K_s.dot(R1_Moon_t1);
    double w2_m = R2_Moon_t1.norm();
    double denom2 = w2_m + K_s.dot(R2_Moon_t1);
    double delta_t_grav_Moon = C_Moon * std::log(numer2 / denom2);

    // ===================== SHAPIRO EARTH =====================
    double C_Earth_grav = (1.0 + Gamma) * cnst::GEARTH / C3;
    double w1_e = xsta[0].norm();
    double numer4 = w1_e + K_s.dot(xsta[0]);
    double w2_e = xsta[1].norm();
    double denom4 = w2_e + K_s.dot(xsta[1]);
    double delta_t_grav_Earth = C_Earth_grav * std::log(numer4 / denom4);

    double delta_t_grav_Pl = 0.0; // планеты не учитываются
    double delta_t_grav = delta_t_grav_Sun + delta_t_grav_Moon + delta_t_grav_Pl + delta_t_grav_Earth;
    double term1 = delta_t_grav + add_grav_Sun1;

    // ===================== DTERM1 EARTH =====================
    double dR1_Earth = xsta[0].dot(vsta[0]) / w1_e;
    double dR2_Earth = xsta[1].dot(vsta[1]) / w2_e;
    double dnum4 = dR1_Earth + K_s.dot(vsta[0]);
    double dden4 = dR2_Earth + K_s.dot(vsta[1]);
    double dterm1_Earth = C_Earth_grav * (dnum4 / numer4 - dden4 / denom4);

    // ===================== DTERM1 SUN =====================
    Eigen::Vector3d dR1_Sun = dotX_1 - Sun.col(1);
    Eigen::Vector3d dR2_Sun = dotX_2 - Sun.col(1);
    double dw1_s = R1_Sun_t1.dot(dR1_Sun) / w1_s;
    double dw2_s = R2_Sun_t1.dot(dR2_Sun) / w2_s;
    double dnum1 = dw1_s + K_s.dot(dR1_Sun);
    double dden1 = dw2_s + K_s.dot(dR2_Sun);
    double dterm1_Sun = C_Sun * (dnum1 / numer1 - dden1 / denom1);

    // ===================== DTERM1 MOON =====================
    Eigen::Vector3d dR1_Moon = dotX_1 - Moon.col(1);
    Eigen::Vector3d dR2_Moon = dotX_2 - Moon.col(1);
    double dw1_m = R1_Moon_t1.dot(dR1_Moon) / w1_m;
    double dw2_m = R2_Moon_t1.dot(dR2_Moon) / w2_m;
    double dnum2 = dw1_m + K_s.dot(dR1_Moon);
    double dden2 = dw2_m + K_s.dot(dR2_Moon);
    double dterm1_Moon = C_Moon * (dnum2 / numer2 - dden2 / denom2);

    double dterm1_Pl = 0.0;
    double ddelta_t_grav = dterm1_Sun + dterm1_Moon + dterm1_Pl + dterm1_Earth;
    double dterm1 = ddelta_t_grav;

    // ===================== TERM 2 =====================
    double term2a = K_starB / C;
    double dterm2a = dotK_starB / C;

    double term2b = 1.0 - (1.0 + Gamma) * U;
    double dterm2b = -(1.0 + Gamma) * dU;

    double term2c = Earth.col(1).dot(Earth.col(1)) / (2.0 * C2);
    double dterm2c = Earth.col(1).dot(Earth.col(2)) / C2;

    double term2d = Earth.col(1).dot(vsta[1]) / C2;
    double dterm2d = (Earth.col(2).dot(vsta[1]) + Earth.col(1).dot(asta[1])) / C2;

    // --- Новые члены для Радиоастрона (THEOR_DELAYcorr) ---
    double term2e = Earth.col(2).dot(xsta[1]) / C2;
    double dterm2e = Earth.col(2).dot(vsta[1]) / C2;

    double bracket2f = K_s.dot(asta[1]) + K_s.dot(Earth.col(2));
    double term2f = term2a * bracket2f / (2.0 * C);
    double dterm2f = dterm2a * bracket2f / (2.0 * C);

    double term2bcd = term2b - term2c - term2d - term2e + term2f;
    double dterm2bcd = dterm2b - dterm2c - dterm2d - dterm2e + dterm2f;

    double term2 = term2a * term2bcd;
    double dterm2 = term2a * dterm2bcd + dterm2a * term2bcd;

    // ===================== TERM 3 =====================
    double term3a = Earth.col(1).dot(b) / C2;
    double term3b = 1.0 + K_s.dot(Earth.col(1)) / (2.0 * C);
    double term3 = term3a * term3b;

    double dterm3a = (Earth.col(2).dot(b) + Earth.col(1).dot(db)) / C2;
    double dterm3b = K_s.dot(Earth.col(2)) / (2.0 * C);
    double dterm3 = dterm3a * term3b + term3a * dterm3b;

    // ===================== EQ 9 =====================
    double numer_Eq9 = term1 - term2 - term3;
    double dnumer_Eq9 = dterm1 - dterm2 - dterm3;

    Eigen::Vector3d vec_sum = Earth.col(1) + vsta[1];
    Eigen::Vector3d dvec_sum = Earth.col(2) + asta[1];

    double den_Eq9 = 1.0 + K_s.dot(vec_sum) / C;
    double dden_Eq9 = K_s.dot(dvec_sum) / C;

    double tv2_tv1 = numer_Eq9 / den_Eq9;
    double dtv2_tv1 = dnumer_Eq9 / den_Eq9 - numer_Eq9 * dden_Eq9 / (den_Eq9 * den_Eq9);

    // ===================== ИНСТРУМЕНТАЛЬНАЯ ЗАДЕРЖКА (ось монтировки) =====================
    // Тропосфера здесь НЕ добавляется (в corr Datmc обнулена) — её прибавляет
    // оркестрация через trop_delay. Складываем задержку смещения оси обеих станций.
    t2_t1  = tv2_tv1  + dtau_off(0, 0) + dtau_off(1, 0);
    dt2_t1 = dtv2_tv1 + dtau_off(0, 1) + dtau_off(1, 1);
}

} // namespace ariadna
