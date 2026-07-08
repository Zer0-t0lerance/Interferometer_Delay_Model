// test_therm_def40.cpp
//
// Юнит-тест THERM_DEF40 (термодеформация антенны, Nothnagel 2009). Проверяет:
//   1) величину вертикального смещения d_up = (gamma_hf*hf + gamma_hp*hp)*(T - T0);
//   2) сохранение длины при VEN->ITRF->J2000 (r2000, vw ортонормированы): |dx_temp| = |d_up|;
//   3) направление: dx_temp || r2000*vw.col(0) (вертикаль в J2000);
//   4) скорость dv_temp от dTdt при r2000=I, dr2000=0 -> vw*(kappa*dTdt, 0, 0);
//   5) физичность: для HART15M (dT~9.8°C) смещение ~ доли мм.

#include "../src/functions.h"
#include <cstdio>
#include <cmath>

using namespace ariadna;

int main() {
    printf("=====================================================================\n");
    printf("  UNIT: THERM_DEF40 -- термодеформация антенны (Nothnagel 2009)\n");
    printf("---------------------------------------------------------------------\n");

    // HART15M: t_0=16.1, hf=6.32/gamma=0.8e-5, hp=3.36/gamma=1.2e-5.
    DefPar dp;
    dp.t_0 = 16.1; dp.hf = 6.32; dp.gamma_hf = 0.8e-5; dp.hp = 3.36; dp.gamma_hp = 1.2e-5;
    const double tC = 25.884;      // локальная температура наблюдения [°C]
    const double dTdt = 0.0;

    // Простая геодезия для vw (широта/долгота HART15M примерно).
    const double lat = -0.4520, lon = 0.4859; // рад
    double cla = std::cos(lat), sla = std::sin(lat), clo = std::cos(lon), slo = std::sin(lon);
    Eigen::Matrix3d vw;
    vw.col(0) << cla * clo, cla * slo, sla;   // Up
    vw.col(1) << -slo, clo, 0.0;              // East
    vw.col(2) << -sla * clo, -sla * slo, cla; // North

    // (1)(2)(3) при r2000 = I, dr2000 = 0.
    Eigen::Matrix3d I = Eigen::Matrix3d::Identity(), Z = Eigen::Matrix3d::Zero();
    Eigen::Vector3d dx, dv;
    THERM_DEF40(tC, dTdt, dp, vw, I, Z, dx, dv);

    const double kappa = dp.gamma_hf * dp.hf + dp.gamma_hp * dp.hp;
    const double d_up = kappa * (tC - dp.t_0);
    printf("  kappa = %.6e м/°C, dT = %.3f °C -> d_up = %.6e м (%.4f мм)\n",
           kappa, tC - dp.t_0, d_up, d_up * 1e3);

    double len_err = std::fabs(dx.norm() - std::fabs(d_up));
    Eigen::Vector3d up_j2000 = vw.col(0); // r2000=I
    double dir_err = (dx - d_up * up_j2000).norm();
    printf("  |dx_temp| = %.6e (ожидание |d_up| = %.6e), |Δ|=%.2e\n", dx.norm(), std::fabs(d_up), len_err);
    printf("  направление вдоль Up: |dx - d_up*Up| = %.2e\n", dir_err);

    // (4) скорость: dTdt != 0, r2000=I, dr2000=0 -> dv = vw*(kappa*dTdt/SECDAY,0,0).
    Eigen::Vector3d dx2, dv2;
    double dTdt2 = 5.0; // °C/сут
    THERM_DEF40(tC, dTdt2, dp, vw, I, Z, dx2, dv2);
    Eigen::Vector3d dv_expect = vw * Eigen::Vector3d(kappa * dTdt2 / cnst::SECDAY, 0.0, 0.0);
    double v_err = (dv2 - dv_expect).norm();
    printf("  скорость: |dv_temp - vw*(kappa*dTdt,0,0)| = %.2e\n", v_err);

    // (5) физичность
    bool phys = std::fabs(d_up) < 2.0e-3 && std::fabs(d_up) > 1.0e-5; // 0.01..2 мм

    bool ok = len_err < 1e-15 && dir_err < 1e-15 && v_err < 1e-18 && phys;
    printf("---------------------------------------------------------------------\n");
    printf("  РЕЗУЛЬТАТ: %s\n", ok ? "PASS" : "FAIL");
    return ok ? 0 : 1;
}
