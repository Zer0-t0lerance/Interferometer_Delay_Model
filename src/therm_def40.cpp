// therm_def40.cpp
//
// Порт подпрограммы THERM_DEF40 (SUBROUTINES/THERM_DEF40.FOR) — термодеформация
// антенны: смещение и скорость опорной точки станции в J2000.0 из-за теплового
// расширения бетонного фундамента и колонны (модель Nothnagel 2009, J. Geodesy 83).
//
// Оригинальный THERM_DEF40 использует ТОЛЬКО высоты фундамента hf и колонны hp с
// их коэффициентами расширения gamma_hf/gamma_hp (в исходнике члены hv/hs/hd/AO
// закомментированы). Опорная точка (пересечение осей) поднимается по вертикали на
//   d_up = (gamma_hf*hf + gamma_hp*hp) * (T - T0),
// где T — локальная температура, T0 — опорная (ANTENNA_INFO). Разность в °C и K
// одинакова, поэтому переводить в K не требуется. Мобильным антеннам (нет в
// каталоге) заданы hf=1, hp=5, gamma_f=1e-5, gamma_a=1.2e-5 (см. DefPar по умолч.).
//
// Смещение задаётся в топоцентрической системе VEN (Vertical, East, North) как
// чисто вертикальное, затем VEN -> ITRF (vw_i) -> J2000 (r2000), для скорости — с
// учётом вращения Земли (dr2000_dt) и производной температуры dTdt (из DMETEO2_DT).
// Конвенция полностью совпадает с SITE_ATM40 / POLE_TIDE / SITE_TIDE_OC.
// ЕДИНИЦЫ: dTdt в °C/СУТ (как dPdt в SITE_ATM40) — вертикальная скорость делится
// на SECDAY для перевода в м/с.
//
// Частные производные (для МНК) не вычисляются — вне области геометрических задержек.

#include "functions.h"

namespace ariadna {

void THERM_DEF40(double tC, double dTdt, const DefPar& dp,
                 const Eigen::Matrix3d& vw_i,
                 const Eigen::Matrix3d& r2000, const Eigen::Matrix3d& dr2000_dt,
                 Eigen::Vector3d& dx_temp, Eigen::Vector3d& dv_temp)
{
    // Эффективный коэффициент вертикального расширения [м/°C].
    const double kappa = dp.gamma_hf * dp.hf + dp.gamma_hp * dp.hp;

    // Вертикальное смещение опорной точки [м] и его скорость [м/с]
    // (dTdt в °C/сут -> делим на SECDAY).
    const double d_up  = kappa * (tC - dp.t_0);
    const double dv_up = kappa * dTdt / cnst::SECDAY;

    // Топоцентрический вектор VEN: только вертикаль (col 0 = Vertical/Up).
    Eigen::Vector3d dr_ven(d_up, 0.0, 0.0);
    Eigen::Vector3d dv_ven(dv_up, 0.0, 0.0);

    // VEN -> ITRF -> J2000 (та же конвенция, что в SITE_ATM40).
    Eigen::Vector3d dx_itrf = vw_i * dr_ven;
    Eigen::Vector3d dv_itrf = vw_i * dv_ven;

    dx_temp = r2000 * dx_itrf;
    dv_temp = r2000 * dv_itrf + dr2000_dt * dx_itrf;
}

} // namespace ariadna
