// test_orbit_space.cpp
//
// Структурный тест наземно-космического пути (без орбитального файла — «вслепую»):
//   1) orbit_interp восстанавливает аналитическую орбиту (положение/скорость/ускорение);
//   2) site_pair для космической станции (is_space) берёт координаты НАПРЯМУЮ из орбиты,
//      минуя наземные поправки (приливы/нагрузки/термо).
//
// Валидации против ARIADNA нет (нет эфемериды орбиты .scf); проверяется корректность
// интеграции космического пути в конвейер.

#include "../src/functions.h"
#include <cstdio>
#include <cmath>

using namespace ariadna;

int main() {
    printf("=====================================================================\n");
    printf("  UNIT: orbit_interp + site_pair(is_space) -- наземно-космический путь\n");
    printf("---------------------------------------------------------------------\n");
    bool ok = true;

    // --- 1. Аналитическая орбита: r(t) = r0 + v0*(t-t0) + 0.5*a0*(t-t0)^2 [км, сутки] ---
    // Проверяем, что сплайн восстановит положение/скорость/ускорение в СИ.
    const double t0 = 58120.0;
    Eigen::Vector3d r0(7000.0, -3000.0, 1500.0);       // км
    Eigen::Vector3d v0_pd(120.0, 80.0, -50.0);         // км/сут
    Eigen::Vector3d a0_pd(3.0, -2.0, 1.0);             // км/сут^2
    std::vector<SpaceStation> orbit(9);
    for (int i = 0; i < 9; ++i) {
        double dt = (i - 4) * 0.02;                    // ±0.08 сут вокруг t0
        orbit[i].mjd = 58120; orbit[i].utc = 0.5 + dt; // t = t0+0.5+dt (сутки)
        double s = orbit[i].utc - 0.5;                 // (t - t0центр)
        orbit[i].xyz = r0 + v0_pd * s + 0.5 * a0_pd * s * s;
    }
    double t_q = 58120.5 + 0.017;                      // произвольный момент внутри
    double s_q = 0.017;
    Eigen::Vector3d x, v, a;
    orbit_interp(orbit, t_q, x, v, a);

    Eigen::Vector3d x_exp = (r0 + v0_pd * s_q + 0.5 * a0_pd * s_q * s_q) * 1000.0;    // м
    Eigen::Vector3d v_exp = (v0_pd + a0_pd * s_q) * 1000.0 / cnst::SECDAY;            // м/с
    Eigen::Vector3d a_exp = a0_pd * 1000.0 / (cnst::SECDAY * cnst::SECDAY);           // м/с^2
    double ex = (x - x_exp).norm(), ev = (v - v_exp).norm(), ea = (a - a_exp).norm();
    printf("  orbit_interp: |dx|=%.3e м  |dv|=%.3e м/с  |da|=%.3e м/с^2\n", ex, ev, ea);
    // Пороги физичные: натуральный сплайн на СИНТЕТИЧЕСКОЙ разреженной квадратичной
    // орбите даёт краевую ошибку ~мм (на реальной плотной .scf — на порядки меньше).
    ok = ok && ex < 5e-3 && ev < 1e-5 && ea < 5e-8;

    // --- 2. site_pair: обе станции космические -> координаты = орбита (поправки минуются) ---
    SitePrep s1, s2;
    s1.is_space = true; s1.x_orbit = Eigen::Vector3d(1e7, 2e7, 3e7);
    s1.v_orbit = Eigen::Vector3d(100, 200, 300); s1.a_orbit = Eigen::Vector3d(1, 2, 3);
    s2.is_space = true; s2.x_orbit = Eigen::Vector3d(-4e6, 5e6, -6e6);
    s2.v_orbit = Eigen::Vector3d(-40, 50, -60); s2.a_orbit = Eigen::Vector3d(-4, 5, -6);

    // Аргументы среды не используются (космическая ветвь site_one возвращается сразу).
    Eigen::VectorXd f(5), fd(5); f.setZero(); fd.setZero();
    Eigen::Vector2d gast(0.0, cnst::OMEGA_EARTH);
    Eigen::Matrix<double, 3, 2> sun_geo = Eigen::Matrix<double, 3, 2>::Zero(), moon_geo = sun_geo;
    Eigen::Matrix3d I = Eigen::Matrix3d::Identity(), Z = Eigen::Matrix3d::Zero();
    std::vector<Eigen::Vector3d> xs, vs, as_;
    site_pair(s1, s2, 58120, 0.5, 2458120.5, 43200.0, 0.16, f, fd, gast, sun_geo, moon_geo,
              0.0, 0.0, 0.0, 0.0, I, Z, Z, xs, vs, as_);

    double d1 = (xs[0] - s1.x_orbit).norm() + (vs[0] - s1.v_orbit).norm() + (as_[0] - s1.a_orbit).norm();
    double d2 = (xs[1] - s2.x_orbit).norm() + (vs[1] - s2.v_orbit).norm() + (as_[1] - s2.a_orbit).norm();
    printf("  site_pair(is_space): |станция1 - орбита|=%.3e  |станция2 - орбита|=%.3e\n", d1, d2);
    ok = ok && d1 < 1e-12 && d2 < 1e-12;

    printf("---------------------------------------------------------------------\n");
    printf("  РЕЗУЛЬТАТ: %s\n", ok ? "PASS" : "FAIL");
    return ok ? 0 : 1;
}
