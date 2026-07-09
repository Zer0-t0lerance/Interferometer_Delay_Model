// orbit_interp.cpp
//
// Интерполяция орбиты космического телескопа (RASTRON) на момент наблюдения.
// Орбита приходит извне готовой эфемеридой (.scf) — таблицей точек с некоторым шагом
// (время, положение, скорость в GCRS/J2000). Между точками положение/скорость/ускорение
// получаем кубическим сплайном (интегратор орбиты Эверхарта НЕ нужен — так решено).
//
// Кадр: орбита задана в GCRS ≈ J2000 (небесная система), поворот r2000 НЕ применяется —
// положение космической станции берётся напрямую (в отличие от наземной ITRF->J2000).
// Единицы входа: SpaceStation.xyz [км], .vel [км/с]; на выход — метры и м/с, м/с^2.

#include "functions.h"
#include "../external/spline.h"

namespace ariadna {

void orbit_interp(const std::vector<SpaceStation>& orbit, double mjd_utc,
                  Eigen::Vector3d& x_j2000, Eigen::Vector3d& v_j2000, Eigen::Vector3d& a_j2000) {
    const int n = static_cast<int>(orbit.size());
    if (n == 0) { x_j2000.setZero(); v_j2000.setZero(); a_j2000.setZero(); return; }
    if (n == 1) {
        x_j2000 = orbit[0].xyz * 1000.0;
        v_j2000 = orbit[0].vel * 1000.0;
        a_j2000 = orbit[0].acc * 1000.0;
        return;
    }

    // Ось времени — СЕКУНДЫ от первой точки орбиты (НЕ mjd+utc): иначе большие значения
    // (~56770 сут) при шаге 1с теряют точность в разностях сплайна -> ~1e-7 отн. ошибка,
    // что на радиусе орбиты ~2e8 м даёт сотни метров. Аргумент запроса — в тех же секундах.
    const double t0 = static_cast<double>(orbit[0].mjd) + orbit[0].utc; // сутки
    std::vector<double> t(n), px(n), py(n), pz(n);
    for (int i = 0; i < n; ++i) {
        t[i]  = ((static_cast<double>(orbit[i].mjd) - orbit[0].mjd) + (orbit[i].utc - orbit[0].utc)) * cnst::SECDAY;
        px[i] = orbit[i].xyz.x() * 1000.0; // км -> м
        py[i] = orbit[i].xyz.y() * 1000.0;
        pz[i] = orbit[i].xyz.z() * 1000.0;
    }

    tk::spline sx, sy, sz;
    sx.set_points(t, px);
    sy.set_points(t, py);
    sz.set_points(t, pz);

    const double q = (mjd_utc - t0) * cnst::SECDAY; // момент запроса, сек от старта орбиты
    x_j2000 << sx(q), sy(q), sz(q);
    v_j2000 << sx.deriv(1, q), sy.deriv(1, q), sz.deriv(1, q);         // м/с (ось уже в сек)
    a_j2000 << sx.deriv(2, q), sy.deriv(2, q), sz.deriv(2, q);         // м/с^2
}

} // namespace ariadna
