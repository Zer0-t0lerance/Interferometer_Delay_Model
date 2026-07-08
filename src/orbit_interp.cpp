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
#include "..\\external\\spline.h"

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

    // Аргумент — время в сутках; узлы должны быть отсортированы по возрастанию.
    std::vector<double> t(n);
    std::vector<double> px(n), py(n), pz(n);
    for (int i = 0; i < n; ++i) {
        t[i]  = static_cast<double>(orbit[i].mjd) + orbit[i].utc;
        px[i] = orbit[i].xyz.x() * 1000.0; // км -> м
        py[i] = orbit[i].xyz.y() * 1000.0;
        pz[i] = orbit[i].xyz.z() * 1000.0;
    }

    tk::spline sx, sy, sz;
    sx.set_points(t, px);
    sy.set_points(t, py);
    sz.set_points(t, pz);

    const double SPD = cnst::SECDAY; // перевод производных: /сут -> /с
    x_j2000 << sx(mjd_utc), sy(mjd_utc), sz(mjd_utc);
    v_j2000 << sx.deriv(1, mjd_utc) / SPD, sy.deriv(1, mjd_utc) / SPD, sz.deriv(1, mjd_utc) / SPD;
    a_j2000 << sx.deriv(2, mjd_utc) / (SPD * SPD), sy.deriv(2, mjd_utc) / (SPD * SPD),
               sz.deriv(2, mjd_utc) / (SPD * SPD);
}

} // namespace ariadna
