// rotation.cpp
//
// Матрица перехода ITRF <-> J2000 (GCRS) по стандарту IAU 2006/2000A через SOFA.
//
// Заменяет прежнюю упрощённую реализацию по IAU 1980. Используется официальная SOFA
// (external/sofa/.../c): iauC2t06a даёт матрицу GCRS -> ITRS (celestial-to-terrestrial),
// поэтому R2000 (ITRF -> GCRS) = транспонирование этой матрицы.
//
// ЕДИНИЦЫ ВРЕМЕНИ: iauC2t06a ожидает TT для прецессии/нутации и UT1 для вращения Земли.
// На вход подаётся JD_TDB; разница TDB-TT < 1.7 мс влияет на матрицу на уровне ~1e-11 и
// здесь принимается пренебрежимо малой (TDB используется как TT). Двухчастное задание JD
// (2400000.5, MJD) сохраняет точность и воспроизводит эталонные тесты SOFA бит-в-бит.
//
// Производная dR2000/dt считается центральной разностью по времени (шаг 1 с): доминирует
// суточное вращение Земли; движение полюса (xp, yp) считается постоянным на этом интервале.

#include "functions.h"

extern "C" {
#include "sofa.h"
}

namespace ariadna {

// R2000 = ITRF -> GCRS(J2000) на момент (jd_tt, jd_ut1) при движении полюса (xp, yp) [рад].
static Eigen::Matrix3d compute_r2000(double jd_tt, double jd_ut1, double xp, double yp) {
    const double djmjd0 = 2400000.5; // опорная точка двухчастного JD (как в тестах SOFA)
    double rc2t[3][3];
    iauC2t06a(djmjd0, jd_tt - djmjd0, djmjd0, jd_ut1 - djmjd0, xp, yp, rc2t); // GCRS -> ITRS

    Eigen::Matrix3d R; // R2000 = transpose(rc2t): ITRS -> GCRS
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
            R(i, j) = rc2t[j][i];
    return R;
}

void get_r2000_matrices(double JD_TDB, double JD_UT1, double xp, double yp,
                        Eigen::Matrix3d& R2000, Eigen::Matrix3d& dR2000_dt) {
    R2000 = compute_r2000(JD_TDB, JD_UT1, xp, yp);

    // Первая производная по времени — центральная разность, шаг 1 секунда.
    const double dt_sec = 1.0;
    const double dday = dt_sec / cnst::SECDAY;
    Eigen::Matrix3d Rp = compute_r2000(JD_TDB + dday, JD_UT1 + dday, xp, yp);
    Eigen::Matrix3d Rm = compute_r2000(JD_TDB - dday, JD_UT1 - dday, xp, yp);
    dR2000_dt = (Rp - Rm) / (2.0 * dt_sec); // [1/с]
}

} // namespace ariadna
