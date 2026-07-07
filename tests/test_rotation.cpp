// test_rotation.cpp
//
// Unit-тест цепочки вращения ITRF <-> J2000(GCRS) по IAU 2006/2000A (SOFA).
//
// Сборка (из корня репозитория):
//   g++ -std=c++17 -I.\external\eigen -I.\external\sofa\20190722\c\src ^
//       src\rotation.cpp tests\test_rotation.cpp ^
//       .\external\sofa\20190722\c\build\libsofa.a -o .\tests\test_rotation.exe
//
// Эталон R2000 = transpose(rc2t), где rc2t — эталонные значения из собственного
// теста SOFA t_c2t06a (t_sofa_c.c) для входа:
//   TT = UT1 = MJD 53736.0, xp=2.55060238e-7, yp=1.860359247e-6 [рад].
// Это внешняя истина (SOFA), совпадение проверяется до 1e-12.

#include "../src/functions.h"
#include <cstdio>
#include <cmath>

using namespace ariadna;

static int g_fail = 0;
static void check(const char* what, double got, double exp, double tol) {
    double d = std::fabs(got - exp);
    bool ok = d <= tol;
    if (!ok) ++g_fail;
    printf("  %-22s got=% .16e exp=% .16e |d|=%.1e %s\n",
           what, got, exp, d, ok ? "OK" : "FAIL");
}

int main() {
    printf("=====================================================================\n");
    printf("  UNIT: get_r2000_matrices (IAU 2006/2000A via SOFA)\n");
    printf("---------------------------------------------------------------------\n");

    // Вход как в тесте SOFA t_c2t06a.
    const double jd = 2400000.5 + 53736.0; // TT = UT1
    const double xp = 2.55060238e-7;       // рад
    const double yp = 1.860359247e-6;      // рад

    Eigen::Matrix3d R2000, dR2000_dt;
    get_r2000_matrices(jd, jd, xp, yp, R2000, dR2000_dt);

    // Эталон R2000 = transpose(rc2t) из SOFA t_c2t06a.
    Eigen::Matrix3d Rref;
    Rref << -0.1810332128305897282,  -0.9834768134136214897,   0.5773474024748545878e-3,
             0.9834769806938592296,  -0.1810332203649130832,   0.3961816829632690581e-4,
             0.6555550962998436505e-4, 0.5749800844905594110e-3, 0.9999998325501747785;

    printf("  [R2000 vs SOFA reference, tol=1e-12]\n");
    const char* nm[3][3] = {{"R(0,0)","R(0,1)","R(0,2)"},
                            {"R(1,0)","R(1,1)","R(1,2)"},
                            {"R(2,0)","R(2,1)","R(2,2)"}};
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
            check(nm[i][j], R2000(i, j), Rref(i, j), 1e-12);

    // Ортонормальность и det=1.
    double orth = (R2000 * R2000.transpose() - Eigen::Matrix3d::Identity()).norm();
    double det = R2000.determinant();
    printf("\n  [Ортонормальность]\n");
    printf("  ||R*R^T - I|| = %.3e  %s\n", orth, orth < 1e-12 ? "OK" : "FAIL");
    printf("  det(R) = %.15f      %s\n", det, std::fabs(det - 1.0) < 1e-12 ? "OK" : "FAIL");
    if (orth >= 1e-12) ++g_fail;
    if (std::fabs(det - 1.0) >= 1e-12) ++g_fail;

    // Физическая проверка dR/dt: S = dR/dt * R^T должна быть антисимметричной,
    // а вектор угловой скорости omega = (S21, S02, S10) ~ угловой скорости Земли.
    Eigen::Matrix3d S = dR2000_dt * R2000.transpose();
    double antisym = (S + S.transpose()).norm();
    double omega = std::sqrt(S(2,1)*S(2,1) + S(0,2)*S(0,2) + S(1,0)*S(1,0));
    const double OMEGA_EARTH = 7.292115e-5; // рад/с
    printf("\n  [Производная dR/dt]\n");
    printf("  ||S + S^T|| = %.3e  (антисимметрия)  %s\n", antisym, antisym < 1e-9 ? "OK" : "FAIL");
    printf("  |omega| = %.9e рад/с  (Земля ~%.6e)  %s\n", omega, OMEGA_EARTH,
           std::fabs(omega - OMEGA_EARTH) < 5e-8 ? "OK" : "FAIL");
    if (antisym >= 1e-9) ++g_fail;
    if (std::fabs(omega - OMEGA_EARTH) >= 5e-8) ++g_fail;

    printf("=====================================================================\n");
    printf("  РЕЗУЛЬТАТ: %s (провалов: %d)\n", g_fail == 0 ? "PASS" : "FAIL", g_fail);
    return g_fail == 0 ? 0 : 1;
}
