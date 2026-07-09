// test_ra_model.cpp
//
// Диагностика модели RASTRON: сравнить с эталоном RA_L.TXT ЧИСТУЮ геометрию -(R·K)/c
// (R — положение из орбиты в J2000, K — единичный вектор на источник) в 13:00 и 13:10.
// Если чистая геометрия совпадает в ОБЕИХ точках — эталон геометрический (без релятивизма),
// и полный конвейер (с аберрацией/Шапиро) для быстрой орбиты даёт расхождение.

#ifdef _WIN32
#define WIN32_LEAN_AND_MEAN
#define NOMINMAX
#include <windows.h>
#endif
#include "../src/functions.h"
#include <cstdio>
#include <cmath>

using namespace ariadna;

int main() {
#ifdef _WIN32
    SetConsoleOutputCP(CP_UTF8);
#endif
    printf("=====================================================================\n");
    printf("  RASTRON: чистая геометрия -(R·K)/c vs эталон RA_L.TXT\n");
    printf("---------------------------------------------------------------------\n");

    std::vector<SpaceStation> orbit;
    if (!read_scf_orbit("example/RA140423_1200_v02.scf", orbit)) { printf("FAIL: scf\n"); return 1; }

    // K из RA/Dec (J2000), единичный вектор.
    double ra = 137.2920483250 * cnst::CDEGRAD, dec = 1.3598938083 * cnst::CDEGRAD;
    Eigen::Vector3d K(std::cos(dec) * std::cos(ra), std::cos(dec) * std::sin(ra), std::sin(dec));

    struct P { const char* lbl; double utc_sec; double ref; };
    P pts[] = {
        {"13:00:00 (блок 0)",  46800.0, -4.74842393132049e-1},
        {"13:10:00 (блок 10)", 47400.0, -4.77120881905411e-1},
    };
    int mjd = 56770;
    for (auto& p : pts) {
        double mjd_utc = mjd + p.utc_sec / cnst::SECDAY;
        Eigen::Vector3d x, v, a; orbit_interp(orbit, mjd_utc, x, v, a);
        double tau_geo = -(x.dot(K)) / cnst::C;                // чистая геометрия
        // Ретардация 1-го порядка: положение в момент прихода фронта t + tau.
        Eigen::Vector3d x_ret, vv, aa; orbit_interp(orbit, mjd_utc + tau_geo / cnst::SECDAY, x_ret, vv, aa);
        double tau_ret = -(x_ret.dot(K)) / cnst::C;
        printf("  %s:\n", p.lbl);
        printf("    чистая -(R·K)/c   = % .12e  |Δэт|=%.2e с (%.1f м)\n", tau_geo, std::fabs(tau_geo - p.ref), std::fabs(tau_geo - p.ref) * 3e8);
        printf("    с ретардацией     = % .12e  |Δэт|=%.2e с (%.1f м)\n", tau_ret, std::fabs(tau_ret - p.ref), std::fabs(tau_ret - p.ref) * 3e8);
        printf("    эталон P0         = % .12e\n", p.ref);
    }
    printf("---------------------------------------------------------------------\n");
    return 0;
}
