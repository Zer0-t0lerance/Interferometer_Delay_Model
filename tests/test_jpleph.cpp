// test_jpleph.cpp
//
// Unit-тест для обёрток jpleph / jpl_eph (соглашение JPLEPH_421):
//   Земля — барицентрическая (SSB); Солнце и Луна — ГЕОЦЕНТРИЧЕСКИЕ.
// Отдельного числового эталона нет (эфемериды проверяются официальным пакетом),
// поэтому тест самосогласованный + санитарные проверки масштаба.
//
// Сборка (из корня репозитория):
//   g++ -std=c++17 -Iexternal\eigen -Iexternal ^
//       src\jpleph.cpp src\ephemeris.cpp tests\test_jpleph.cpp -o tests\test_jpleph.exe
// Запуск: из корня репозитория (путь к файлу эфемерид относительный).

#include "../src/functions.h"
#include <cstdio>
#include <cmath>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>

using namespace ariadna;

// Находит файл вне зависимости от рабочей папки (корень репо или tests/).
static std::string resolve_path(const char* argv_override, std::vector<std::string> candidates) {
    if (argv_override) candidates.insert(candidates.begin(), argv_override);
    for (const auto& c : candidates) {
        std::ifstream f(c, std::ios::binary);
        if (f.good()) return c;
    }
    return candidates.back();
}

static int g_fail = 0;
static void expect(const char* what, bool ok, const char* detail = "") {
    if (!ok) ++g_fail;
    printf("  %-42s %s  %s\n", what, ok ? "OK" : "FAIL", detail);
}

int main(int argc, char** argv) {
    printf("=====================================================================\n");
    printf("  UNIT: jpleph / jpl_eph (frame convention + sanity)\n");
    printf("---------------------------------------------------------------------\n");

    std::string eph = resolve_path(argc > 1 ? argv[1] : nullptr,
        {"external/dephem-master/linux_p1550p2650.440t",
         "../external/dephem-master/linux_p1550p2650.440t"});
    printf("  Эфемериды: %s\n", eph.c_str());

    try {
        init_ephemeris(eph);
    } catch (const std::exception& e) {
        std::cerr << "SKIP/FAIL: не удалось загрузить эфемериды: " << e.what() << "\n";
        return 2;
    }

    const double jd = 2458120.5;
    const double ct = 0.7113110196885;

    Eigen::Matrix3d E, S, M;
    get_celestial_bodies(jd, ct, E, S, M); // все три в SSB

    Eigen::Matrix3d earth;
    Eigen::MatrixXd sun, moon;
    jpl_eph(jd, ct, earth, sun, moon);

    Eigen::Vector3d e_v, s_v, m_v;
    jpleph(jd, ct, e_v, s_v, m_v);

    const double tol = 1e-6; // м — обёртки должны совпадать с get_celestial_bodies бит-в-бит (до арифметики)

    // 1. Земля в jpl_eph = SSB Земля из get_celestial_bodies (поз/скор/ускор)
    expect("jpl_eph earth == SSB Earth (pos)", (earth.col(0) - E.col(0)).norm() < tol);
    expect("jpl_eph earth == SSB Earth (vel)", (earth.col(1) - E.col(1)).norm() < tol);
    expect("jpl_eph earth == SSB Earth (acc)", (earth.col(2) - E.col(2)).norm() < tol);

    // 2. Солнце/Луна в jpl_eph — геоцентрические (тело_SSB - Земля_SSB)
    expect("jpl_eph sun  == geocentric Sun (pos)",  (sun.col(0)  - (S.col(0) - E.col(0))).norm() < tol);
    expect("jpl_eph sun  == geocentric Sun (vel)",  (sun.col(1)  - (S.col(1) - E.col(1))).norm() < tol);
    expect("jpl_eph moon == geocentric Moon (pos)", (moon.col(0) - (M.col(0) - E.col(0))).norm() < tol);
    expect("jpl_eph moon == geocentric Moon (vel)", (moon.col(1) - (M.col(1) - E.col(1))).norm() < tol);

    // 3. Векторный jpleph == позиции матричного jpl_eph
    expect("jpleph earth == jpl_eph earth pos", (e_v - earth.col(0)).norm() < tol);
    expect("jpleph sun   == jpl_eph sun pos",   (s_v - sun.col(0)).norm()   < tol);
    expect("jpleph moon  == jpl_eph moon pos",  (m_v - moon.col(0)).norm()  < tol);

    // 4. Санитарные проверки масштаба
    double AU = 1.495978707e11;
    double sun_dist = s_v.norm();
    double moon_dist = m_v.norm();
    double earth_ssb = e_v.norm();
    char d1[64], d2[64], d3[64];
    std::snprintf(d1, sizeof(d1), "|Sun_geo|=%.4e m (~1 AU)", sun_dist);
    std::snprintf(d2, sizeof(d2), "|Moon_geo|=%.4e m (~3.8e8)", moon_dist);
    std::snprintf(d3, sizeof(d3), "|Earth_SSB|=%.4e m (~1 AU)", earth_ssb);
    expect("Sun geocentric distance ~ 1 AU", std::fabs(sun_dist - AU) < 0.05 * AU, d1);
    expect("Moon geocentric distance ~ 3.6..4.1e8 m", moon_dist > 3.4e8 && moon_dist < 4.1e8, d2);
    expect("Earth SSB distance ~ 1 AU", std::fabs(earth_ssb - AU) < 0.1 * AU, d3);

    printf("=====================================================================\n");
    printf("  РЕЗУЛЬТАТ: %s (провалов: %d)\n", g_fail == 0 ? "PASS" : "FAIL", g_fail);
    return g_fail == 0 ? 0 : 1;
}
