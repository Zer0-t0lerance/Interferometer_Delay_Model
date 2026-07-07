// test_pipeline_obs.cpp
//
// Сверка EOP-НЕЗАВИСИМОГО "переда" конвейера для реального наблюдения
// (сессия 18JAN02XA, база FORTLEZA-HART15M, источник 0017+200, эпоха MJD 58120.7088)
// против входов дампа theor_delay:
//   * source_vec  -> K_s
//   * jpleph      -> Earth (SSB), Sun/Moon (геоцентрические)
// Эти величины не зависят от EOP (только от эпохи), поэтому сверяются сразу.
//
// Сборка (из корня репозитория):
//   g++ -std=c++17 -I.\external\eigen -I.\external -I.\external\dephem-master\include ^
//       -I.\external\sofa\20190722\c\src ^
//       src\jpleph.cpp src\ephemeris.cpp src\t_eph40.cpp src\site_functions.cpp src\GEOD.cpp ^
//       .\external\dephem-master\include\dephem\EphemerisRelease.cpp ^
//       .\external\sofa\20190722\c\build\libsofa.a ^
//       tests\test_pipeline_obs.cpp -o tests\test_pipeline_obs.exe

#include "../src/functions.h"
#include <cstdio>
#include <cmath>
#include <fstream>
#include <string>
#include <vector>

using namespace ariadna;

static int g_fail = 0;
static void chk(const char* name, double got, double exp, double tol) {
    double d = std::fabs(got - exp);
    bool ok = d <= tol;
    if (!ok) ++g_fail;
    printf("  %-16s calc=% .10e  dump=% .10e  |d|=%.2e  %s\n", name, got, exp, d, ok ? "OK" : "FAIL");
}

static std::string find_eph() {
    const char* c[] = {"external/dephem-master/linux_p1550p2650.440t",
                       "../external/dephem-master/linux_p1550p2650.440t"};
    for (auto p : c) { std::ifstream f(p, std::ios::binary); if (f.good()) return p; }
    return c[0];
}

int main() {
    printf("=====================================================================\n");
    printf("  UNIT: конвейер (перед) vs дамп theor_delay — 18JAN02XA, 0017+200\n");
    printf("---------------------------------------------------------------------\n");

    // --- source_vec: источник 0017+200 (RA/Dec из каталога сессии) ---
    std::vector<Source> src(1);
    src[0].name = "0017+200";
    src[0].ra = 0.085655996508366;
    src[0].dec = 0.355395794394430;
    src[0].ra_rate = 0.0;
    src[0].dec_rate = 0.0;

    double t_mean = 58120.7088; // MJD (для собственного движения; здесь нулевого)
    std::vector<Eigen::Vector3d> k_star;
    source_vec(src, t_mean, k_star);

    printf("  [source_vec -> K_s]\n");
    chk("K_s.x", k_star[0].x(), 0.9340717157118170, 5e-6);
    chk("K_s.y", k_star[0].y(), 0.08020509321046379, 5e-6);
    chk("K_s.z", k_star[0].z(), 0.3479614532247550, 5e-6);

    // --- Эфемериды: theor_delay ждёт Солнце/Луну в SSB -> get_celestial_bodies. ---
    try { init_ephemeris(find_eph()); }
    catch (const std::exception& e) { printf("SKIP: эфемериды: %s\n", e.what()); return 2; }

    // Валидированная эпоха этого наблюдения (из test_ephemeris): jd0=2458120.5, ct=0.7113110196885.
    double jd0 = 2458120.5;
    double ct = 0.7113110196885;

    Eigen::Matrix3d Earth, Sun, Moon; // все в SSB (поз/скор/ускор)
    get_celestial_bodies(jd0, ct, Earth, Sun, Moon);

    // Проверка КОНВЕНЦИИ (SSB, не геоцентр!) и порядка. Точная эпоха наблюдения
    // неизвестна (~сек), поэтому допуск 1e7 м (~5 мин): геоцентрический вариант
    // отличался бы на ~1.5e11 м. Точная сверка эфемерид — в test_ephemeris (vs Fortran DE).
    printf("\n  [get_celestial_bodies -> Earth/Sun/Moon в SSB; конвенция+порядок, допуск 1e7 м]\n");
    chk("Earth.x", Earth(0,0), -30331171034.03649, 1e7);
    chk("Earth.y", Earth(1,0),  132858910822.4568, 1e7);
    chk("Earth.z", Earth(2,0),   57575534750.98912, 1e7);
    chk("Sun.x(SSB)", Sun(0,0),    268184443.84710, 1e7);
    chk("Sun.y(SSB)", Sun(1,0),    849738216.96790, 1e7);
    chk("Moon.x(SSB)",Moon(0,0),-30457209591.74120, 1e7);
    chk("Moon.y(SSB)",Moon(1,0),133170659155.41185, 1e7);

    printf("=====================================================================\n");
    printf("  РЕЗУЛЬТАТ: %s (провалов: %d)\n", g_fail == 0 ? "PASS" : "FAIL", g_fail);
    printf("  [геометрия (site/baseline/r2000) требует EOP на MJD 58120 — следующий шаг]\n");
    return g_fail == 0 ? 0 : 1;
}
