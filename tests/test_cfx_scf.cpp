// test_cfx_scf.cpp
//
// Проверка входного слоя: разбор задания коррелятора (.cfx) и орбиты (.scf) из example/.
// Печатает разобранные станции/источник/сканы и сводку по орбите; сверяет ключевые
// значения (MJD старта скана = 56770 для 2014-04-23; первая точка орбиты из шапки .scf).

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
    printf("  Разбор задания .cfx и орбиты .scf (example/)\n");
    printf("---------------------------------------------------------------------\n");
    bool ok = true;

    CfxTask task;
    const char* cfx = "example/GVLBI_RAKS01A0_L_20140423T130000_ASC_V1.cfx";
    if (!parse_cfx(cfx, task)) { printf("FAIL: не разобрать %s\n", cfx); return 1; }

    printf("Станции (%zu):\n", task.stations.size());
    for (const auto& s : task.stations) {
        if (s.is_space)
            printf("  %-9s iam=%-3s КОСМОС orb=%s poly=%s\n", s.name.c_str(), s.iam.c_str(),
                   s.orb_file.substr(s.orb_file.find_last_of("/\\")+1).c_str(), s.poly_file.c_str());
        else
            printf("  %-9s iam=%-3s XYZ=(%.1f,%.1f,%.1f) axoff=%.3f эпоха=%.0f %s poly=%s\n",
                   s.name.c_str(), s.iam.c_str(), s.xyz.x(), s.xyz.y(), s.xyz.z(),
                   s.axoff, s.epoch_mjd, s.mount.c_str(), s.poly_file.c_str());
    }
    printf("Источники (%zu):\n", task.sources.size());
    for (const auto& s : task.sources)
        printf("  %-9s RA=%.6f рад  Dec=%.6f рад\n", s.name.c_str(), s.ra, s.dec);
    printf("Сканы (%zu):\n", task.scans.size());
    for (const auto& sc : task.scans) {
        printf("  MJD %d utc=%.6f (%.0f с)  src=%-9s станции:", sc.mjd, sc.utc, sc.dur_sec, sc.source.c_str());
        for (const auto& t : sc.tel_iam) printf(" %s", t.c_str());
        printf("\n");
    }
    printf("Выходная папка (%%W): %s\n", task.out_path.c_str());

    // Сверки .cfx
    if (task.scans.empty() || task.scans[0].mjd != 56770) { printf("FAIL: MJD старта скана != 56770\n"); ok = false; }

    // Орбита .scf
    std::vector<SpaceStation> orbit;
    const char* scf = "example/RA140423_1200_v02.scf";
    if (!read_scf_orbit(scf, orbit)) { printf("FAIL: не прочитать %s\n", scf); return 1; }
    printf("\nОрбита: %zu точек\n", orbit.size());
    if (!orbit.empty()) {
        const auto& p0 = orbit.front();
        printf("  первая: MJD %d utc=%.6f  X=%.6f Y=%.6f Z=%.6f км  V=(%.6f,%.6f,%.6f) км/с\n",
               p0.mjd, p0.utc, p0.xyz.x(), p0.xyz.y(), p0.xyz.z(), p0.vel.x(), p0.vel.y(), p0.vel.z());
        // Шапка .scf: 2014-04-23T12:00:00 -3898.581953 201293.478887 -48207.641896 ...
        if (p0.mjd != 56770 || std::fabs(p0.xyz.x() - (-3898.581953)) > 1e-4) {
            printf("FAIL: первая точка орбиты не совпала с шапкой .scf\n"); ok = false;
        }
    }

    printf("---------------------------------------------------------------------\n");
    printf("  РЕЗУЛЬТАТ: %s\n", ok ? "PASS" : "FAIL");
    return ok ? 0 : 1;
}
