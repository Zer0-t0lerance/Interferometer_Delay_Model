// test_delay_poly.cpp
//
// Проверка слоя полиномов: считаем полиномы задержки BADARY по всему сеансу (13:00-14:00)
// и сверяем P0 (и P1) КАЖДОГО 1-мин блока с эталоном example/BADARY.TXT. Ожидаем P0 в
// пределах ~0.3 м от эталона (как в одноточечной сверке test_geodelay); полиномы сшиты.

#ifdef _WIN32
#define WIN32_LEAN_AND_MEAN
#define NOMINMAX
#include <windows.h>
#endif
#include "../src/functions.h"
#include "../src/catalog_bridge.h"
#include "../src/READ_CAT.h"
#include <cstdio>
#include <cmath>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>

using namespace ariadna;

static std::string find_eph() {
    const char* c[] = {"external/dephem-master/linux_p1550p2650.440t", "../external/dephem-master/linux_p1550p2650.440t"};
    for (auto p : c) { std::ifstream f(p, std::ios::binary); if (f.good()) return p; } return c[0];
}

// Считать P0 каждого блока из эталонного .TXT.
static std::vector<double> read_ref_P0(const std::string& path) {
    std::vector<double> p0; std::ifstream f(path); std::string line;
    while (std::getline(f, line)) {
        if (line.rfind("P0", 0) == 0) { size_t e = line.find('='); if (e != std::string::npos) p0.push_back(std::stod(line.substr(e + 1))); }
    }
    return p0;
}

int main() {
#ifdef _WIN32
    SetConsoleOutputCP(CP_UTF8);
#endif
    printf("=====================================================================\n");
    printf("  Полиномы задержки BADARY по сеансу vs эталон BADARY.TXT (P0 блоков)\n");
    printf("---------------------------------------------------------------------\n");
    try { init_ephemeris(find_eph()); } catch (const std::exception& e) { printf("SKIP: %s\n", e.what()); return 2; }

    CfxTask task;
    if (!parse_cfx("example/GVLBI_RAKS01A0_L_20140423T130000_ASC_V1.cfx", task)) { printf("FAIL: cfx\n"); return 1; }
    // BADARY и источник из задания.
    const CfxStation* bad = nullptr; for (auto& s : task.stations) if (s.name == "BADARY") bad = &s;
    if (!bad) { printf("FAIL: нет BADARY\n"); return 1; }
    const CfxSource& src = task.sources[0];

    // Сеанс: от первого скана до конца последнего (13:00..14:00).
    int mjd0 = task.scans.front().mjd; double utc0 = task.scans.front().utc;
    double t_end = 0; for (auto& sc : task.scans) { double e = (sc.mjd - mjd0) * 86400.0 + sc.utc * 86400.0 + sc.dur_sec; if (e > t_end) t_end = e; }
    double dur = t_end - utc0 * 86400.0;

    // EOP (56767..56773), leap 2014 = 35.
    const double rows[7][4] = {
        {56767, 0.060944, 0.436953, -0.2295574}, {56768, 0.062152, 0.437984, -0.2308630},
        {56769, 0.063753, 0.438706, -0.2322999}, {56770, 0.065273, 0.439403, -0.2338988},
        {56771, 0.066868, 0.439919, -0.2356298}, {56772, 0.068557, 0.440457, -0.2374814},
        {56773, 0.070301, 0.441033, -0.2394680}};
    std::vector<EOPData> eop(7);
    for (int i = 0; i < 7; ++i) { eop[i].mjd = rows[i][0]; eop[i].x = rows[i][1]; eop[i].y = rows[i][2]; eop[i].ut1_utc = rows[i][3]; eop[i].ut1_tai = rows[i][3] - 35.0; }

    printf("  сеанс: MJD %d utc0=%.6f dur=%.0f с; блок 60 с, степень 5, сетка 6 с\n", mjd0, utc0, dur);
    StationPoly poly = compute_station_poly(*bad, src, mjd0, utc0, dur, eop, 60.0, 5, 6.0);
    write_station_poly("example/BADARY_OUR.TXT", poly);

    std::vector<double> ref = read_ref_P0("example/BADARY.TXT");
    printf("  блоков: наши=%zu эталон=%zu\n", poly.blocks.size(), ref.size());

    double maxerr = 0; int n = static_cast<int>(std::min(poly.blocks.size(), ref.size()));
    for (int i = 0; i < n; ++i) {
        double our = poly.blocks[i].coef[0], r = ref[i], e = std::fabs(our - r);
        if (e > maxerr) maxerr = e;
        if (i < 3 || i == n - 1)
            printf("  блок %2d: P0 наш=% .12e эталон=% .12e |Δ|=%.2e с (%.2f м)\n", i, our, r, e, e * 3e8);
    }
    printf("  МАКС |ΔP0| по всем блокам = %.3e с (%.2f м)\n", maxerr, maxerr * 3e8);
    bool ok = (n > 0) && (maxerr < 5e-9); // ~1.5 м допуск
    printf("---------------------------------------------------------------------\n");
    printf("  РЕЗУЛЬТАТ: %s\n", ok ? "PASS" : "FAIL");
    return ok ? 0 : 1;
}
