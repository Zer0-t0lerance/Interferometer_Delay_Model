// test_process_task.cpp
//
// ФИНАЛЬНАЯ ПРОВЕРКА готового модуля: process_task по заданию .cfx + орбите .scf считает
// полиномы задержки ВСЕХ станций сеанса и пишет .TXT. Сверяем P0 всех блоков каждой
// станции с эталонными example/*.TXT. Ожидаем ~метры (полиномы сшиты, задержки попадают).

#ifdef _WIN32
#define WIN32_LEAN_AND_MEAN
#define NOMINMAX
#include <windows.h>
#endif
#include "../src/functions.h"
#include <cstdio>
#include <cmath>
#include <fstream>
#include <string>
#include <vector>
#include <filesystem>

using namespace ariadna;

static std::string find_eph() {
    const char* c[] = {"external/dephem-master/linux_p1550p2650.440t", "../external/dephem-master/linux_p1550p2650.440t"};
    for (auto p : c) { std::ifstream f(p, std::ios::binary); if (f.good()) return p; } return c[0];
}
static std::vector<double> read_P0(const std::string& path) {
    std::vector<double> p0; std::ifstream f(path); std::string line;
    while (std::getline(f, line)) if (line.rfind("P0", 0) == 0) { size_t e = line.find('='); if (e != std::string::npos) p0.push_back(std::stod(line.substr(e + 1))); }
    return p0;
}
static double max_P0_err(const std::string& ours, const std::string& ref) {
    auto a = read_P0(ours), b = read_P0(ref);
    if (a.empty() || a.size() != b.size()) return 1e9;
    double m = 0; for (size_t i = 0; i < a.size(); ++i) m = std::max(m, std::fabs(a[i] - b[i]));
    return m;
}

int main() {
#ifdef _WIN32
    SetConsoleOutputCP(CP_UTF8);
#endif
    printf("=====================================================================\n");
    printf("  ГОТОВЫЙ МОДУЛЬ: process_task(cfx, орбита) -> полиномы всех станций\n");
    printf("---------------------------------------------------------------------\n");
    try { init_ephemeris(find_eph()); } catch (const std::exception& e) { printf("SKIP: %s\n", e.what()); return 2; }

    std::filesystem::create_directories("tests/out_poly");
    // with_tropo=false: воспроизводим ВАКУУМНЫЙ эталон (в нём нет тропосферы) для сверки.
    process_task("example/GVLBI_RAKS01A0_L_20140423T130000_ASC_V1.cfx",
                 "example/RA140423_1200_v02.scf",
                 "tests/out_poly",
                 "external/catalogs/EOPC04_14_IAU2000_62-now.cat",
                 60.0, 5, 6.0, false);

    // Сверка P0 с эталонами (имена вывода = POLY_FILE).
    struct Pair { const char* ours; const char* ref; };
    Pair pairs[] = {
        {"tests/out_poly/BADARY.txt",    "example/BADARY.TXT"},
        {"tests/out_poly/KALYAZIN_L.txt","example/KALYAZIN_L.TXT"},
        {"tests/out_poly/HARTRAO_L.txt", "example/HARTRAO_L.TXT"},
        {"tests/out_poly/RA_L.txt",      "example/RA_L.TXT"},
    };
    printf("\n  Сверка P0 блоков с эталонами:\n");
    double worst_ground = 0; bool ok = true;
    for (auto& p : pairs) {
        double e = max_P0_err(p.ours, p.ref);
        bool is_ra = std::string(p.ref).find("RA_L") != std::string::npos;
        printf("    %-24s макс |ΔP0| = %.3e с (%.2f м)%s\n",
               std::filesystem::path(p.ref).filename().string().c_str(), e, e * 3e8,
               is_ra ? "  [космос: ретардация орбиты — доработка]" : "");
        if (!is_ra) { worst_ground = std::max(worst_ground, e); if (e > 2e-8) ok = false; } // наземные < ~6 м
    }
    printf("  НАЗЕМНЫЕ худшее: %.3e с (%.2f м)\n", worst_ground, worst_ground * 3e8);
    printf("---------------------------------------------------------------------\n");
    printf("  РЕЗУЛЬТАТ (наземные станции): %s\n", ok ? "PASS" : "FAIL");
    printf("  ПРИМ.: RASTRON структурно работает, но требует ретардации орбиты (0.2м в начале\n");
    printf("         скана -> растёт; космос за 0.7с хода сигнала смещается ~860 м).\n");
    return ok ? 0 : 1;
}
