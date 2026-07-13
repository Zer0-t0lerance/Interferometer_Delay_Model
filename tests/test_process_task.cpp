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
// Читаем блоки как пары (время старта -> P0). Ключ = строка start (сравнение по времени,
// а не по индексу): посканова модель может пропускать сканы (станция участвует не везде),
// поэтому число блоков у нас и в эталоне может отличаться — сверяем по совпадающим стартам.
static std::vector<std::pair<std::string,double>> read_P0(const std::string& path) {
    std::vector<std::pair<std::string,double>> v; std::ifstream f(path); std::string line, cur;
    while (std::getline(f, line)) {
        if (line.rfind("start", 0) == 0) cur = line;
        else if (line.rfind("P0", 0) == 0) { size_t e = line.find('='); if (e != std::string::npos) v.push_back({cur, std::stod(line.substr(e + 1))}); }
    }
    return v;
}
static double max_P0_err(const std::string& ours, const std::string& ref) {
    auto a = read_P0(ours), b = read_P0(ref);
    if (a.empty()) return 1e9;
    double m = 0; int matched = 0;
    for (const auto& pa : a) for (const auto& pb : b)
        if (pa.first == pb.first) { m = std::max(m, std::fabs(pa.second - pb.second)); ++matched; break; }
    return matched > 0 ? m : 1e9;
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
    double worst = 0; bool ok = true;
    for (auto& p : pairs) {
        double e = max_P0_err(p.ours, p.ref);
        bool is_ra = std::string(p.ref).find("RA_L") != std::string::npos;
        printf("    %-24s макс |ΔP0| = %.3e с (%.2f м)%s\n",
               std::filesystem::path(p.ref).filename().string().c_str(), e, e * 3e8,
               is_ra ? "  [космос: геометрия + ход бортовых часов]" : "");
        worst = std::max(worst, e); if (e > 2e-8) ok = false; // все станции (вкл. RASTRON) < ~6 м
    }
    printf("  ХУДШЕЕ по всем станциям: %.3e с (%.2f м)\n", worst, worst * 3e8);
    printf("---------------------------------------------------------------------\n");
    printf("  РЕЗУЛЬТАТ: %s\n", ok ? "PASS" : "FAIL");
    printf("  ПРИМ.: RASTRON — со 167 м до ~1.2 м: добавлен релятивистский ход бортовых\n");
    printf("         часов относительно наземных (L_G + L_orb, борт выше в гравитац. яме).\n");
    return ok ? 0 : 1;
}
