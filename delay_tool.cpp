// delay_tool.cpp
//
// CLI-инструмент модели задержки: по заданию .cfx и орбите .scf считает полиномы
// задержки всех станций сеанса и пишет файлы .TXT в выходной каталог.
//
// Использование (из любого каталога — эфемериды и EOP ищутся автоматически в проекте):
//   delay_tool <cfx> <scf> <out_dir> [block_sec=60] [degree=5]
//
// Пример:
//   ./delay_tool example/GVLBI_RAKS01A0_L_20140423T130000_ASC_V1.cfx \
//                example/RA140423_1200_v02.scf  out_poly
//
// Каталоги/эфемериды лежат в проекте (external/) и находятся автоматически: инструмент
// ищет их относительно текущего каталога и нескольких родительских уровней.

#ifdef _WIN32
#define WIN32_LEAN_AND_MEAN
#define NOMINMAX
#include <windows.h>
#endif
#include "src/functions.h"
#include <cstdio>
#include <fstream>
#include <string>
#include <vector>
#include <cctype>
#include <filesystem>

using namespace ariadna;

// Поиск файла проекта: относительно cwd и до 4 родительских уровней ("", "../", "../../"...).
static std::string find_project_file(const std::string& rel) {
    std::string prefix;
    for (int up = 0; up <= 4; ++up) {
        std::string path = prefix + rel;
        std::ifstream f(path, std::ios::binary);
        if (f.good()) return path;
        prefix += "../";
    }
    return rel; // не найдено — вернём как есть (будет диагностика ниже)
}

int main(int argc, char** argv) {
#ifdef _WIN32
    SetConsoleOutputCP(CP_UTF8);
#endif
    std::vector<std::string> a(argv + 1, argv + argc);
    if (a.size() < 2) {
        std::printf("Использование: %s <cfx> [scf] <out_dir> [block_sec=60] [degree=5] [tropo=1] [recv=PUSHCH22]\n", argv[0]);
        std::printf("  <cfx>     файл задания коррелятора (.cfx)\n");
        std::printf("  [scf]     файл орбиты космической станции (.scf). НЕОБЯЗАТЕЛЕН: если не указан\n");
        std::printf("            (или \"-\"), орбита берётся из cfx (ORB_FILE). Указывайте, чтобы\n");
        std::printf("            подставить свой файл (напр. при проверке на своём компьютере).\n");
        std::printf("  <out_dir> каталог для выходных полиномов (.TXT)\n");
        std::printf("  [tropo]   1 = с тропосферой (по умолчанию), 0 = только геометрия в вакууме\n");
        std::printf("  [recv]    имя пункта приёма в ITRF2005_2.CAT (по умолч. PUSHCH22; напр. GBT_VLBA)\n");
        std::printf("            при наличии космоса пишется <cfx>_p.cfx с пересчитанными TIMEOFS\n");
        return 1;
    }
    // <scf> опционален. Слот присутствует, если 2-й аргумент оканчивается на .scf или равен "-".
    // "-" или отсутствие -> орбиту берёт process_task из cfx (ORB_FILE); иначе -> из этого файла.
    auto is_scf_slot = [](const std::string& s) {
        if (s == "-") return true;
        if (s.size() < 4) return false;
        std::string e = s.substr(s.size() - 4);
        for (char& c : e) c = (char)std::tolower((unsigned char)c);
        return e == ".scf";
    };
    std::string cfx = a[0];
    std::string scf; size_t p = 1;
    if (a.size() >= 3 && is_scf_slot(a[1])) { scf = (a[1] == "-") ? "" : a[1]; p = 2; }
    if (p >= a.size()) { std::fprintf(stderr, "Не указан <out_dir>.\n"); return 1; }
    std::string outdir = a[p];
    size_t o = p + 1; // индекс первого необязательного аргумента (block_sec)
    double block_sec = (a.size() > o + 0) ? std::stod(a[o + 0]) : 60.0;
    int degree       = (a.size() > o + 1) ? std::stoi(a[o + 1]) : 5;
    bool with_tropo  = (a.size() > o + 2) ? (std::stoi(a[o + 2]) != 0) : true;
    std::string recv = (a.size() > o + 3) ? a[o + 3] : "PUSHCH22";

    // Эфемериды и каталог EOP — автоматически из проекта.
    std::string eph = find_project_file("external/dephem-master/linux_p1550p2650.440t");
    std::string eop = find_project_file("external/catalogs/EOPC04_14_IAU2000_62-now.cat");

    try { init_ephemeris(eph); }
    catch (const std::exception& e) {
        std::fprintf(stderr, "Не удалось загрузить эфемериды (%s): %s\n", eph.c_str(), e.what());
        std::fprintf(stderr, "Положите файл эфемерид в external/dephem-master/ (см. BUILD.md).\n");
        return 2;
    }

    std::error_code ec; std::filesystem::create_directories(outdir, ec);
    std::printf("Задание: %s\nОрбита:  %s\nВыход:   %s\nБлок=%.0f с, степень=%d, тропосфера=%s\n---\n",
                cfx.c_str(), scf.empty() ? "из cfx (ORB_FILE)" : scf.c_str(), outdir.c_str(), block_sec, degree,
                with_tropo ? "вкл" : "выкл (только геометрия в вакууме)");

    process_task(cfx, scf, outdir, eop, block_sec, degree, 6.0, with_tropo, recv);

    std::printf("---\nГотово. Полиномы записаны в %s\n", outdir.c_str());
    return 0;
}
