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
    if (argc < 4) {
        std::printf("Использование: %s <cfx> <scf> <out_dir> [block_sec=60] [degree=5]\n", argv[0]);
        std::printf("  <cfx>     файл задания коррелятора (.cfx)\n");
        std::printf("  <scf>     файл орбиты космической станции (.scf); \"-\" если нет космоса\n");
        std::printf("  <out_dir> каталог для выходных полиномов (.TXT)\n");
        return 1;
    }
    std::string cfx = argv[1];
    std::string scf = (std::string(argv[2]) == "-") ? "" : argv[2];
    std::string outdir = argv[3];
    double block_sec = (argc > 4) ? std::stod(argv[4]) : 60.0;
    int degree = (argc > 5) ? std::stoi(argv[5]) : 5;

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
    std::printf("Задание: %s\nОрбита:  %s\nВыход:   %s\nБлок=%.0f с, степень=%d\n---\n",
                cfx.c_str(), scf.empty() ? "(нет)" : scf.c_str(), outdir.c_str(), block_sec, degree);

    process_task(cfx, scf, outdir, eop, block_sec, degree, 6.0);

    std::printf("---\nГотово. Полиномы записаны в %s\n", outdir.c_str());
    return 0;
}
