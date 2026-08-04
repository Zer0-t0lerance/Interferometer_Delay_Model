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
#include <cstdlib>
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
    if (a.size() < 1) {
        std::printf("Использование: %s <cfx> [scf] [out_dir] [block_sec=60] [degree=5] [tropo=1] [recv=PUSHCH22]\n", argv[0]);
        std::printf("  <cfx>      файл задания коррелятора (.cfx)\n");
        std::printf("  [scf]      файл орбиты (.scf). НЕОБЯЗАТЕЛЕН: если не указан (или \"-\"), орбита\n");
        std::printf("             берётся из cfx (ORB_FILE). Указывайте, чтобы подставить свой файл.\n");
        std::printf("  [out_dir]  каталог для выходных полиномов. НЕОБЯЗАТЕЛЕН: если не указан (или \"-\"),\n");
        std::printf("             берётся из cfx (%%W); имена файлов — из cfx (POLY_FILE).\n");
        std::printf("  [tropo]    1 = с тропосферой (по умолчанию), 0 = только геометрия в вакууме\n");
        std::printf("  [recv]     имя пункта приёма в ITRF2005_2.CAT. По умолчанию auto: определяется\n");
        std::printf("             по префиксу файлов космоса (PUSH->PUSHCH22, GBT->GBT_VLBA). Можно задать явно.\n");
        std::printf("             при наличии космоса пишется <cfx>_p.cfx с пересчитанными TIMEOFS.\n");
        std::printf("  Для каждой станции пишется файл полиномов задержки и файл координат *_uvw.\n");
        return 1;
    }
    // <scf> опционален: слот есть, если аргумент оканчивается на .scf или равен "-".
    auto is_scf_slot = [](const std::string& s) {
        if (s == "-") return true;
        if (s.size() < 4) return false;
        std::string e = s.substr(s.size() - 4);
        for (char& c : e) c = (char)std::tolower((unsigned char)c);
        return e == ".scf";
    };
    // Число (block/degree/tropo) — чтобы отличить их от [out_dir] (нечисловой) при разборе.
    auto is_number = [](const std::string& s) {
        if (s.empty()) return false;
        char* end = nullptr; std::strtod(s.c_str(), &end); return end != s.c_str() && *end == '\0'; };

    std::string cfx = a[0];
    size_t i = 1;
    std::string scf;    // "" = орбита из cfx (ORB_FILE)
    if (i < a.size() && is_scf_slot(a[i])) { scf = (a[i] == "-") ? "" : a[i]; ++i; }
    std::string outdir; // "" = каталог из cfx (%W)
    if (i < a.size() && !is_number(a[i])) { outdir = (a[i] == "-") ? "" : a[i]; ++i; }
    // Безопасный разбор числовых аргументов: при нечисле — понятная ошибка вместо краха stod.
    auto need_num = [&](size_t idx, const char* what) -> std::string {
        if (!is_number(a[idx])) {
            std::fprintf(stderr, "Ошибка: аргумент #%zu '%s' должен быть числом (%s).\n", idx + 2, a[idx].c_str(), what);
            std::fprintf(stderr, "Порядок аргументов: <cfx> [scf.scf] [out_dir] [block] [degree] [tropo] [recv].\n");
            std::fprintf(stderr, "Похоже, указан лишний путь. Пример:\n  delay_tool task.cfx orbit.scf out_dir\n");
            std::exit(1);
        }
        return a[idx];
    };
    double block_sec = (i < a.size()) ? std::stod(need_num(i++, "block_sec")) : 60.0;
    int degree       = (i < a.size()) ? std::stoi(need_num(i++, "degree")) : 5;
    bool with_tropo  = (i < a.size()) ? (std::stoi(need_num(i++, "tropo")) != 0) : true;
    std::string recv = (i < a.size()) ? a[i++] : "auto"; // auto = определить по префиксу файлов (PUSH/GBT)

    // Эфемериды и каталог EOP — автоматически из проекта.
    std::string eph = find_project_file("external/dephem-master/linux_p1550p2650.440t");
    std::string eop = find_project_file("external/catalogs/EOPC04_14_IAU2000_62-now.cat");

    try { init_ephemeris(eph); }
    catch (const std::exception& e) {
        std::fprintf(stderr, "Не удалось загрузить эфемериды (%s): %s\n", eph.c_str(), e.what());
        std::fprintf(stderr, "Положите файл эфемерид в external/dephem-master/ (см. BUILD.md).\n");
        return 2;
    }

    std::error_code ec; if (!outdir.empty()) std::filesystem::create_directories(outdir, ec);
    std::printf("Задание: %s\nОрбита:  %s\nВыход:   %s\nБлок=%.0f с, степень=%d, тропосфера=%s\n---\n",
                cfx.c_str(), scf.empty() ? "из cfx (ORB_FILE)" : scf.c_str(),
                outdir.empty() ? "из cfx (%W)" : outdir.c_str(), block_sec, degree,
                with_tropo ? "вкл" : "выкл (только геометрия в вакууме)");

    process_task(cfx, scf, outdir, eop, block_sec, degree, 6.0, with_tropo, recv);

    std::printf("---\nГотово. Полиномы записаны в %s\n", outdir.c_str());
    return 0;
}
