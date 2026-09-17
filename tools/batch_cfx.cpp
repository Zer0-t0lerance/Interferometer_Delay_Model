// batch_cfx.cpp
//
// Пакетный прогон заданий коррелятора: по дереву «корень/<эксперимент>/*.cfx» считает
// полиномы задержки и u,v,w и готовит новые задания *_p.cfx.
//
// Что делает для каждого входного .cfx:
//   1) берёт МЕТКУ ВРЕМЕНИ для расчёта TIMEOFS из САМОГО ЭТОГО ЗАДАНИЯ — из второго поля
//      строки «TIMEOFS<xx> = <задержка>, <метка>». Это делается через штатную точку
//      расширения модели (set_time_mark_hook / set_time_mark_from_cfx), сама модель
//      не меняется;
//   2) кладёт полиномы в <эксперимент>/new_poly (рядом с имеющейся папкой poly);
//   3) переносит готовый <имя>_p.cfx из new_poly в папку эксперимента, рядом с исходным;
//   4) в этом _p.cfx правит строку вывода коррелятора: <...>.uvx -> <...>_p.uvx.
//
// ПОЧЕМУ МЕТКА ИМЕННО ИЗ ЗАДАНИЯ, А НЕ ИЗ ИМЕНИ ФАЙЛА. Метку в строку TIMEOFS записывает
// отдельная программа коррелятора: она читает сами файлы данных и вычисляет фактическое
// начало записи. Оно может не совпадать с временем в имени файла — данные бывает начинаются
// не с первого отсчёта (сбой записи, повреждение файла). Расхождение меньше секунды, но для
// задержки это очень много. Поэтому метка из задания точнее имени файла и берётся первой.
//
// Если строки TIMEOFS для файла данных нет (та программа отработала не до конца), время
// берётся из имени файла данных (кодировка YYYYDDDHHMMSS) с предупреждением. Отключается
// ключом --no-fallback: тогда без метки TIMEOFS для этого файла просто не пишется.
//
// Сборка и запуск (из корня репозитория):
//   sh build.sh tools/batch_cfx.cpp          # соберёт tools/batch_cfx.exe
//   ./tools/batch_cfx.exe Correlator_Tests
//
// Ключи:
//   --no-fallback     без строки TIMEOFS не писать TIMEOFS вовсе (не брать время из имени)
//   --only=<подстрока> обработать только задания, чей путь содержит подстроку
//   --dry-run         только разобрать задания и показать план, ничего не писать

#ifdef _WIN32
#define WIN32_LEAN_AND_MEAN
#define NOMINMAX
#include <windows.h>
#endif
#include "../src/functions.h"
#include <cstdio>
#include <cctype>
#include <filesystem>
#include <fstream>
#include <map>
#include <regex>
#include <string>
#include <vector>

using namespace ariadna;
namespace fs = std::filesystem;

// ------------------------------------------------------------------ метки текущего задания
static std::map<std::string, double> g_marks;   // индекс FILExx -> метка UTC в MJD
static bool g_fallback_name = true;             // нет метки -> брать время из имени файла
static int  g_used = 0, g_missing = 0;

// Обработчик метки времени: отдаёт модели метку из строки TIMEOFS текущего задания.
static bool mark_from_task(const TimeMarkRequest& req, double& mjd_utc) {
    auto it = g_marks.find(req.file_index);
    if (it != g_marks.end()) { mjd_utc = it->second; ++g_used; return true; }
    ++g_missing;
    if (g_fallback_name) {
        // Штатный разбор имени файла: YYYYDDDHHMMSS. Переключаем источник на .cfx,
        // спрашиваем модель и возвращаем переключатель обратно.
        set_time_mark_from_cfx(true);
        bool ok = get_time_mark(req, mjd_utc);
        set_time_mark_from_cfx(false);
        std::printf("      ВНИМАНИЕ: FILE%s — строки TIMEOFS в задании нет, время взято "
                    "из имени файла%s.\n      Метку по данным пишет отдельная программа "
                    "коррелятора; здесь она её не записала.\n",
                    req.file_index.c_str(), ok ? "" : " (не разобрано)");
        return ok;
    }
    std::printf("      FILE%s: строки TIMEOFS нет -> TIMEOFS для этого файла не пишется\n",
                req.file_index.c_str());
    return false;
}

// -------------------------------------------------------------------- разбор задания
// Метки берём из КОСМИЧЕСКОГО блока [$TLSC]. Признак блока — строка ORB_FILE, и она стоит
// ПОСЛЕ строк FILExx, поэтому блок разбирается целиком, а не построчно на лету.
static bool read_marks(const fs::path& cfx, std::map<std::string, double>& marks, int& nfiles) {
    std::ifstream in(cfx);
    if (!in) return false;
    std::vector<std::vector<std::string>> blocks;
    std::vector<std::string>* cur = nullptr;
    std::string line;
    while (std::getline(in, line)) {
        while (!line.empty() && (line.back() == '\r' || line.back() == '\n')) line.pop_back();
        size_t a = line.find_first_not_of(" \t");
        std::string t = (a == std::string::npos) ? "" : line.substr(a);
        if (t.rfind("[$TLSC]", 0) == 0) { blocks.emplace_back(); cur = &blocks.back(); continue; }
        if (t.rfind("[$", 0) == 0) { cur = nullptr; continue; }
        if (cur) cur->push_back(t);
    }
    static const std::regex re_file(R"(FILE(\d+)\s*=\s*(.*))");
    static const std::regex re_tofs(R"(TIMEOFS(\d+)\s*=\s*[^,]+,\s*([-\d.eE+]+))");
    for (const auto& b : blocks) {
        bool space = false;
        for (const auto& t : b) if (t.rfind("ORB_FILE", 0) == 0) { space = true; break; }
        if (!space) continue;
        nfiles = 0;
        std::smatch m;
        for (const auto& t : b) {
            if (std::regex_match(t, m, re_file)) ++nfiles;
            else if (std::regex_search(t, m, re_tofs)) marks[m[1].str()] = std::stod(m[2].str());
        }
        return true;
    }
    return false;
}

// ------------------------------------------------- правка имени выходного файла коррелятора
// Строка вида «OUT FILE = %W:ИМЯ.uvx» -> «... ИМЯ_p.uvx». Идемпотентно.
static bool patch_out_file(const fs::path& cfx) {
    std::ifstream in(cfx);
    if (!in) return false;
    std::vector<std::string> lines;
    std::string line;
    bool patched = false;
    while (std::getline(in, line)) {
        while (!line.empty() && (line.back() == '\r' || line.back() == '\n')) line.pop_back();
        size_t a = line.find_first_not_of(" \t");
        std::string t = (a == std::string::npos) ? "" : line.substr(a);
        // «OUT FILE» может быть записано с разным числом пробелов — сравниваем без них
        std::string key;
        for (char c : t.substr(0, 10)) if (!std::isspace((unsigned char)c)) key += (char)std::toupper((unsigned char)c);
        if (key.rfind("OUTFILE", 0) == 0) {
            size_t dot = line.rfind(".uvx");
            if (dot == std::string::npos) dot = line.rfind(".UVX");
            if (dot != std::string::npos) {
                std::string stem = line.substr(0, dot);
                if (stem.size() < 2 || stem.compare(stem.size() - 2, 2, "_p") != 0) {
                    line = stem + "_p" + line.substr(dot);
                    patched = true;
                }
            }
        }
        lines.push_back(line);
    }
    in.close();
    std::ofstream out(cfx);
    if (!out) return false;
    for (const auto& l : lines) out << l << "\n";
    return patched;
}

int main(int argc, char** argv) {
#ifdef _WIN32
    SetConsoleOutputCP(CP_UTF8);
#endif
    std::vector<std::string> a(argv + 1, argv + argc);
    bool dry = false;
    std::string root, only;
    for (const auto& s : a) {
        if (s == "--no-fallback") g_fallback_name = false;
        else if (s.rfind("--only=", 0) == 0) only = s.substr(7);
        else if (s == "--dry-run") dry = true;
        else if (s.rfind("--", 0) == 0) { std::fprintf(stderr, "Неизвестный ключ: %s\n", s.c_str()); return 1; }
        else root = s;
    }
    if (root.empty()) {
        std::printf("Использование: %s <корень> [--no-fallback] [--only=<подстрока>] [--dry-run]\n", argv[0]);
        std::printf("  <корень>           папка с подпапками экспериментов (напр. Correlator_Tests)\n");
        std::printf("  По умолчанию метка берётся из строки TIMEOFS задания (её пишет отдельная\n");
        std::printf("  программа коррелятора по самим данным — она точнее имени файла). Если строки\n");
        std::printf("  нет, время берётся из имени файла данных с предупреждением.\n");
        std::printf("  --no-fallback      без строки TIMEOFS не писать TIMEOFS вовсе\n");
        std::printf("  --only=<подстрока> обработать только задания, чей путь содержит подстроку\n");
        std::printf("  --dry-run          показать план, ничего не писать\n");
        std::printf("\nДля каждого <корень>/<эксперимент>/*.cfx: полиномы -> <эксперимент>/new_poly,\n");
        std::printf("новое задание -> <эксперимент>/<имя>_p.cfx, вывод коррелятора -> *_p.uvx.\n");
        return 0;
    }
    if (!fs::is_directory(root)) { std::fprintf(stderr, "Нет такой папки: %s\n", root.c_str()); return 1; }

    const std::string eph = "external/dephem-master/linux_p1550p2650.440t";
    const std::string eop = "external/catalogs/EOPC04_14_IAU2000_62-now.cat";
    if (!dry) {
        try { init_ephemeris(eph); }
        catch (const std::exception& e) {
            std::fprintf(stderr, "Не удалось загрузить эфемериды (%s): %s\n", eph.c_str(), e.what());
            std::fprintf(stderr, "Запускать из корня репозитория (см. BUILD.md).\n");
            return 2;
        }
        set_time_mark_hook(&mark_from_task);
        set_time_mark_from_cfx(false);   // метку даёт обработчик, а не имя файла
    }

    int tasks = 0, done = 0, failed = 0, files_tot = 0, marks_tot = 0, no_mark = 0, uvx = 0;
    std::vector<std::string> gaps;

    std::vector<fs::path> exps;
    for (const auto& e : fs::directory_iterator(root))
        if (e.is_directory() && e.path().filename() != "new") exps.push_back(e.path());
    std::sort(exps.begin(), exps.end());

    for (const auto& exp : exps) {
        std::vector<fs::path> cfxs, scfs;
        for (const auto& f : fs::directory_iterator(exp)) {
            if (!f.is_regular_file()) continue;
            std::string n = f.path().filename().string();
            std::string ext = f.path().extension().string();
            for (char& c : ext) c = (char)std::tolower((unsigned char)c);
            if (ext == ".cfx" && n.size() > 6 && n.compare(n.size() - 6, 6, "_p.cfx") != 0) {
                if (only.empty() || f.path().string().find(only) != std::string::npos)
                    cfxs.push_back(f.path());
            } else if (ext == ".scf") scfs.push_back(f.path());
        }
        if (cfxs.empty()) continue;
        std::sort(cfxs.begin(), cfxs.end());
        std::string scf = scfs.empty() ? "" : scfs.front().string();
        std::printf("\n=== %s === (заданий: %zu, орбита: %s)\n", exp.filename().string().c_str(),
                    cfxs.size(), scf.empty() ? "из cfx (ORB_FILE)" : scfs.front().filename().string().c_str());

        fs::path outdir = exp / "new_poly";
        for (const auto& cfx : cfxs) {
            ++tasks;
            g_marks.clear(); g_used = g_missing = 0;
            int nfiles = 0;
            if (!read_marks(cfx, g_marks, nfiles)) {
                std::printf("  %s: космический блок не найден — пропуск\n", cfx.filename().string().c_str());
                ++failed; continue;
            }
            files_tot += nfiles; marks_tot += (int)g_marks.size();
            if ((int)g_marks.size() < nfiles) {
                no_mark += nfiles - (int)g_marks.size();
                gaps.push_back(exp.filename().string() + "/" + cfx.filename().string() +
                               " — файлов " + std::to_string(nfiles) + ", меток " + std::to_string(g_marks.size()));
            }
            std::printf("  %s: файлов данных %d, меток в задании %zu\n",
                        cfx.filename().string().c_str(), nfiles, g_marks.size());
            if (dry) continue;

            std::error_code ec; fs::create_directories(outdir, ec);
            try {
                process_task(cfx.string(), scf, outdir.string(), eop, 60.0, 5, 6.0, true, "auto");
            } catch (const std::exception& e) {
                std::fprintf(stderr, "  ОШИБКА на %s: %s\n", cfx.filename().string().c_str(), e.what());
                ++failed; continue;
            }
            // Готовый _p.cfx модель кладёт рядом с полиномами — переносим к исходному заданию.
            std::string stem = cfx.filename().string();
            stem = stem.substr(0, stem.size() - 4) + "_p.cfx";
            fs::path made = outdir / stem, dest = exp / stem;
            if (fs::exists(made)) {
                fs::remove(dest, ec);
                fs::rename(made, dest, ec);
                if (ec) { fs::copy_file(made, dest, fs::copy_options::overwrite_existing, ec); fs::remove(made, ec); }
                if (patch_out_file(dest)) ++uvx;
                std::printf("  задание -> %s\n", dest.filename().string().c_str());
                ++done;
            } else {
                std::fprintf(stderr, "  %s: _p.cfx не создан (нет космоса или не найден пункт приёма)\n",
                             cfx.filename().string().c_str());
                ++failed;
            }
        }
    }

    std::printf("\n---------------------------------------------------------------\n");
    std::printf("Заданий обработано: %d из %d (сбоев: %d)\n", done, tasks, failed);
    std::printf("Файлов данных космоса: %d, меток взято из заданий: %d\n", files_tot, marks_tot);
    std::printf("Файлов без метки в задании: %d (%s)\n", no_mark,
                g_fallback_name ? "время взято из имени файла" : "TIMEOFS для них не записан");
    std::printf("Имя вывода коррелятора исправлено на *_p.uvx: в %d заданиях\n", uvx);
    if (!gaps.empty()) {
        std::printf("\nЗадания, где меток меньше, чем файлов данных:\n");
        for (const auto& g : gaps) std::printf("  %s\n", g.c_str());
    }
    return failed ? 1 : 0;
}
