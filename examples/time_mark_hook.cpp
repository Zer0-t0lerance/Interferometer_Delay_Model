// time_mark_hook.cpp
//
// Пример: подача СВОИХ меток времени для расчёта TIMEOFS.
//
// TIMEOFS — задержка сброса сигнала с космического телескопа на наземный пункт приёма.
// Она считается на определённый момент времени; этот момент здесь называется МЕТКОЙ ВРЕМЕНИ.
// По умолчанию метка читается из задания .cfx (строка FILExx космической станции). Если
// метки надо брать откуда-то ещё — из своей таблицы, журнала сеанса, базы данных, — их
// отдаёт ваш обработчик.
//
// Три шага:
//   1) написать функцию с сигнатурой TimeMarkHook;
//   2) set_time_mark_hook(&ваша_функция);
//   3) set_time_mark_from_cfx(false)  — переключить источник на обработчик.
// Дальше process_task вызывается как обычно.
//
// Сборка (из корня репозитория):
//   Windows PowerShell:  powershell -ExecutionPolicy Bypass -File .\build.ps1 examples\time_mark_hook.cpp
//   Linux / Git Bash:    sh build.sh examples/time_mark_hook.cpp
// Запуск — из корня репозитория (нужны каталог example/ и эфемериды в external/).

#ifdef _WIN32
#define WIN32_LEAN_AND_MEAN
#define NOMINMAX
#include <windows.h>
#endif
#include "../src/functions.h"
#include <cstdio>
#include <map>
#include <string>

using namespace ariadna;

// ------------------------------------------------------------------------------------
// Ваш источник меток. Здесь для примера — таблица «индекс файла -> метка UTC в MJD».
// В рабочем коде это может быть чтение своего файла, запрос в базу, разбор журнала.
// Таблицу удобно заполнить ОДИН РАЗ до вызова process_task: обработчик дёргается на
// каждую строку FILExx и не должен лезть в медленные источники.
// ------------------------------------------------------------------------------------
static std::map<std::string, double> g_marks;

// Календарная дата UTC -> MJD (сутки). Григорианский календарь.
// jdn — юлианский день для ПОЛУДНЯ даты, поэтому MJD на 00:00 этой даты = jdn - 2400001.
static double ymd_hms_to_mjd(int Y, int M, int D, int hh, int mm, double ss) {
    int a = (14 - M) / 12, y = Y + 4800 - a, m = M + 12 * a - 3;
    long jdn = D + (153 * m + 2) / 5 + 365L * y + y / 4 - y / 100 + y / 400 - 32045;
    return (jdn - 2400001L) + (hh * 3600.0 + mm * 60.0 + ss) / 86400.0;
}

// ------------------------------------------------------------------------------------
// Обработчик метки времени. Вызывается моделью на каждую строку FILExx космической
// станции при записи <cfx>_p.cfx.
//
//   req.cfx_path   — путь к текущему заданию;
//   req.station    — имя космической станции в задании (напр. RASTRON);
//   req.file_index — индекс строки: "00", "01", ...;
//   req.file_value — значение строки FILExx как оно записано в задании.
//
// Вернуть true и заполнить mjd_utc — метка принята. Вернуть false — метки нет, TIMEOFS
// для этого файла данных не будет записан (остальные файлы обрабатываются как обычно).
// ------------------------------------------------------------------------------------
static bool my_time_mark(const TimeMarkRequest& req, double& mjd_utc) {
    auto it = g_marks.find(req.file_index);
    if (it == g_marks.end()) {
        std::printf("  [метка] %s FILE%s — своей метки нет, TIMEOFS пропущен\n",
                    req.station.c_str(), req.file_index.c_str());
        return false;
    }
    mjd_utc = it->second;
    std::printf("  [метка] %s FILE%s -> MJD %.9f  (из %s)\n",
                req.station.c_str(), req.file_index.c_str(), mjd_utc, req.file_value.c_str());
    return true;
}

int main() {
#ifdef _WIN32
    SetConsoleOutputCP(CP_UTF8);
#endif
    // Метки для примера: сеанс 23 апреля 2014, четыре файла по 15 минут с 13:00 UTC.
    g_marks["00"] = ymd_hms_to_mjd(2014, 4, 23, 13,  0, 0.0);
    g_marks["01"] = ymd_hms_to_mjd(2014, 4, 23, 13, 15, 0.0);
    g_marks["02"] = ymd_hms_to_mjd(2014, 4, 23, 13, 30, 0.0);
    g_marks["03"] = ymd_hms_to_mjd(2014, 4, 23, 13, 45, 0.0);

    init_ephemeris("external/dephem-master/linux_p1550p2650.440t");

    // Шаги 2 и 3: поставить обработчик и переключить на него источник метки.
    set_time_mark_hook(&my_time_mark);
    set_time_mark_from_cfx(false);
    std::printf("источник метки времени: %s\n\n",
                time_mark_from_cfx() ? "задание .cfx" : "обработчик");

    process_task("example/GVLBI_RAKS01A0_L_20140423T130000_ASC_V1.cfx",
                 "example/RA140423_1200_v02.scf",
                 "out_poly_hook",
                 "external/catalogs/EOPC04_14_IAU2000_62-now.cat");

    // Вернуть штатный источник (полезно, если в одном процессе считается несколько заданий
    // и часть из них должна брать метку из своего .cfx).
    set_time_mark_from_cfx(true);
    std::printf("\nГотово. TIMEOFS в out_poly_hook/..._p.cfx посчитаны по вашим меткам.\n");
    return 0;
}
