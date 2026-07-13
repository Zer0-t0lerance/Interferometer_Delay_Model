// use_as_library.cpp
//
// Пример встраивания модели задержки как БИБЛИОТЕКИ в свой C/C++ код.
// Весь публичный API — в src/functions.h (namespace ariadna); структуры — в structures.h,
// константы — в constants.h. Никаких правок модели не нужно: подключаешь заголовок,
// собираешь со всеми модулями src/*.cpp (+ SOFA + dephem) и вызываешь функции.
//
// Сборка (из корня репозитория — build-скрипт собирает ЛЮБОЙ main с полным набором модулей):
//   Windows PowerShell:  powershell -ExecutionPolicy Bypass -File .\build.ps1 examples\use_as_library.cpp
//   Linux / Git Bash:    sh build.sh examples/use_as_library.cpp
// Запуск — из корня репозитория (нужны каталог example/ и эфемериды в external/).

#ifdef _WIN32
#define WIN32_LEAN_AND_MEAN
#define NOMINMAX
#include <windows.h>
#endif
#include "../src/functions.h"
#include <cstdio>

using namespace ariadna;

int main() {
#ifdef _WIN32
    SetConsoleOutputCP(CP_UTF8);
#endif
    // 1) Эфемериды JPL — инициализировать ОДИН раз за процесс.
    init_ephemeris("external/dephem-master/linux_p1550p2650.440t");

    // 2) Высокоуровневый вызов: по заданию .cfx и орбите .scf посчитать полиномы задержки
    //    и координат (u,v,w) ВСЕХ станций и записать файлы (<станция>.txt, <станция>_uvw.txt).
    //    Это ровно то, что делает CLI delay_tool — одна функция.
    process_task("example/GVLBI_RAKS01A0_L_20140423T130000_ASC_V1.cfx",
                 "example/RA140423_1200_v02.scf",          // "" -> орбита из cfx (ORB_FILE)
                 "out_lib",                                 // "" -> каталог из cfx (%W)
                 "external/catalogs/EOPC04_14_IAU2000_62-now.cat",
                 60.0, 5, 6.0, /*with_tropo=*/true, "PUSHCH22");
    std::printf("Готово: полиномы записаны в out_lib/ (файлы *.txt и *_uvw.txt).\n");

    // 3) В ПАМЯТИ, без файлов:
    //    - compute_station_poly(st, src, ...): полиномы ОДНОЙ станции -> StationPoly (задержка)
    //      + StationUvw (координаты). Коэффициенты: poly.blocks[b].coef[k] (задержка),
    //      uvw.blocks[b].u/v/w[k] (координаты). Нужны узлы EOP (select_eop_nodes) и физика
    //      станции из каталогов — как process_task готовит внутри (см. src/delay_poly.cpp).
    //    - process_ariadna(...) / compute_delay_obs(...): сырые ЗНАЧЕНИЯ задержки (не полиномы)
    //      на заданные моменты -> std::vector<DelayResult> / одно наблюдение.
    //    То есть можно получить и полиномы, и сами задержки — по потребности.
    return 0;
}
