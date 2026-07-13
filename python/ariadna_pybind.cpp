// ariadna_pybind.cpp
//
// pybind11-обёртка для использования модели задержки из Python. Экспортирует адаптер
// compute_task_polys (полиномы всех станций в память) и структуры результата, а также
// init_ephemeris. НЕ изменяет модель — только связывает публичный API с Python.
//
// Сборка: см. python/build_pybind.sh / python/build_pybind.ps1. Пример: python/example.py.

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>       // авто-конвертация std::vector <-> list
#include "../src/functions.h"
#include "../src/delay_api.h"

namespace py = pybind11;
using namespace ariadna;

PYBIND11_MODULE(ariadna, m) {
    m.doc() = "Модель задержки VLBI (перенос ARIADNA): полиномы задержки и координат (u,v,w).";

    // Блок полинома задержки: коэффициенты P0..Pn (по возрастанию степени), в секундах.
    py::class_<DelayPolyBlock>(m, "DelayBlock")
        .def_readonly("mjd", &DelayPolyBlock::mjd)
        .def_readonly("utc_start", &DelayPolyBlock::utc_start)
        .def_readonly("utc_stop", &DelayPolyBlock::utc_stop)
        .def_readonly("source", &DelayPolyBlock::source)
        .def_readonly("coef", &DelayPolyBlock::coef);   // список float [P0..Pn]

    py::class_<StationPoly>(m, "StationPoly")
        .def_readonly("telescope", &StationPoly::telescope)
        .def_readonly("source", &StationPoly::source)
        .def_readonly("order", &StationPoly::order)
        .def_readonly("blocks", &StationPoly::blocks);  // список DelayBlock

    // Блок полинома координат: коэффициенты по осям u, v, w (в секундах).
    py::class_<UvwPolyBlock>(m, "UvwBlock")
        .def_readonly("mjd", &UvwPolyBlock::mjd)
        .def_readonly("utc_start", &UvwPolyBlock::utc_start)
        .def_readonly("utc_stop", &UvwPolyBlock::utc_stop)
        .def_readonly("source", &UvwPolyBlock::source)
        .def_readonly("u", &UvwPolyBlock::u)
        .def_readonly("v", &UvwPolyBlock::v)
        .def_readonly("w", &UvwPolyBlock::w);

    py::class_<StationUvw>(m, "StationUvw")
        .def_readonly("telescope", &StationUvw::telescope)
        .def_readonly("order", &StationUvw::order)
        .def_readonly("blocks", &StationUvw::blocks);   // список UvwBlock

    py::class_<TaskPolys>(m, "TaskPolys")
        .def_readonly("delay", &TaskPolys::delay)        // список StationPoly
        .def_readonly("uvw", &TaskPolys::uvw);           // список StationUvw

    m.def("init_ephemeris", &init_ephemeris, py::arg("eph_file"),
          "Загрузить эфемериды JPL (один раз за процесс, до compute_task_polys).");

    m.def("compute_task_polys",
          [](const std::string& cfx, const std::string& orbit, const std::string& eop,
             double block, int degree, double sample, bool tropo) {
              TaskPolys out;
              bool ok = compute_task_polys(cfx, orbit, eop, block, degree, sample, tropo, out);
              if (!ok) throw std::runtime_error("compute_task_polys: ошибка (см. stderr)");
              return out;
          },
          py::arg("cfx"), py::arg("orbit") = "", py::arg("eop"),
          py::arg("block") = 60.0, py::arg("degree") = 5, py::arg("sample") = 6.0,
          py::arg("tropo") = true,
          "Полиномы задержки и координат (u,v,w) всех станций сеанса В ПАМЯТИ.\n"
          "orbit=\"\" -> орбита берётся из cfx (ORB_FILE). Возвращает TaskPolys.");
}
