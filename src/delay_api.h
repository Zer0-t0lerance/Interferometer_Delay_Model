// delay_api.h
//
// ОТДЕЛЬНЫЙ модуль-адаптер для встраивания (в т.ч. из Python): считает полиномы задержки и
// координат (u,v,w) ВСЕХ станций сеанса В ПАМЯТИ, без записи файлов. Не изменяет существующий
// код модели — только оркеструет уже публичные функции (parse_cfx, read_scf_orbit, EOP-узлы,
// каталоги нагрузок, compute_station_poly). Логика расчёта — та же, что в process_task.

#pragma once
#include <string>
#include <vector>
#include "structures.h"

namespace ariadna {

// Результат: полиномы всех станций сеанса в памяти. Векторы delay и uvw параллельны
// (i-я станция), имя станции — в delay[i].telescope / uvw[i].telescope.
struct TaskPolys {
    std::vector<StationPoly> delay; // полиномы задержки по станциям (порядок cfx)
    std::vector<StationUvw>  uvw;   // полиномы координат u,v,w по станциям
};

/**
 * @brief Посчитать полиномы задержки и координат (u,v,w) ВСЕХ станций сеанса В ПАМЯТИ.
 *
 * Делает то же, что process_task (посканно, свой источник у каждого скана, участие станций,
 * ход бортовых часов для космоса), но НЕ пишет файлы — возвращает результат в out. Эфемериды
 * должны быть инициализированы заранее (init_ephemeris).
 *
 * @param[in]  cfx_path   Путь к файлу задания .cfx.
 * @param[in]  orbit_path Путь к файлу орбиты .scf; "" — взять из cfx (ORB_FILE).
 * @param[in]  eop_path   Путь к каталогу EOP (EOPC04); рядом ищутся каталоги нагрузок/антенн.
 * @param[in]  block_sec  Длительность блока полинома [с] (эталон = 60).
 * @param[in]  degree     Степень полинома (эталон = 5 -> order = 6).
 * @param[in]  sample_sec Шаг сетки расчёта [с] (плотнее блока; типично 6).
 * @param[in]  with_tropo Подключать тропосферу (true) или только геометрия в вакууме (false).
 * @param[out] out        Полиномы всех станций (delay + uvw).
 * @return true при успехе; false при ошибке разбора cfx / отсутствии сканов / нехватке EOP.
 */
bool compute_task_polys(const std::string& cfx_path, const std::string& orbit_path,
                        const std::string& eop_path, double block_sec, int degree,
                        double sample_sec, bool with_tropo, TaskPolys& out);

} // namespace ariadna
