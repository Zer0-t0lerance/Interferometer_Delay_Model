#pragma once
#include <vector>
#include "structures.h"
#include "READ_CAT.h"

namespace ariadna {

/**
 * @brief Переносит данные океанической нагрузки из сырого каталога в структуры станций.
 */
void map_ocean_tides_to_stations(const std::vector<oc_record>& raw_oc_data, std::vector<Station>& stations);

/**
 * @brief Переносит данные атмосферной нагрузки из сырого каталога в структуры станций.
 */
void map_atm_loading_to_stations(const std::vector<atm_record>& raw_atm_data, std::vector<Station>& stations);

/**
 * @brief Строит станции по именам из каталога VTRF (координаты ITRF на эпоху 2000.0 + скорости).
 * @param[in]  raw_vtrf Записи каталога VTRF (ReadVTRF).
 * @param[in]  names    Имена нужных станций (IVS-имена).
 * @param[out] stations Заполненные Station (name, xyz [м], vel [м/год]); ненайденные пропускаются.
 */
void build_stations_from_vtrf(const std::vector<ant_vtrf_record>& raw_vtrf,
                              const std::vector<std::string>& names,
                              std::vector<Station>& stations);

/**
 * @brief Строит источники по именам из каталога ICRF (RA h/m/s, Dec d/m/s -> радианы).
 * @param[in]  raw_icrf Записи каталога ICRF (ReadICRF).
 * @param[in]  names    Имена нужных источников (J2000).
 * @param[out] sources  Заполненные Source (name, ra, dec [рад]); ненайденные пропускаются.
 */
void build_sources_from_icrf(const std::vector<src_icrf_record>& raw_icrf,
                             const std::vector<std::string>& names,
                             std::vector<Source>& sources);

/**
 * @brief Выбирает N узлов EOP вокруг даты наблюдения и переводит в EOPData (для interp_iers/interp_eop).
 * @param[in]  raw_eop  Записи каталога EOP (ReadEOP), отсортированы по MJD.
 * @param[in]  mjd_obs  MJD наблюдения.
 * @param[in]  n        Число узлов (обычно 7).
 * @param[out] nodes    N узлов EOPData вокруг mjd_obs (mjd, ut1_utc, ut1_tai, x, y, dpsi, deps).
 */
void select_eop_nodes(const std::vector<eop_record>& raw_eop, int mjd_obs, int n,
                      std::vector<EOPData>& nodes);

} // namespace ariadna