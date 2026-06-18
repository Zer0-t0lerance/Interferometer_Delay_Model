#pragma once
#include <vector>
#include "structures.h"
#include "READ_CAT.h"

namespace ariadna {

/**
 * @brief Переносит данные океанической нагрузки из сырого каталога в структуры станций.
 */
void map_ocean_tides_to_stations(const std::vector<oc_record>& raw_oc_data, std::vector<Station>& stations);

} // namespace ariadna