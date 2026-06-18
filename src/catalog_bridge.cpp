#include "catalog_bridge.h"
#include <iostream>
#include <string>

namespace ariadna {

void map_ocean_tides_to_stations(const std::vector<oc_record>& raw_oc_data, 
                                 std::vector<Station>& stations) 
{
    for (Station& st : stations) {
        bool found = false;
        for (const oc_record& rec : raw_oc_data) {
            // Сравниваем имена станций
            if (st.name == std::string(rec.telescope)) {
                for (int axis = 0; axis < 3; ++axis) {
                    for (int wave = 0; wave < 11; ++wave) {
                        int flat_index = axis * 11 + wave;
                        st.tide_data.amplitudes(axis, wave) = rec.coeff1[flat_index];
                        st.tide_data.phases(axis, wave)     = rec.coeff2[flat_index];
                    }
                }
                found = true;
                break;
            }
        }
        if (!found) {
            // Если станция не найдена в каталоге океана, смещений нет
            st.tide_data.amplitudes.setZero();
            st.tide_data.phases.setZero();
        }
    }
}

} // namespace ariadna