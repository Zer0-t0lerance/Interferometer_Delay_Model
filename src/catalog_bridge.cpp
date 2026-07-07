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

void map_atm_loading_to_stations(const std::vector<atm_record>& raw_atm_data,
                                 std::vector<Station>& stations)
{
    // Раскладка atm_record.coeff[30] (по 10 значений на строку каталога V/E/N):
    //   строка V -> coeff[0..8], строка E -> coeff[10..18], строка N -> coeff[20..28];
    //   в каждой строке 8 модельных коэффициентов A1 B1 A2 B2 A3 B3 b0 b1, 9-й = Sigma (не используется).
    for (Station& st : stations) {
        st.atm_load.coef.setZero();
        st.atm_load.has_data = false;
        for (const atm_record& rec : raw_atm_data) {
            if (st.name == std::string(rec.telescope)) {
                for (int comp = 0; comp < 3; ++comp) {       // 0=V, 1=E, 2=N
                    for (int k = 0; k < 8; ++k) {            // A1 B1 A2 B2 A3 B3 b0 b1
                        st.atm_load.coef(comp, k) = rec.coeff[comp * 10 + k];
                    }
                }
                st.atm_load.has_data = true;
                break;
            }
        }
        // p_0 (опорное давление) заполняется отдельно из каталога ANTENNA_INFO.
    }
}

} // namespace ariadna