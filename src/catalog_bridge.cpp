#include "catalog_bridge.h"
#include "functions.h"
#include <iostream>
#include <string>
#include <cstring>
#include <cmath>

namespace ariadna {

// Имена в каталогах — поля фиксированной ширины (strncpy без нуль-терминатора).
// Берём первый токен (до пробела), не длиннее maxlen.
static std::string trim_name(const char* s, size_t maxlen) {
    std::string r;
    for (size_t i = 0; i < maxlen; ++i) {
        char c = s[i];
        if (c == '\0' || c == ' ' || c == '\t' || c == '\r' || c == '\n') break;
        r.push_back(c);
    }
    return r;
}

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

void build_stations_from_vtrf(const std::vector<ant_vtrf_record>& raw_vtrf,
                              const std::vector<std::string>& names,
                              std::vector<Station>& stations) {
    stations.clear();
    for (const std::string& nm : names) {
        for (const ant_vtrf_record& rec : raw_vtrf) {
            if (nm == trim_name(rec.ant_name, 9)) {
                Station st;
                st.name = nm;
                st.xyz << rec.x, rec.y, rec.z;       // ITRF на эпоху 2000.0 [м]
                st.vel << rec.vx, rec.vy, rec.vz;    // скорость [м/год]
                st.domes = rec.domes_nb;
                stations.push_back(st);
                break;
            }
        }
    }
}

void build_sources_from_icrf(const std::vector<src_icrf_record>& raw_icrf,
                             const std::vector<std::string>& names,
                             std::vector<Source>& sources) {
    sources.clear();
    for (const std::string& nm : names) {
        for (const src_icrf_record& rec : raw_icrf) {
            if (nm == trim_name(rec.src_name_j2000, 16)) {
                Source s;
                s.name = nm;
                s.icrf_name = rec.src_name_1950;
                // RA (часы) -> радианы; Dec (градусы) -> радианы (знак по dec_d).
                s.ra = (rec.ra_h + rec.ra_m / 60.0 + rec.ra_s / 3600.0) * 15.0 * cnst::CDEGRAD;
                double sign = (rec.dec_d < 0.0 || (rec.dec_d == 0.0 && (rec.dec_m < 0 || rec.dec_s < 0))) ? -1.0 : 1.0;
                s.dec = sign * (std::fabs(rec.dec_d) + std::fabs(rec.dec_m) / 60.0 + std::fabs(rec.dec_s) / 3600.0) * cnst::CDEGRAD;
                s.ra_rate = 0.0; s.dec_rate = 0.0;
                sources.push_back(s);
                break;
            }
        }
    }
}

void build_def_par_from_ant_info(const std::vector<ant_record>& raw_ant,
                                 std::vector<Station>& stations) {
    // Маппинг ant_record (ANTENNA_INFO, Nothnagel) -> DefPar станции + опорное давление.
    // Соответствует READ_CAT42corr def_par: hf=ant_h, gamma_hf=th_exp_f, hp=f_a_len (fixed axis),
    // gamma_hp=th_exp_f_a, AO=a_of_len, hv=dist_mov_a, hs=sub_h.
    for (Station& st : stations) {
        for (const ant_record& rec : raw_ant) {
            if (st.name == trim_name(rec.telescope, 8)) {
                DefPar& d = st.def_par;
                d.name     = st.name;
                d.t_0      = rec.ref_temp;
                d.p_0      = rec.ref_press;
                d.ant_diam = rec.ant_diam;
                d.hf = rec.ant_h;      d.gamma_hf = rec.th_exp_f;
                d.hp = rec.f_a_len;    d.gamma_hp = rec.th_exp_f_a;
                d.AO = rec.a_of_len;   d.gamma_AO = rec.a_of_thermal;
                d.hv = rec.dist_mov_a; d.gamma_hv = rec.therm_exp;
                d.hs = rec.sub_h;      d.gamma_hs = rec.sub_exp;
                st.atm_load.p_0 = rec.ref_press; // опорное давление для SITE_ATM40
                break;
            }
        }
    }
}

void select_eop_nodes(const std::vector<eop_record>& raw_eop, int mjd_obs, int n,
                      std::vector<EOPData>& nodes) {
    nodes.clear();
    // Индекс записи с MJD == mjd_obs (0h дня наблюдения).
    int center = -1;
    for (size_t i = 0; i < raw_eop.size(); ++i) {
        if (static_cast<int>(raw_eop[i].MJD + 0.5) == mjd_obs) { center = static_cast<int>(i); break; }
    }
    if (center < 0) return;
    int half = n / 2;
    int start = center - half;
    if (start < 0) start = 0;
    if (start + n > static_cast<int>(raw_eop.size())) start = static_cast<int>(raw_eop.size()) - n;

    for (int k = start; k < start + n; ++k) {
        const eop_record& r = raw_eop[k];
        double idelt;
        nsec(r.MJD, idelt);              // TAI-UTC (високосные секунды)
        EOPData e{};
        e.mjd     = r.MJD;
        e.ut1_utc = r.ut1_utc;
        e.ut1_tai = r.ut1_utc - idelt;   // UT1-TAI = (UT1-UTC) - (TAI-UTC)
        e.x       = r.x;
        e.y       = r.y;
        e.dpsi    = r.dpsi;              // в EOPC04_IAU2000 это dX (используется как есть)
        e.deps    = r.deps;             // dY
        nodes.push_back(e);
    }
}

} // namespace ariadna