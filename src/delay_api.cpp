// delay_api.cpp
//
// Реализация адаптера compute_task_polys: полиномы всех станций в память, без файлов.
// Оркеструет уже существующие публичные функции; НИЧЕГО в модели не меняет. Порядок действий
// повторяет process_task (src/delay_poly.cpp): орбита -> узлы EOP -> физика из каталогов ->
// посканный расчёт с учётом участия станций и своего источника у каждого скана.

#include "delay_api.h"
#include "functions.h"
#include "catalog_bridge.h"
#include "READ_CAT.h"
#include "constants.h"
#include <cstring>
#include <set>

namespace ariadna {

bool compute_task_polys(const std::string& cfx_path, const std::string& orbit_path,
                        const std::string& eop_path, double block_sec, int degree,
                        double sample_sec, bool with_tropo, TaskPolys& out) {
    out.delay.clear(); out.uvw.clear();

    CfxTask task;
    if (!parse_cfx(cfx_path, task)) { std::fprintf(stderr, "compute_task_polys: не разобрать %s\n", cfx_path.c_str()); return false; }
    if (task.scans.empty()) { std::fprintf(stderr, "compute_task_polys: нет сканов\n"); return false; }

    // Орбита космической станции: из параметра либо из cfx (ORB_FILE).
    std::vector<SpaceStation> orbit;
    bool has_space = false; std::string cfx_orb;
    for (const auto& s : task.stations) if (s.is_space) { has_space = true; if (!s.orb_file.empty()) cfx_orb = s.orb_file; }
    std::string orb = orbit_path.empty() ? cfx_orb : orbit_path;
    if (has_space && !orb.empty()) read_scf_orbit(orb, orbit);

    // Узлы EOP (7 шт.) вокруг начала сеанса.
    int mjd0 = task.scans.front().mjd;
    std::vector<eop_record> raw_eop; char epath[256]; std::strncpy(epath, eop_path.c_str(), 255); epath[255] = 0;
    std::vector<EOPData> eop;
    if (ReadEOP(raw_eop, epath) == 0 || !raw_eop.empty()) select_eop_nodes(raw_eop, mjd0, cnst::EOP_NDATA, eop);
    if (eop.size() < (size_t)cnst::EOP_NDATA) { std::fprintf(stderr, "compute_task_polys: мало узлов EOP\n"); return false; }

    // Физика станций из каталогов (та же папка, что EOP): океан/атм нагрузка + термопараметры.
    std::string catdir = eop_path;
    size_t sl = catdir.find_last_of("/\\"); catdir = (sl == std::string::npos) ? "" : catdir.substr(0, sl + 1);
    std::vector<Station> phys(task.stations.size());
    for (size_t i = 0; i < task.stations.size(); ++i) phys[i].name = task.stations[i].name;
    {
        char p_oc[256], p_at[256], p_an[256];
        std::snprintf(p_oc, 256, "%sVLBI_ocload_40.cat", catdir.c_str());
        std::snprintf(p_at, 256, "%sVLBI_atmload4_12.cat", catdir.c_str());
        std::snprintf(p_an, 256, "%santenna-info.cat", catdir.c_str());
        std::vector<oc_record> oc; std::vector<atm_record> at; std::vector<ant_record> an;
        ReadOC(oc, p_oc);        map_ocean_tides_to_stations(oc, phys);
        ReadATM(at, p_at);       map_atm_loading_to_stations(at, phys);
        ReadANT_INFO(an, p_an);  build_def_par_from_ant_info(an, phys);
    }

    auto find_source = [&](const std::string& nm) -> CfxSource {
        for (const auto& s : task.sources) if (s.name == nm) return s;
        CfxSource z; z.name = nm; return z;
    };
    auto participates = [&](const CfxStation& st, const CfxScan& sc) -> bool {
        if (sc.tel_iam.empty()) return true;
        for (const auto& t : sc.tel_iam) if (t == st.iam) return true;
        return false;
    };

    // Посканный расчёт: у каждого скана свой источник; станция считается только в сканах, где
    // участвует. Полиномы задержки и координат (u,v,w) копятся отдельно по каждой станции.
    for (size_t i = 0; i < task.stations.size(); ++i) {
        const CfxStation& st = task.stations[i];
        StationPoly poly; poly.telescope = st.name; poly.order = degree + 1;
        StationUvw uvw; uvw.telescope = st.name; uvw.order = degree + 1;
        for (const auto& sc : task.scans) {
            if (!participates(st, sc)) continue;
            CfxSource src = find_source(sc.source);
            StationPoly sp = compute_station_poly(st, src, sc.mjd, sc.utc, sc.dur_sec,
                                                  eop, block_sec, degree, sample_sec, orbit, with_tropo, phys[i], &uvw);
            for (const auto& b : sp.blocks) poly.blocks.push_back(b);
            if (poly.source.empty()) poly.source = src.name;
        }
        out.delay.push_back(std::move(poly));
        out.uvw.push_back(std::move(uvw));
    }
    return true;
}

} // namespace ariadna
