// test_final_build.cpp
//
// ИТОГОВЫЙ БИЛД: сквозной тест всего конвейера задержки через process_ariadna,
// собранный "по эталону" tests/Вывод.txt (реальный вывод ARIADNA сессии 18JAN02XA).
//
// Из Вывод.txt берутся ТОЧНО: каталог 8 станций (координаты/скорости/монтировка/
// смещение оси), источник 0017+200, метеоусловия эталонного наблюдения FORTLEZA-
// HART15M и тропосферные величины (E = 38.6469°/39.7425°). Океаническая/атмосферная
// нагрузка и термопараметры антенн читаются из реальных каталогов. EOP — реальные
// узлы EOPC04 (INTERP_EOP40). Сверяется итоговая теоретическая задержка с дампом.
//
// Тропосферные промежутки (Datmc, Z_d/Z_w, E) из Вывод.txt покомпонентно проверены
// отдельными тестами (test_theor_delay, test_delay_pipeline, test_trop_delay).

#include "../src/functions.h"
#include "../src/catalog_bridge.h"
#include "../src/READ_CAT.h"
#include <cstdio>
#include <cmath>
#include <fstream>
#include <sstream>
#include <string>

using namespace ariadna;

static std::string find_eph() {
    const char* c[] = {"external/dephem-master/linux_p1550p2650.440t",
                       "../external/dephem-master/linux_p1550p2650.440t"};
    for (auto p : c) { std::ifstream f(p, std::ios::binary); if (f.good()) return p; }
    return c[0];
}

// Станция из Вывод.txt: имя, координаты ITRF@2000.0, скорости, монтировка, смещение оси.
static Station mk(const char* name, double x, double y, double z,
                  double vx, double vy, double vz, double offs) {
    Station s; s.name = name; s.axsty = "AZEL"; s.offs = offs;
    s.xyz << x, y, z; s.vel << vx, vy, vz;
    return s;
}

int main() {
    printf("=====================================================================\n");
    printf("  ИТОГОВЫЙ БИЛД: process_ariadna по эталону Вывод.txt (18JAN02XA)\n");
    printf("---------------------------------------------------------------------\n");

    try { init_ephemeris(find_eph()); }
    catch (const std::exception& e) { printf("SKIP: эфемериды недоступны: %s\n", e.what()); return 2; }

    // --- 8 станций ТОЧНО из Вывод.txt (строки 3..10) ---
    std::vector<Station> st;
    st.push_back(mk("FORTLEZA", 4985370.025, -3955020.358,  -428472.184, -0.0024, -0.0039, 0.0126, 0.00637));
    st.push_back(mk("HART15M",  5085490.783,  2668161.315, -2768692.721,  0.0019,  0.0216, 0.0133, 1.49500));
    st.push_back(mk("NYALES20", 1202462.642,   252734.460,  6237766.127, -0.0144,  0.0075, 0.0108, 0.52420));
    st.push_back(mk("WETTZ13N", 4075627.777,   931774.140,  4801552.281, -0.0159,  0.0170, 0.0102, 0.00000));
    st.push_back(mk("KATH12M", -4147354.357,  4581542.476, -1573303.689, -0.0363, -0.0097, 0.0588, 0.00000));
    st.push_back(mk("YARRA12M",-2388895.733,  5043349.924, -3078591.260, -0.0494,  0.0091, 0.0505, 0.00000));
    st.push_back(mk("KOKEE",   -5543837.699, -2054567.352,  2387852.201, -0.0093,  0.0629, 0.0324, 0.52081));
    st.push_back(mk("ISHIOKA", -3959635.864,  3296825.493,  3747042.533, -0.0193, -0.0017, 0.0053, 0.00000));

    // --- Нагрузки и термопараметры из реальных каталогов ---
    {
        char p_oc[256]  = "external/catalogs/VLBI_ocload_40.cat";
        char p_atm[256] = "external/catalogs/VLBI_atmload4_12.cat";
        char p_ant[256] = "external/catalogs/antenna-info.cat";
        std::vector<oc_record> oc; std::vector<atm_record> atm; std::vector<ant_record> ant;
        ReadOC(oc, p_oc);         map_ocean_tides_to_stations(oc, st);
        ReadATM(atm, p_atm);      map_atm_loading_to_stations(atm, st);
        ReadANT_INFO(ant, p_ant); build_def_par_from_ant_info(ant, st);
        int oc_n = 0, atm_n = 0, ant_n = 0;
        for (auto& s : st) {
            if (s.tide_data.amplitudes.cwiseAbs().sum() > 0) ++oc_n;
            if (s.atm_load.has_data) ++atm_n;
            if (s.def_par.hf != 1.0) ++ant_n; // не дефолт мобильной антенны
        }
        printf("  каталоги подхвачены: океан %d/8, атм %d/8, antenna-info %d/8 станций\n", oc_n, atm_n, ant_n);
    }

    // --- Источник 0017+200 (Вывод.txt строка 11) ---
    std::vector<Source> src(1);
    src[0].name = "0017+200"; src[0].ra = 0.085655996508366; src[0].dec = 0.355395794394430;

    // --- Эталонное наблюдение FORTLEZA-HART15M (метео из Вывод.txt) ---
    // UTC полной точности: TT-доля 0.709643333 (jd_tt=2458121.209643333) - 69.184с.
    const double utc_obs = 0.708842592;
    std::vector<Observation> obs(1);
    obs[0].mjd = 58120; obs[0].utc = utc_obs;
    obs[0].sta1 = 0; obs[0].sta2 = 1; obs[0].sou = 0;
    obs[0].t1 = 30.427; obs[0].p1 = 1006.7; obs[0].e1 = 60.973;   // FORTLEZA
    obs[0].t2 = 25.884; obs[0].p2 = 861.233; obs[0].e2 = 43.395;  // HART15M

    // --- EOP: реальные узлы EOPC04 вокруг 58120 ---
    const double eop_rows[7][4] = {
        {58117, 0.063071, 0.245448, 0.2182464}, {58118, 0.061214, 0.246580, 0.2172353},
        {58119, 0.059224, 0.247646, 0.2163654}, {58120, 0.057406, 0.248566, 0.2156078},
        {58121, 0.055667, 0.249547, 0.2148691}, {58122, 0.054388, 0.250409, 0.2140616},
        {58123, 0.053497, 0.251568, 0.2131861},
    };
    std::vector<EOPData> eop(7);
    for (int i = 0; i < 7; ++i) {
        eop[i].mjd = eop_rows[i][0]; eop[i].x = eop_rows[i][1]; eop[i].y = eop_rows[i][2];
        eop[i].ut1_utc = eop_rows[i][3];
        double idelt; nsec(eop_rows[i][0], idelt); eop[i].ut1_tai = eop_rows[i][3] - idelt;
        eop[i].dpsi = 0.0; eop[i].deps = 0.0;
    }

    // --- Полный прогон: результат В ПАМЯТИ (без промежуточного файла) ---
    std::vector<SpaceStation> space; std::vector<OrbitData> orb;
    std::vector<DelayResult> res;
    double mjd_mean = 58120, utc_mean = utc_obs, t_mean = 58120.0 + utc_obs;
    process_ariadna(st, src, obs, space, orb, 1, 2, 4, 0.0, "", eop, mjd_mean, utc_mean, t_mean, res);
    if (res.empty()) { printf("FAIL: пустой результат\n"); return 1; }
    double tau = res[0].tau;

    const double ref = -0.3450839807711632e-03;
    const double d = std::fabs(tau - ref);
    printf("  process_ariadna: tau = % .12e с\n", tau);
    printf("  дамп ARIADNA:    ref = % .12e с\n", ref);
    printf("  |tau - ref| = %.3e с  (~%.2f см)\n", d, d * 3e8 * 100);
    bool ok = d < 5e-9;
    printf("---------------------------------------------------------------------\n");
    printf("  Итоговый билд воспроизводит задержку ARIADNA (<5e-9 с): %s\n", ok ? "OK" : "FAIL");
    printf("  РЕЗУЛЬТАТ: %s\n", ok ? "PASS" : "FAIL");
    return ok ? 0 : 1;
}
