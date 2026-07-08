// test_process_ariadna.cpp
//
// Проверка верхнего слоя оркестрации process_ariadna: полный пакетный проход по
// наблюдениям (site + source_vec один раз, затем цикл: время -> EOP(interp_iers)
// -> эфемериды -> r2000 -> compute_delay_obs -> файл). Прогон на эталонном
// наблюдении сессии 18JAN02XA (FORTLEZA-HART15M, 0017+200) с РЕАЛЬНЫМИ узлами
// каталога EOPC04 (MJD 58117..58123). Сверка итоговой задержки из выходного файла
// с дампом ARIADNA.
//
// Отличие от test_compute_delay: EOP интерполируются внутри process_ariadna из
// каталожных узлов (INTERP_EOP40, как ARIADNA), а океаническая/атмосферная нагрузки
// и термопараметры антенн читаются из РЕАЛЬНЫХ каталогов (VLBI_ocload/atmload,
// antenna-info) — проверяется вся цепочка оркестрации и адаптеры каталогов.
//
// ОСТАТОК ~5 см обусловлен: (1) суб-суточной поправкой UT1 (~5e-5 с -> 24 мм), которую
// оба независимых интерполятора (interp_eop terms_71/lib и interp_iers PMUT1/GRAVI)
// дают согласованно, но которая уходит от «сырой» интерполяции дампа ARIADNA;
// (2) аргументом времени эфемерид (TDB). Оба требуют дампов ARIADNA для бит-сверки.

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

int main() {
    printf("=====================================================================\n");
    printf("  UNIT: process_ariadna -- пакетный проход vs дамп (18JAN02XA)\n");
    printf("---------------------------------------------------------------------\n");

    try { init_ephemeris(find_eph()); }
    catch (const std::exception& e) { printf("SKIP: эфемериды недоступны: %s\n", e.what()); return 2; }

    // --- Станции (координаты каталога @2000.0 + скорости) ---
    std::vector<Station> st(2);
    st[0].name = "FORTLEZA"; st[0].axsty = "AZEL"; st[0].offs = 0.00637;
    st[0].xyz  = Eigen::Vector3d(4985370.025, -3955020.358, -428472.184);
    st[0].vel  = Eigen::Vector3d(-0.0024, -0.0039, 0.0126);
    st[1].name = "HART15M";  st[1].axsty = "AZEL"; st[1].offs = 1.49500;
    st[1].xyz  = Eigen::Vector3d(5085490.783, 2668161.315, -2768692.721);
    st[1].vel  = Eigen::Vector3d(0.0019, 0.0216, 0.0133);

    // --- Нагрузки и термопараметры — из РЕАЛЬНЫХ каталогов ---
    {
        char p_oc[256]  = "external/catalogs/VLBI_ocload_40.cat";
        char p_atm[256] = "external/catalogs/VLBI_atmload4_12.cat";
        char p_ant[256] = "external/catalogs/antenna-info.cat";
        std::vector<oc_record> oc; std::vector<atm_record> atm; std::vector<ant_record> ant;
        ReadOC(oc, p_oc);         map_ocean_tides_to_stations(oc, st);
        ReadATM(atm, p_atm);      map_atm_loading_to_stations(atm, st);
        ReadANT_INFO(ant, p_ant); build_def_par_from_ant_info(ant, st);
        printf("  каталоги: oc=%zu atm=%zu ant=%zu\n", oc.size(), atm.size(), ant.size());
        printf("  FORTLEZA: hf=%.3f gamma_hf=%.2e hp=%.3f t_0=%.1f p_0=%.1f oc_amp[V,M2]=%.5f atm_has=%d\n",
               st[0].def_par.hf, st[0].def_par.gamma_hf, st[0].def_par.hp, st[0].def_par.t_0,
               st[0].atm_load.p_0, st[0].tide_data.amplitudes(0, 0), (int)st[0].atm_load.has_data);
        printf("  HART15M : hf=%.3f gamma_hf=%.2e hp=%.3f t_0=%.1f p_0=%.1f oc_amp[V,M2]=%.5f atm_has=%d\n",
               st[1].def_par.hf, st[1].def_par.gamma_hf, st[1].def_par.hp, st[1].def_par.t_0,
               st[1].atm_load.p_0, st[1].tide_data.amplitudes(0, 0), (int)st[1].atm_load.has_data);
    }

    // --- Источник 0017+200 ---
    std::vector<Source> src(1);
    src[0].name = "0017+200"; src[0].ra = 0.085655996508366; src[0].dec = 0.355395794394430;
    src[0].ra_rate = 0.0; src[0].dec_rate = 0.0;

    // --- Наблюдение ---
    // ВНИМАНИЕ: UTC берём с полной точностью. Дамп даёт TT-долю 0.709643333
    // (jd_tt=2458121.209643333); UTC = TT - 69.184с/86400 = 0.708842592.
    // Грубое 0.7088 теряет ~3.7 с -> поворот Земли -> ошибка задержки ~км.
    const double utc_obs = 0.708842592;
    std::vector<Observation> obs(1);
    obs[0].mjd = 58120; obs[0].utc = utc_obs;
    obs[0].sta1 = 0; obs[0].sta2 = 1; obs[0].sou = 0;
    obs[0].t1 = 30.427; obs[0].p1 = 1006.7; obs[0].e1 = 60.973;
    obs[0].t2 = 25.884; obs[0].p2 = 861.233; obs[0].e2 = 43.395;

    // --- EOP: реальные узлы EOPC04 (MJD, x", y", UT1-UTC s) вокруг 58120 ---
    const double eop_rows[7][4] = {
        {58117, 0.063071, 0.245448, 0.2182464},
        {58118, 0.061214, 0.246580, 0.2172353},
        {58119, 0.059224, 0.247646, 0.2163654},
        {58120, 0.057406, 0.248566, 0.2156078},
        {58121, 0.055667, 0.249547, 0.2148691},
        {58122, 0.054388, 0.250409, 0.2140616},
        {58123, 0.053497, 0.251568, 0.2131861},
    };
    std::vector<EOPData> eop(7);
    for (int i = 0; i < 7; ++i) {
        eop[i].mjd = eop_rows[i][0];
        eop[i].x = eop_rows[i][1]; eop[i].y = eop_rows[i][2];
        eop[i].ut1_utc = eop_rows[i][3];
        double idelt; nsec(eop_rows[i][0], idelt);
        eop[i].ut1_tai = eop_rows[i][3] - idelt; // UT1-TAI = (UT1-UTC) - (TAI-UTC); нужно interp_eop
        eop[i].dpsi = 0.0; eop[i].deps = 0.0;
    }

    // --- Прогон: результат В ПАМЯТИ (без промежуточного файла) ---
    std::vector<SpaceStation> space; std::vector<OrbitData> orb;
    std::vector<DelayResult> res;
    double mjd_mean = 58120, utc_mean = utc_obs, t_mean = 58120.0 + utc_obs;
    process_ariadna(st, src, obs, space, orb, 1, 0, 0, 0.0, "", eop,
                    mjd_mean, utc_mean, t_mean, res);
    if (res.empty()) { printf("FAIL: пустой результат\n"); return 1; }
    double tau = res[0].tau;

    const double ref = -0.3450839807711632e-03;
    const double d = std::fabs(tau - ref);
    printf("  process_ariadna: tau  = % .12e\n", tau);
    printf("  дамп ARIADNA:    ref  = % .12e\n", ref);
    printf("  |tau - ref| = %.3e с  (~%.2f см)\n", d, d * 3e8 * 100);
    bool ok = d < 5e-9;
    printf("---------------------------------------------------------------------\n");
    printf("  Оркестрация (site+source+EOP interp+цикл) даёт задержку близко к ARIADNA (<5e-9 с): %s\n", ok ? "OK" : "FAIL");
    printf("  РЕЗУЛЬТАТ: %s\n", ok ? "PASS" : "FAIL");
    return ok ? 0 : 1;
}
