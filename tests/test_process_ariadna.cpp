// test_process_ariadna.cpp
//
// Проверка верхнего слоя оркестрации process_ariadna: полный пакетный проход по
// наблюдениям (site + source_vec один раз, затем цикл: время -> EOP(interp_iers)
// -> эфемериды -> r2000 -> compute_delay_obs -> файл). Прогон на эталонном
// наблюдении сессии 18JAN02XA (FORTLEZA-HART15M, 0017+200) с РЕАЛЬНЫМИ узлами
// каталога EOPC04 (MJD 58117..58123). Сверка итоговой задержки из выходного файла
// с дампом ARIADNA.
//
// Отличие от test_compute_delay: здесь EOP не задаются готовыми числами, а
// интерполируются внутри process_ariadna из каталожных узлов — проверяется вся
// цепочка оркестрации, включая выбор окна EOP и interp_iers.

#include "../src/functions.h"
#include "../src/catalog_bridge.h"
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

    // --- Станции (координаты каталога @2000.0 + скорости; нагрузки обнулены) ---
    std::vector<Station> st(2);
    st[0].name = "FORTLEZA"; st[0].axsty = "AZEL"; st[0].offs = 0.00637;
    st[0].xyz  = Eigen::Vector3d(4985370.025, -3955020.358, -428472.184);
    st[0].vel  = Eigen::Vector3d(-0.0024, -0.0039, 0.0126);
    // ANTENNA_INFO (Nothnagel): t_0, hf/gamma_hf, hp/gamma_hp для термодеформации
    st[0].def_par.t_0 = 26.7; st[0].def_par.hf = 4.44; st[0].def_par.gamma_hf = 1.0e-5;
    st[0].def_par.hp = 3.71;  st[0].def_par.gamma_hp = 1.2e-5;
    st[1].name = "HART15M";  st[1].axsty = "AZEL"; st[1].offs = 1.49500;
    st[1].xyz  = Eigen::Vector3d(5085490.783, 2668161.315, -2768692.721);
    st[1].vel  = Eigen::Vector3d(0.0019, 0.0216, 0.0133);
    st[1].def_par.t_0 = 16.1; st[1].def_par.hf = 6.32; st[1].def_par.gamma_hf = 0.8e-5;
    st[1].def_par.hp = 3.36;  st[1].def_par.gamma_hp = 1.2e-5;
    for (auto& s : st) {
        s.tide_data.amplitudes.setZero(); s.tide_data.phases.setZero();
        s.atm_load.coef.setZero(); s.atm_load.has_data = false;
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
        eop[i].ut1_utc = eop_rows[i][3]; eop[i].ut1_tai = 0.0;
        eop[i].dpsi = 0.0; eop[i].deps = 0.0;
    }

    // --- Прогон ---
    std::string out_path = "tests/out_process_ariadna.txt";
    std::vector<SpaceStation> space; std::vector<OrbitData> orb;
    double mjd_mean = 58120, utc_mean = utc_obs, t_mean = 58120.0 + utc_obs;
    process_ariadna(st, src, obs, space, orb, 1, 0, 0, 0.0, out_path, eop,
                    mjd_mean, utc_mean, t_mean);

    // --- Чтение результата ---
    std::ifstream fin(out_path);
    if (!fin) { printf("FAIL: нет выходного файла %s\n", out_path.c_str()); return 1; }
    double tau = 0.0; bool got = false;
    std::string line;
    while (std::getline(fin, line)) {
        if (line.empty() || line[0] == '#') continue;
        std::istringstream ss(line);
        int mjd, s1, s2, so; double utc, t, dt;
        if (ss >> mjd >> utc >> s1 >> s2 >> so >> t >> dt) { tau = t; got = true; }
    }
    if (!got) { printf("FAIL: не удалось прочитать задержку из файла\n"); return 1; }

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
