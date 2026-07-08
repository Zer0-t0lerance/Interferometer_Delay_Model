// verify.cpp
//
// ПОДРОБНАЯ СВЕРКА полного конвейера с эталонным дампом ARIADNA (tests/Вывод.txt),
// сессия 18JAN02XA, база FORTLEZA-HART15M, источник 0017+200. Показывает по всему
// циклу: РЕФЕРЕНС (из дампа) | ПОСЧИТАНО | |ОШИБКА| — в консоль И в файл-отчёт.
//
// Данные идут В ПАМЯТИ: process_ariadna возвращает results + debug (промежуточные
// величины), БЕЗ промежуточных файлов (files только для человекочитаемого отчёта).
//
// Референсные числа взяты из tests/Вывод.txt (Debug output TROP_DELAY / itог задержки).

#include "../src/functions.h"
#include "../src/catalog_bridge.h"
#include "../src/READ_CAT.h"
#include <cstdio>
#include <cmath>
#include <fstream>
#include <string>
#include <vector>

using namespace ariadna;

static std::string find_eph() {
    const char* c[] = {"external/dephem-master/linux_p1550p2650.440t",
                       "../external/dephem-master/linux_p1550p2650.440t"};
    for (auto p : c) { std::ifstream f(p, std::ios::binary); if (f.good()) return p; }
    return c[0];
}

// Одна строка отчёта: имя | референс | посчитано | |ошибка| (в консоль + файл).
static int g_fail = 0;
static void row(std::FILE* fp, const char* name, const char* unit,
                double ref, double calc, double tol) {
    double err = std::fabs(calc - ref);
    const char* verdict = (err <= tol) ? "OK" : "БОЛЬШЕ ДОП.";
    if (err > tol) ++g_fail;
    char line[256];
    std::snprintf(line, sizeof(line), "  %-30s %-6s | % .12e | % .12e | %.3e  [%s]\n",
                  name, unit, ref, calc, err, verdict);
    std::fputs(line, stdout);
    if (fp) std::fputs(line, fp);
}
static void hdr(std::FILE* fp, const char* s) {
    std::fputs(s, stdout); std::fputs("\n", stdout);
    if (fp) { std::fputs(s, fp); std::fputs("\n", fp); }
}

int main() {
    std::FILE* fp = std::fopen("tests/verify_report.txt", "w");

    hdr(fp, "=====================================================================================");
    hdr(fp, "  ПОДРОБНАЯ СВЕРКА КОНВЕЙЕРА С ДАМПОМ ARIADNA (Вывод.txt)");
    hdr(fp, "  Сессия 18JAN02XA | база FORTLEZA-HART15M | источник 0017+200 | MJD 58120");
    hdr(fp, "  Колонки: величина | референс (ARIADNA) | посчитано (наш конвейер) | |ошибка|");
    hdr(fp, "=====================================================================================");

    try { init_ephemeris(find_eph()); }
    catch (const std::exception& e) { std::printf("SKIP: эфемериды недоступны: %s\n", e.what()); if (fp) std::fclose(fp); return 2; }

    // --- Станции (координаты из Вывод.txt) + нагрузки/термо из каталогов ---
    std::vector<Station> st(2);
    st[0].name = "FORTLEZA"; st[0].axsty = "AZEL"; st[0].offs = 0.00637;
    st[0].xyz << 4985370.025, -3955020.358, -428472.184; st[0].vel << -0.0024, -0.0039, 0.0126;
    st[1].name = "HART15M";  st[1].axsty = "AZEL"; st[1].offs = 1.49500;
    st[1].xyz << 5085490.783, 2668161.315, -2768692.721; st[1].vel << 0.0019, 0.0216, 0.0133;
    {
        char p_oc[256] = "external/catalogs/VLBI_ocload_40.cat";
        char p_atm[256] = "external/catalogs/VLBI_atmload4_12.cat";
        char p_ant[256] = "external/catalogs/antenna-info.cat";
        std::vector<oc_record> oc; std::vector<atm_record> atm; std::vector<ant_record> ant;
        ReadOC(oc, p_oc); map_ocean_tides_to_stations(oc, st);
        ReadATM(atm, p_atm); map_atm_loading_to_stations(atm, st);
        ReadANT_INFO(ant, p_ant); build_def_par_from_ant_info(ant, st);
    }

    std::vector<Source> src(1);
    src[0].name = "0017+200"; src[0].ra = 0.085655996508366; src[0].dec = 0.355395794394430;

    const double utc_obs = 0.708842592;
    std::vector<Observation> obs(1);
    obs[0].mjd = 58120; obs[0].utc = utc_obs; obs[0].sta1 = 0; obs[0].sta2 = 1; obs[0].sou = 0;
    obs[0].t1 = 30.427; obs[0].p1 = 1006.7; obs[0].e1 = 60.973;
    obs[0].t2 = 25.884; obs[0].p2 = 861.233; obs[0].e2 = 43.395;

    const double eop_rows[7][4] = {
        {58117, 0.063071, 0.245448, 0.2182464}, {58118, 0.061214, 0.246580, 0.2172353},
        {58119, 0.059224, 0.247646, 0.2163654}, {58120, 0.057406, 0.248566, 0.2156078},
        {58121, 0.055667, 0.249547, 0.2148691}, {58122, 0.054388, 0.250409, 0.2140616},
        {58123, 0.053497, 0.251568, 0.2131861}};
    std::vector<EOPData> eop(7);
    for (int i = 0; i < 7; ++i) {
        eop[i].mjd = eop_rows[i][0]; eop[i].x = eop_rows[i][1]; eop[i].y = eop_rows[i][2];
        eop[i].ut1_utc = eop_rows[i][3];
        double idelt; nsec(eop_rows[i][0], idelt); eop[i].ut1_tai = eop_rows[i][3] - idelt;
    }

    // --- Прогон: результаты и промежутки В ПАМЯТИ (output_path="" -> файла нет) ---
    std::vector<SpaceStation> space; std::vector<OrbitData> orb;
    std::vector<DelayResult> res; std::vector<CompDebug> dbg;
    process_ariadna(st, src, obs, space, orb, 1, 2, 4, 0.0, "", eop,
                    58120, utc_obs, 58120.0 + utc_obs, res, &dbg);
    if (res.empty() || dbg.empty()) { std::printf("FAIL: пустой результат\n"); if (fp) std::fclose(fp); return 1; }
    const CompDebug& d = dbg[0];

    // Наша раскладка: E(станция, значение); Datmc_d(станция, задержка). Станция 0=FORTLEZA, 1=HART15M.
    hdr(fp, "\n-- Геометрия/тропосфера (наблюдение FORTLEZA-HART15M) ------------------------------");
    row(fp, "Элевация FORTLEZA",  "рад", 0.6745164463503079, d.E(0, 0), 5e-6);
    row(fp, "Элевация HART15M",   "рад", 0.6936373986295212, d.E(1, 0), 5e-6);
    // Зенитные задержки — в СЕКУНДАХ (Вывод.txt Zen_dry/Zen_wet); Z_tot в метрах = Zen*c.
    row(fp, "Зенит.сухая FORTLEZA", "с", 0.7665726395296646e-08, d.Zd(0), 5e-13);
    row(fp, "Зенит.сухая HART15M",  "с", 0.6554083158236512e-08, d.Zd(1), 5e-13);
    row(fp, "Зенит.влажн FORTLEZA", "с", 0.8426630615738624e-09, d.Zw(0), 5e-13);
    row(fp, "Зенит.влажн HART15M",  "с", 0.4673903691945521e-09, d.Zw(1), 5e-13);
    row(fp, "Тропо сухая FORTLEZA", "с", -0.1225039869341167e-07, d.Datmc_d(0, 0), 5e-11);
    row(fp, "Тропо сухая HART15M",  "с",  0.1023313654522243e-07, d.Datmc_d(1, 0), 5e-11);
    row(fp, "Тропо влажн FORTLEZA", "с", -0.1348080932573042e-08, d.Datmc_w(0, 0), 5e-11);
    row(fp, "Тропо влажн HART15M",  "с",  0.7304528775933357e-09, d.Datmc_w(1, 0), 5e-11);

    hdr(fp, "\n-- ИТОГ: теоретическая задержка ---------------------------------------------------");
    row(fp, "Задержка t2-t1", "с", -0.3450839807711632e-03, res[0].tau, 5e-9);

    double dtau_err_cm = std::fabs(res[0].tau - (-0.3450839807711632e-03)) * 3e8 * 100;
    char sline[160];
    std::snprintf(sline, sizeof(sline), "\n  Итоговая ошибка задержки: ~%.2f см\n", dtau_err_cm);
    std::fputs(sline, stdout); if (fp) std::fputs(sline, fp);

    hdr(fp, "=====================================================================================");
    const char* verdict = (g_fail == 0) ? "  РЕЗУЛЬТАТ: PASS (все величины в допуске)"
                                         : "  РЕЗУЛЬТАТ: есть расхождения выше допуска (см. таблицу)";
    hdr(fp, verdict);
    hdr(fp, "  Полный отчёт также записан в tests/verify_report.txt");
    hdr(fp, "=====================================================================================");

    if (fp) std::fclose(fp);
    return g_fail == 0 ? 0 : 1;
}
