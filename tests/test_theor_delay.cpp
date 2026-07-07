// test_theor_delay.cpp
//
// Сверка theor_delay с ПОЛНЫМ отладочным дампом Fortran-подпрограммы TEOR_DELAY
// (операционный прогон ARIADNA). Проверяются как промежуточные величины
// (независимый пересчёт из входных данных), так и финальный выход самой theor_delay.
//
// Раскладка матриц 2x2 (Datmc_d/Datmc_w/dtau_off/dt_temp): строка = задержка(0)/скорость(1),
// столбец = станция 1(0)/2(1).
//
// Сборка (из корня репозитория):
//   g++ -std=c++17 -I.\external\eigen -I.\external ^
//       src\theor_delay.cpp tests\test_theor_delay.cpp -o tests\test_theor_delay.exe

#include "../src/functions.h"
#include <cstdio>
#include <cmath>

using namespace ariadna;

static int g_fail = 0;
static void chk(const char* name, double got, double exp, double tol) {
    double d = std::fabs(got - exp);
    bool ok = d <= tol;
    if (!ok) ++g_fail;
    printf("  %-22s calc=% .16e  dump=% .16e  |d|=%.2e  %s\n",
           name, got, exp, d, ok ? "OK" : "FAIL");
}

int main() {
    printf("=====================================================================\n");
    printf("  UNIT: theor_delay vs Fortran TEOR_DELAY dump\n");
    printf("---------------------------------------------------------------------\n");

    const double C = cnst::C, C2 = C * C;

    Eigen::Vector3d K_s(0.9340717157118170, 0.08020509321046379, 0.3479614532247550);

    Eigen::Matrix<double, 3, 2> base_line;
    base_line.col(0) << 415219.8984698132, 6610745.296371628, -2340691.244571441;
    base_line.col(1) << -482.0557489968165, 30.57309281858454, 0.8338814000745651;

    std::vector<Eigen::Vector3d> xsta(2), vsta(2), asta(2);
    xsta[0] << 4788183.14731, -4190717.19517, -436901.22171;
    xsta[1] << 5203403.04578, 2420028.10120, -2777592.46628;

    Eigen::Matrix3d Earth, Sun, Moon;
    Earth.col(0) << -30331171034.03649, 132858910822.4568, 57575534750.98912;
    Earth.col(1) << -29619.14420142996, -5774.970368033236, -2504.268390466165;
    Earth.col(2) << -0.02643228993058815, -0.01833876281565307, -0.002325872820870766;

    Sun.col(0) << 268184443.84710, 849738216.96790, 348845058.27214;
    Sun.col(1) << -10.14974, 7.86276, 3.67712;
    Sun.col(2).setZero();

    Moon.col(0) << -30457209591.74120, 133170659155.41185, 57696076560.55701;
    Moon.col(1) << -30654.27062, -6147.62886, -2567.30178;
    Moon.col(2).setZero();

    vsta[0] = Eigen::Vector3d(-29313.55154, -5425.75602, -2504.78267) - Earth.col(1);
    vsta[1] = Eigen::Vector3d(-29795.60728676801, -5395.182928854237, -2503.948787164999) - Earth.col(1);
    asta[1] << 0, 0, 0;
    asta[0] = asta[1] - Eigen::Vector3d(-0.002229426718217218, -0.03515211179078134, 0.000002485990340675793);

    // Раскладка: строка0 = задержки станций 1,2; строка1 = скорости станций 1,2.
    Eigen::Matrix2d Datmc_d, Datmc_w, dtau_off, dt_temp;
    Datmc_d << -0.1225039869341168e-07, 0.1023313654522243e-07, 0.9572258037255259e-12, 0.3346159616445429e-12;
    Datmc_w << -0.1348080932573043e-08, 0.7304528775933354e-09, 0.1057065200668694e-12, 0.2396256388890582e-13;
    dtau_off << 0.1658888657089142e-10, -0.3833381489451837e-08, -0.8349692241350816e-15, -0.8724128831589889e-13;
    dt_temp << 0.8286346437956392e-10, 0.1445568241005169e-09, 0, 0;

    // -------- Независимый пересчёт промежуточных величин, сверка с дампом --------
    printf("  [Промежуточные величины: независимый пересчёт vs дамп]\n");
    Eigen::Vector3d b = base_line.col(0), db = base_line.col(1);
    double K_starB = K_s.dot(b);
    double dotK_starB = K_s.dot(db);
    chk("K_starB",   K_starB,   0.1035902786359673e+06, 1e-4);
    chk("dotK_starB",dotK_starB,-0.4475323641911676e+03, 1e-10);

    Eigen::Vector3d R_geocen = Earth.col(0) - Sun.col(0);
    double R_geo = R_geocen.norm();
    double U = (cnst::GSUN / C2) / R_geo;
    chk("R_geo", R_geo, 0.1470973698637544e+12, 10.0);
    chk("U",     U,     0.1003841904089365e-07, 1e-18);

    double V_Earth = Earth.col(1).norm();
    chk("V_Earth", V_Earth, 0.3028061006895753e+05, 1e-8);

    double term2a = K_starB / C;
    double term2b = 1.0 - 2.0 * U;
    double term2c = Earth.col(1).dot(Earth.col(1)) / (2.0 * C2);
    double term2d = Earth.col(1).dot(vsta[1]) / C2;
    chk("term2a", term2a, 0.3455399756453089e-03, 1e-15);
    chk("term2b", term2b, 0.9999999799231619e+00, 1e-14);
    chk("term2c", term2c, 0.5101029556441379e-08, 1e-18);
    chk("term2d", term2d, 0.3374249253578403e-10, 1e-20);

    // term2 в дампе — БЕЗ term2e/f (операционный прогон), поэтому сверяем без них.
    double term2_noef = term2a * (term2b - term2c - term2d);
    chk("term2 (no e/f)", term2_noef, 0.3455399669336898e-03, 1e-15);

    double term3a = Earth.col(1).dot(b) / C2;
    double term3b = 1.0 + K_s.dot(Earth.col(1)) / (2.0 * C);
    double term3 = term3a * term3b;
    chk("term3a", term3a, -0.4963932136927889e-06, 1e-18);
    chk("term3b", term3b,  0.9999516315788251e+00, 1e-14);
    chk("term3",  term3,  -0.4963692039367605e-06, 1e-18);

    double den_Eq9 = 1.0 + K_s.dot(Earth.col(1) + vsta[1]) / C;
    chk("den_Eq9", den_Eq9, 0.9999028153242351e+00, 1e-14);

    // tv2_tv1 через term1 из дампа (Шапиро проверяется транзитивно финалом).
    double term1_dump = -0.3334437597304269e-09;
    double numer_Eq9 = term1_dump - term2_noef - term3;
    double tv2_tv1 = numer_Eq9 / den_Eq9;
    chk("numer_Eq9", numer_Eq9, -0.3450439311735128e-03, 1e-15);
    chk("tv2_tv1",   tv2_tv1,   -0.3450774674152973e-03, 1e-15);

    // w2_w1 = vsta2 - vsta1; дамп подтверждает w2_w1 == base_line.col(1) (и K_dotw2w1 == dotK_starB).
    // Сверяем через base_line: dotX_1 в дампе напечатан лишь с ~8 знаками, поэтому
    // vsta[1]-vsta[0] в тесте теряет точность (эффект на финал ~1e-22, ниже шума 2.8e-14).
    double K_dotw2w1 = K_s.dot(base_line.col(1));
    chk("K_dotw2w1", K_dotw2w1, -0.4475323641911676e+03, 1e-10);

    double Datmc11 = Datmc_d(0, 0) + Datmc_w(0, 0);
    double Datmc21 = Datmc_d(0, 1) + Datmc_w(0, 1);
    chk("Datmc(1) delay", Datmc11, -0.1359847962598472e-07, 1e-20);
    chk("Datmc(2) delay", Datmc21,  0.1096358942281577e-07, 1e-20);

    double tg2_tg1 = tv2_tv1 + Datmc11 * K_dotw2w1 / C;
    double t2_t1_a = tg2_tg1 + Datmc11 + Datmc21;
    chk("tg2_tg1",  tg2_tg1,  -0.3450774673949974e-03, 1e-15);
    chk("t2_t1_a",  t2_t1_a,  -0.3450801022852005e-03, 1e-15);

    // -------- Финальный выход самой theor_delay (с term2e/f, ~1e-15 ниже шума) --------
    printf("\n  [Финал theor_delay vs дамп]\n");
    double t2_t1, dt2_t1;
    theor_delay(base_line, xsta, vsta, asta, K_s, Earth, Sun, Moon,
                Datmc_d, Datmc_w, dtau_off, dt_temp, t2_t1, dt2_t1);
    chk("t2_t1",  t2_t1,  -0.3450839807711632e-03, 1e-13);
    chk("dt2_t1", dt2_t1,  0.1492797110326536e-05, 5e-12); // ~1e-12 предсущ. шум (по согласованию оставлен)

    printf("=====================================================================\n");
    printf("  РЕЗУЛЬТАТ: %s (провалов: %d)\n", g_fail == 0 ? "PASS" : "FAIL", g_fail);
    return g_fail == 0 ? 0 : 1;
}
