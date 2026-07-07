// test_delay_pipeline.cpp
//
// Интеграция "низа" конвейера для реального наблюдения (18JAN02XA, FORTLEZA-HART15M,
// 0017+200): baseline -> trop_delay -> theor_delay, сверка ИТОГОВОЙ задержки с дампом.
// Все геометрические/небесные входы взяты из логов старых тестов и дампа theor_delay
// (dtau_off и dt_temp пока логированные: mount_tel/THERM_DEF подключим следующими).
//
// Сборка (из корня репозитория):
//   g++ -std=c++17 -I.\external\eigen -I.\external ^
//       src\baseline.cpp src\trop_delay.cpp src\nhmf2.cpp src\nwmf2.cpp src\sast_dry.cpp ^
//       src\sast_wet.cpp src\theor_delay.cpp tests\test_delay_pipeline.cpp ^
//       -o tests\test_delay_pipeline.exe

#include "../src/functions.h"
#include <cstdio>
#include <cmath>

using namespace ariadna;

static int g_fail = 0;
static void chk(const char* name, double got, double exp, double tol) {
    double d = std::fabs(got - exp);
    bool ok = d <= tol;
    if (!ok) ++g_fail;
    printf("  %-18s calc=% .12e  dump=% .12e  |d|=%.2e  %s\n", name, got, exp, d, ok ? "OK" : "FAIL");
}

int main() {
    printf("=====================================================================\n");
    printf("  UNIT: конвейер baseline->trop_delay->theor_delay vs дамп (18JAN02XA)\n");
    printf("---------------------------------------------------------------------\n");

    // --- Логированные J2000-координаты станций (test_baseline / дамп) ---
    Eigen::MatrixXd xsta(3, 2), vsta(3, 2), asta(3, 2);
    xsta.col(0) << 4788183.147308413, -4190717.195172735, -436901.2217094282;
    xsta.col(1) << 5203403.045778226,  2420028.101198892, -2777592.466280870;
    vsta.col(0) <<  305.5926636587654, 349.2143463604144, -0.5142780989077815;
    vsta.col(1) << -176.4630853380511, 379.7874391789990,  0.3196033011667836;
    asta.setZero();

    // r2000 наблюдения (test_baseline, ITRF->J2000)
    Eigen::Matrix3d r2000;
    r2000 << 0.9988360617015190E+00,  0.4820309780423414E-01,  0.1727196188789892E-02,
            -0.4820310447537553E-01,  0.9988375540050903E+00, -0.3778973045624359E-04,
            -0.1727009998570988E-02, -0.4551047279605037E-04,  0.9999985076815172E+00;

    // --- 1. BASELINE ---
    Eigen::MatrixXd base_line_m(3, 2);
    Eigen::Vector3d b_cfs;
    baseline(r2000, xsta, vsta, asta, base_line_m, b_cfs);
    printf("  [baseline -> вектор базы J2000]\n");
    chk("base_line.x", base_line_m(0,0), 0.4152198984698132E+06, 1e-4);
    chk("base_line.y", base_line_m(1,0), 0.6610745296371628E+07, 1e-4);
    chk("base_line.z", base_line_m(2,0), -0.2340691244571441E+07, 1e-4);

    // --- 2. TROP_DELAY (метео/элевации/азимуты из логов) ---
    Observation obs{};
    obs.t1 = 30.427; obs.p1 = 1006.7; obs.e1 = 60.973;
    obs.t2 = 25.884; obs.p2 = 861.233; obs.e2 = 43.395;
    Station sta1, sta2;
    sta1.lat_geod = -0.06768136568262291; sta1.h_geod = 23.77689785504481;
    sta2.lat_geod = -0.4518611056170612;  sta2.h_geod = 1410.125546193682;

    Eigen::MatrixXd e(2,2), az(2,2);
    e  << 0.6745164463503075, -6.524628647339084e-5,
          0.6936373986295217, -6.196985962630866e-5;
    az << 1.04390201532089,  -1.499847631644231e-5,
          5.85310838945029,  -2.811691949922670e-5;

    Eigen::MatrixXd datmc_d, datmc_w, datmp_hmf, datmp_wmf, dgrad_n, dgrad_e, zen_dry, zen_wet;
    double jd = 2458121.209643333, ct = 0.0;
    trop_delay(obs, jd, ct, sta1, sta2, e, az, datmc_d, datmc_w, datmp_hmf, datmp_wmf,
               dgrad_n, dgrad_e, zen_dry, zen_wet);
    // Допуск 1e-11 как в test_trop_delay: остаточные ~3e-12 — от реверс-инжиниринга
    // производных e/az (собственная точность trop_delay, не ошибка конвейера).
    printf("\n  [trop_delay -> Datmc (сухая/влажная)]\n");
    chk("Datmc_d st1", datmc_d(0,0), -0.1225039869341168e-07, 1e-11);
    chk("Datmc_d st2", datmc_d(1,0),  0.1023313654522243e-07, 1e-11);

    // --- 3. THEOR_DELAY (Солнце/Луна в SSB!) ---
    Eigen::Matrix<double,3,2> base_line;
    base_line.col(0) = base_line_m.col(0);
    base_line.col(1) = base_line_m.col(1);

    std::vector<Eigen::Vector3d> xs = {xsta.col(0), xsta.col(1)};
    std::vector<Eigen::Vector3d> vs = {vsta.col(0), vsta.col(1)};
    std::vector<Eigen::Vector3d> as = {asta.col(0), asta.col(1)};
    Eigen::Vector3d K_s(0.9340717157118170, 0.08020509321046379, 0.3479614532247550);

    Eigen::Matrix3d Earth, Sun, Moon;
    Earth.col(0) << -30331171034.03649, 132858910822.4568, 57575534750.98912;
    Earth.col(1) << -29619.14420142996, -5774.970368033236, -2504.268390466165;
    Earth.col(2) << -0.02643228993058815, -0.01833876281565307, -0.002325872820870766;
    Sun.col(0) << 268184443.84710, 849738216.96790, 348845058.27214;   // SSB
    Sun.col(1) << -10.14974, 7.86276, 3.67712; Sun.col(2).setZero();
    Moon.col(0) << -30457209591.74120, 133170659155.41185, 57696076560.55701; // SSB
    Moon.col(1) << -30654.27062, -6147.62886, -2567.30178; Moon.col(2).setZero();

    // ВАЖНО: trop_delay отдаёт Datmc в раскладке [строка=станция, столбец=задержка/скорость],
    // а theor_delay ждёт [строка=задержка/скорость, столбец=станция] -> транспонируем на стыке.
    Eigen::Matrix2d Datmc_d = Eigen::Matrix2d(datmc_d).transpose();
    Eigen::Matrix2d Datmc_w = Eigen::Matrix2d(datmc_w).transpose();
    Eigen::Matrix2d dtau_off, dt_temp;
    dtau_off << 0.1658888657089142e-10, -0.3833381489451837e-08,
               -0.8349692241350816e-15, -0.8724128831589889e-13;   // из mount_tel (лог)
    dt_temp  << 0.8286346437956392e-10,  0.1445568241005169e-09,
                0.0, 0.0;                                            // из THERM_DEF (лог)

    double t2_t1, dt2_t1;
    theor_delay(base_line, xs, vs, as, K_s, Earth, Sun, Moon, Datmc_d, Datmc_w, dtau_off, dt_temp, t2_t1, dt2_t1);

    // t2_t1 совпадает с дампом в пределах точности Datmc (~1e-12): вся цепочка
    // baseline->trop_delay->theor_delay скомпонована верно.
    printf("\n  [theor_delay -> ИТОГОВАЯ задержка]\n");
    chk("t2_t1",  t2_t1,  -0.3450839807711632e-03, 5e-12);
    chk("dt2_t1", dt2_t1,  0.1492797110326536e-05, 5e-12);

    printf("=====================================================================\n");
    printf("  РЕЗУЛЬТАТ: %s (провалов: %d)\n", g_fail == 0 ? "PASS" : "FAIL", g_fail);
    printf("  [осталось подключить: mount_tel(dtau_off), aber_source(e/az), site(j1,j2), EOP]\n");
    return g_fail == 0 ? 0 : 1;
}
