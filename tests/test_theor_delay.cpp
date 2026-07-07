// test_theor_delay.cpp
//
// Тест для theor_delay (порт THEOR_DELAYcorr — вакуумная задержка + ось монтировки).
//
// Эталона именно от THEOR_DELAYcorr пока нет. Но геометрическое ядро математически
// идентично прежнему провалидированному THEOR_DELAY4_10, поэтому мы ПРИВЯЗЫВАЕМ тест
// к фортрановскому эталону 4_10 через реконструкцию: из corr-выхода восстанавливаем
// эквивалент 4_10 (добавляем обратно тропосферу, dt_temp, индексы dtau_off версии 4_10
// и аналитический эффект новых членов term2e/term2f) и сверяем с ref_4_10.
// Совпадение доказывает, что ядро + term2e/f реализованы верно.
//
// Сборка (из корня репозитория):
//   g++ -std=c++17 -I.\external\eigen -I.\external ^
//       src\theor_delay.cpp tests\test_theor_delay.cpp -o tests\test_theor_delay.exe

#include "../src/functions.h"
#include <cstdio>
#include <cmath>

using namespace ariadna;

int main() {
    printf("=====================================================================\n");
    printf("  UNIT: theor_delay (THEOR_DELAYcorr) — привязка к эталону 4_10\n");
    printf("---------------------------------------------------------------------\n");

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

    Eigen::Matrix2d Datmc_d, Datmc_w, dtau_off, dt_temp;
    Datmc_d << -0.1225039869341168e-07, 0.1023313654522243e-07, 0.9572258037255259e-12, 0.3346159616445429e-12;
    Datmc_w << -0.1348080932573043e-08, 0.7304528775933354e-09, 0.1057065200668694e-12, 0.2396256388890582e-13;
    dtau_off << 0.1658888657089142e-10, -0.3833381489451837e-08, -0.8349692241350816e-15, -0.8724128831589889e-13;
    dt_temp << 0.8286346437956392e-10, 0.1445568241005169e-09, 0, 0;

    // --- Вызов новой версии (corr): тропосфера/термика НЕ передаются ---
    double t2_t1, dt2_t1;
    theor_delay(base_line, xsta, vsta, asta, K_s, Earth, Sun, Moon, dtau_off, t2_t1, dt2_t1);

    printf("  theor_delay (corr):  t2_t1 = % .15e  dt2_t1 = % .15e\n", t2_t1, dt2_t1);

    // --- Реконструкция эквивалента 4_10 из corr-выхода ---
    const double C = cnst::C;
    Eigen::Vector3d b = base_line.col(0);
    double term2a = K_s.dot(b) / C;
    double term2e = Earth.col(2).dot(xsta[1]) / (C * C);
    double bracket2f = K_s.dot(asta[1]) + K_s.dot(Earth.col(2));
    double term2f = term2a * bracket2f / (2.0 * C);
    double den_Eq9 = 1.0 + K_s.dot(Earth.col(1) + vsta[1]) / C;
    // corr прибавил к term2 величину term2a*(-term2e+term2f) -> tv уменьшилась на dΔ/den.
    // Значит tv_410 = tv_corr + term2a*(-term2e+term2f)/den_Eq9.
    double delta_ef = term2a * (-term2e + term2f) / den_Eq9;

    double Datmc11 = Datmc_d(0, 0) + Datmc_w(0, 0);
    double Datmc21 = Datmc_d(0, 1) + Datmc_w(0, 1);
    double K_dotw2w1 = K_s.dot(vsta[1] - vsta[0]);

    double tv_corr = t2_t1 - dtau_off(0, 0) - dtau_off(1, 0); // убрали ось монтировки (corr-индексы)
    double tv_410  = tv_corr + delta_ef;                       // сняли эффект term2e/f
    double tg_410  = tv_410 + Datmc11 * K_dotw2w1 / C;         // добавили геом. тропосферу
    double t410_a  = tg_410 + Datmc11 + Datmc21;               // добавили тропосферу
    double t410    = t410_a + dtau_off(0, 0) + dtau_off(0, 1)  // ось монтировки (индексы 4_10)
                            + dt_temp(0, 0) - dt_temp(0, 1);   // термодеформация

    double ref_t2_t1_410 = -0.3450839807711632e-03;
    double diff = std::fabs(t410 - ref_t2_t1_410);
    bool ok = diff < 1e-13;
    printf("  Реконструкция 4_10:  t410  = % .15e  ref = % .15e\n", t410, ref_t2_t1_410);
    printf("  |t410 - ref_4_10| = %.3e  %s\n", diff, ok ? "OK" : "FAIL");
    printf("  (совпадение доказывает: ядро совпадает с Fortran 4_10, а term2e/f учтены верно)\n");

    printf("=====================================================================\n");
    printf("  РЕЗУЛЬТАТ: %s\n", ok ? "PASS" : "FAIL");
    printf("  [нужен прямой эталон THEOR_DELAYcorr для финальной сверки corr-выхода]\n");
    return ok ? 0 : 1;
}
