// test_site_atm40.cpp
//
// Unit-тест для SITE_ATM40 (атмосферная нагрузка).
// Эталона из Fortran нет, поэтому проверки самосогласованные: для контролируемых
// наборов коэффициентов результат сверяется с независимым закрытым расчётом.
//
// Сборка (из корня репозитория):
//   g++ -std=c++17 -Iexternal\eigen -Iexternal ^
//       src\site_atm40.cpp tests\test_site_atm40.cpp -o tests\test_site_atm40.exe

#include "../src/functions.h"
#include <cstdio>
#include <cmath>

using namespace ariadna;

static int g_fail = 0;
static void check(const char* what, double got, double exp, double tol = 1e-12) {
    double d = std::fabs(got - exp);
    bool ok = d <= tol;
    if (!ok) ++g_fail;
    printf("  %-40s got=% .12e exp=% .12e %s\n", what, got, exp, ok ? "OK" : "FAIL");
}

int main() {
    printf("=====================================================================\n");
    printf("  UNIT: SITE_ATM40 (atmospheric loading)\n");
    printf("---------------------------------------------------------------------\n");

    const int mjd = 57422;
    const double utc = 0.5;
    const double Dt = (mjd - cnst::ATM_EPOCH_MJD) + utc; // 11722.5 сут
    const double w1 = cnst::TWOPI / cnst::ATM_PERIODS[0]; // годовая
    const Eigen::Matrix3d I = Eigen::Matrix3d::Identity();
    const Eigen::Matrix3d Z = Eigen::Matrix3d::Zero();

    Eigen::Vector3d dx, dv;

    // --- Случай 1: нулевые коэффициенты -> нулевое смещение ---
    {
        AtmLoadData a; // конструктор обнуляет coef, p_0=0
        SITE_ATM40(mjd, utc, 1010.0, 0.5, a, I, I, Z, dx, dv);
        check("zero-coef dx.norm", dx.norm(), 0.0);
        check("zero-coef dv.norm", dv.norm(), 0.0);
    }

    // --- Случай 2: только постоянный член + регрессия по давлению (вертикаль) ---
    // dr_V = b0 + b1*(P - p_0) = 20 + 3*(1010-1000) = 50 мм = 0.05 м; E,N = 0.
    {
        AtmLoadData a;
        a.coef(0, 6) = 20.0;  // b0 [мм]
        a.coef(0, 7) = 3.0;   // b1 [мм/мбар]
        a.p_0 = 1000.0;
        SITE_ATM40(mjd, utc, 1010.0, 0.0, a, I, I, Z, dx, dv);
        check("const+reg  dx(V)=0.05 m", dx(0), 0.05);
        check("const+reg  dx(E)=0",      dx(1), 0.0);
        check("const+reg  dx(N)=0",      dx(2), 0.0);
        check("const+reg  dv=0 (dPdt=0)", dv.norm(), 0.0);
    }

    // --- Случай 3: годовая косинусная гармоника по вертикали, A1=100 мм ---
    // dr_V = 100*cos(w1*Dt) мм = 0.1*cos(w1*Dt) м.
    {
        AtmLoadData a;
        a.coef(0, 0) = 100.0; // A1 [мм]
        SITE_ATM40(mjd, utc, 0.0, 0.0, a, I, I, Z, dx, dv);
        double exp = 0.1 * std::cos(std::fmod(w1 * Dt, cnst::TWOPI));
        check("annual cos  dx(V)", dx(0), exp);
    }

    // --- Случай 4: скорость от годовой гармоники по востоку, A1=100 мм ---
    // dv_E = -A1*w1*sin(w1*Dt) [мм/сут] / SECDAY * 1e-3 [м/с].
    {
        AtmLoadData a;
        a.coef(1, 0) = 100.0; // East A1 [мм]
        SITE_ATM40(mjd, utc, 0.0, 0.0, a, I, I, Z, dx, dv);
        double dv_mm_day = -100.0 * w1 * std::sin(std::fmod(w1 * Dt, cnst::TWOPI));
        double exp = dv_mm_day / cnst::SECDAY * 1e-3;
        check("annual vel  dv(E)", dv(1), exp);
    }

    // --- Случай 5: цепочка VEN->ITRF->J2000 с реальными матрицами ---
    // Чистый постоянный член по вертикали -> dr_ven=(0.05,0,0) м, dv_ven=0.
    // Тогда dx = r2000*vw*(0.05,0,0), dv = dr2000_dt*vw*(0.05,0,0).
    {
        AtmLoadData a;
        a.coef(0, 6) = 50.0; // b0 = 50 мм -> 0.05 м
        // Произвольные, но воспроизводимые матрицы поворота/производной.
        Eigen::Matrix3d vw = Eigen::AngleAxisd(0.3, Eigen::Vector3d::UnitZ()).toRotationMatrix();
        Eigen::Matrix3d r2000 = Eigen::AngleAxisd(0.7, Eigen::Vector3d::UnitX()).toRotationMatrix();
        Eigen::Matrix3d dr2000; dr2000 << 0,1e-4,0, -1e-4,0,0, 0,0,0;
        SITE_ATM40(mjd, utc, 1000.0, 0.0, a, vw, r2000, dr2000, dx, dv);

        Eigen::Vector3d dr_ven(0.05, 0.0, 0.0);
        Eigen::Vector3d dx_itrf = vw * dr_ven;
        Eigen::Vector3d dx_exp = r2000 * dx_itrf;
        Eigen::Vector3d dv_exp = dr2000 * dx_itrf; // dv_ven=0 -> остаётся только вращательный член
        check("frame dx.x", dx(0), dx_exp(0));
        check("frame dx.y", dx(1), dx_exp(1));
        check("frame dx.z", dx(2), dx_exp(2));
        check("frame dv.x (rot term)", dv(0), dv_exp(0));
        check("frame dv.y (rot term)", dv(1), dv_exp(1));
        check("frame dv.z (rot term)", dv(2), dv_exp(2));
    }

    printf("=====================================================================\n");
    printf("  РЕЗУЛЬТАТ: %s (провалов: %d)\n", g_fail == 0 ? "PASS" : "FAIL", g_fail);
    return g_fail == 0 ? 0 : 1;
}
