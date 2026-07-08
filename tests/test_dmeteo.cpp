// test_dmeteo.cpp
//
// Юнит-тест DMETEO1_DT/DMETEO2_DT: метео станции заданы известными полиномами по
// времени; проверяем, что производные восстанавливаются точно.
//   P(x) = 1000 + 5 x + 2 x^2  -> dP/dt = 5 + 4 x
//   T(x) = 20 + 1 x            -> dT/dt = 1
//   H(x) = 50 - 3 x            -> dH/dt = -3
// где x = (mjd+utc) - t_mean [сут].

#include "../src/functions.h"
#include <cstdio>
#include <cmath>
#include <vector>

using namespace ariadna;

int main() {
    printf("=====================================================================\n");
    printf("  UNIT: DMETEO1_DT / DMETEO2_DT -- производные метео (МНК-полином)\n");
    printf("---------------------------------------------------------------------\n");

    const double t_mean = 58120.5;
    const int ndeg = 2;

    // 5 наблюдений станции 0 в моменты x = -0.2..0.2 сут.
    std::vector<Observation> obs;
    const double xs[5] = {-0.2, -0.1, 0.0, 0.1, 0.2};
    for (double x : xs) {
        Observation o{};
        double t = t_mean + x;
        o.mjd = (int)std::floor(t);
        o.utc = t - o.mjd;
        o.sta1 = 0; o.sta2 = 1; o.sou = 0;
        o.p1 = 1000.0 + 5.0 * x + 2.0 * x * x; // давление
        o.t1 = 20.0 + 1.0 * x;                 // температура
        o.e1 = 50.0 - 3.0 * x;                 // влажность
        o.p2 = 900.0; o.t2 = 10.0; o.e2 = 40.0;
        obs.push_back(o);
    }

    std::vector<Eigen::VectorXd> tc, pc, hc;
    dmeteo1_dt(obs, 2, t_mean, ndeg, tc, pc, hc);

    // Проверка в точке x = 0.1 (mjd/utc из t_mean+0.1).
    double xq = 0.1, tq = t_mean + xq;
    int mjd = (int)std::floor(tq); double utc = tq - mjd;
    double dTdt, dPdt, dHumdt;
    dmeteo2_dt(0, ndeg, mjd, utc, t_mean, tc, pc, hc, dTdt, dPdt, dHumdt);

    double dP_exp = 5.0 + 4.0 * xq;  // 5.4
    double dT_exp = 1.0;
    double dH_exp = -3.0;
    printf("  dPdt = %.9f (ожид %.9f)  |Δ|=%.2e\n", dPdt, dP_exp, std::fabs(dPdt - dP_exp));
    printf("  dTdt = %.9f (ожид %.9f)  |Δ|=%.2e\n", dTdt, dT_exp, std::fabs(dTdt - dT_exp));
    printf("  dHumdt = %.9f (ожид %.9f)  |Δ|=%.2e\n", dHumdt, dH_exp, std::fabs(dHumdt - dH_exp));

    // Станция без наблюдений (индекс не участвует) -> нулевые коэффициенты, rate 0.
    // (в этом тесте обе станции участвуют; проверим устойчивость: 1 точка -> rate 0)

    bool ok = std::fabs(dPdt - dP_exp) < 1e-6 &&
              std::fabs(dTdt - dT_exp) < 1e-6 &&
              std::fabs(dHumdt - dH_exp) < 1e-6;
    printf("---------------------------------------------------------------------\n");
    printf("  РЕЗУЛЬТАТ: %s\n", ok ? "PASS" : "FAIL");
    return ok ? 0 : 1;
}
