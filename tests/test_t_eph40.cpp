// test_t_eph40.cpp
//
// Проверка t_eph40 против эталона SOFA iauDtdb (внешняя истина).
// SOFA t_dtdb: iauDtdb(2448939.5, 0.123, 0.76543, 5.0123, 5525.242, 3190.0)
//              = -0.1280368005936998991e-2 сек (TDB - TT).
// Отсюда ct = 0.123 + (TDB-TT)/86400, и (ct - tt_frac)*86400 == TDB-TT.
//
// Сборка (из корня репозитория):
//   g++ -std=c++17 -I.\external\eigen -I.\external\sofa\20190722\c\src ^
//       src\t_eph40.cpp tests\test_t_eph40.cpp ^
//       .\external\sofa\20190722\c\build\libsofa.a -o .\tests\test_t_eph40.exe

#include "../src/functions.h"
#include <cstdio>
#include <cmath>

using namespace ariadna;

int main() {
    printf("=====================================================================\n");
    printf("  UNIT: t_eph40 (координатное время TDB, SOFA iauDtdb)\n");
    printf("---------------------------------------------------------------------\n");

    const double jd = 2448939.5, tt = 0.123, ut1 = 0.76543;
    const double lon = 5.0123, u_km = 5525.242, v_km = 3190.0;

    double ct = t_eph40(jd, tt, ut1, lon, u_km, v_km);

    // Ожидаемое ct вычисляем так же, как t_eph40 (без вычитания близких чисел,
    // чтобы не терять точность): ct = tt + (TDB-TT)/86400.
    const double ref_dtdb = -0.1280368005936998991e-2; // сек (эталон SOFA)
    double ct_expected = tt + ref_dtdb / 86400.0;
    double d = std::fabs(ct - ct_expected);
    bool ok = d < 1e-16;

    printf("  ct (доля суток TDB) = %.17f\n", ct);
    printf("  ожидаемое ct        = %.17f\n", ct_expected);
    printf("  |ct - ct_expected| = %.2e  %s\n", d, ok ? "OK" : "FAIL");

    // Санитарная проверка масштаба поправки (~1.7 мс -> ~2e-8 доли суток).
    bool scale_ok = std::fabs(ct - tt) < 5e-8;
    printf("  |ct - tt| = %.2e (ожид. < 5e-8)  %s\n", std::fabs(ct - tt), scale_ok ? "OK" : "FAIL");

    printf("=====================================================================\n");
    bool pass = ok && scale_ok;
    printf("  РЕЗУЛЬТАТ: %s\n", pass ? "PASS" : "FAIL");
    return pass ? 0 : 1;
}
