// t_eph40.cpp
//
// Координатное время T_eph (≈ TDB) как доля суток — роль подпрограммы T_EPH40:
// даёт аргумент времени для входа в эфемериды JPL DE (которые заданы в T_eph ≈ TDB).
//
// Реализовано на SOFA iauDtdb (TDB − TT, секунды) — авторитетная модель IAU:
//   ct = tt_frac + (TDB − TT) / 86400.
// Топоцентрические аргументы (долгота, расстояние от оси/экватора) дают суточный
// член TDB−TT (уровень нс); для геоцентра их можно задать нулями.

#include "functions.h"

extern "C" {
#include "sofa.h"
}

namespace ariadna {

double t_eph40(double jd, double tt_frac, double ut1_frac,
               double lon_gcen, double u_site_km, double v_site_km) {
    // iauDtdb(TT как две части, UT1-доля, вост. долгота [рад], u,v [км]) -> TDB-TT [сек].
    double dtdb = iauDtdb(jd, tt_frac, ut1_frac, lon_gcen, u_site_km, v_site_km);
    return tt_frac + dtdb / cnst::SECDAY; // доля суток в шкале TDB
}

} // namespace ariadna
