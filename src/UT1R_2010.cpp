#include "functions.h"

namespace ariadna {

void ut1r_2010(const Eigen::VectorXd& f, double& dut, double& dlod, double& domega) {
    double ut = 0.0;
    double lod = 0.0;
    double omg = 0.0;

    for (int j = 0; j < cnst::N_TIDAL_TERMS; ++j) {
        double arg_arcsec = 0.0;
        for (int i = 0; i < 5; ++i) {
            arg_arcsec += cnst::TIDAL_COEFFS[j][i] * f(i);
        }

        // Стандарт IERS: аргумент в радианах
        double arg_rad = std::fmod(arg_arcsec, 1296000.0) * cnst::CARCRAD;

        ut  += cnst::TIDAL_COEFFS[j][5] * std::sin(arg_rad) + cnst::TIDAL_COEFFS[j][6] * std::cos(arg_rad);
        lod += cnst::TIDAL_COEFFS[j][7] * std::cos(arg_rad) + cnst::TIDAL_COEFFS[j][8] * std::sin(arg_rad);
        omg += cnst::TIDAL_COEFFS[j][9] * std::cos(arg_rad) + cnst::TIDAL_COEFFS[j][10] * std::sin(arg_rad);
    }

    // МАСШТАБИРОВАНИЕ (из UT1R_2010.F90.TXT)
    dut    = ut  * 1.0e-4;   // Перевод в секунды
    dlod   = lod * 1.0e-5;   // Перевод в секунды
    domega = omg * 1.0e-14;  // Перевод в рад/сек
}

} // namespace ariadna