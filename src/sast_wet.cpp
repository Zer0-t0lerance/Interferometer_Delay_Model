// sast_wet.cpp
#include "functions.h"

namespace ariadna {

void sast_wet(double rel_hum, double tc, double dot_rel_hum, double dot_tc, double& z_w, double& dot_z_w) {
    // Handle zero humidity
    if (rel_hum == 0.0) {
        rel_hum = cnst::HUMID_DEFAULT;
    }

    // Calculate saturation vapor pressure
    double temp = tc + cnst::ESAT_TEMP_OFFSET;
    double esat = cnst::ESAT_COEFF * std::exp(cnst::ESAT_EXP_COEFF * tc / temp);
    double dot_esat = esat * (cnst::ESAT_EXP_COEFF / temp - cnst::ESAT_EXP_COEFF * tc / (temp * temp)) * dot_tc;

    // Calculate zenith delay
    double work = rel_hum / cnst::PERCENT_TO_FRACTION * esat;
    temp = tc + cnst::KELVIN_OFFSET;
    z_w = cnst::SAST_WET_COEFF * (cnst::SAST_WET_TEMP_COEFF / temp + cnst::SAST_WET_CONST) * work;
    dot_z_w = -cnst::SAST_WET_COEFF * (dot_tc * cnst::SAST_WET_TEMP_COEFF / (temp * temp)) * work +
              cnst::SAST_WET_COEFF * (cnst::SAST_WET_TEMP_COEFF / temp + cnst::SAST_WET_CONST) * dot_rel_hum * esat +
              z_w * dot_esat / esat;

//    printf("SAST WET rel_hum = %.16le\n", rel_hum);
//    printf("SAST WET dot_rel_hum = %.16le\n", dot_rel_hum);
//    printf("SAST WET tc = %.16le\n", tc);
//    printf("SAST WET dot_tc = %.16le\n", dot_tc);
//    printf("SAST WET Z_w = %.16le\n", z_w);
//    printf("SAST WET dot_Z_w = %.16le\n", dot_z_w);
//    printf("SAST WET esat = %.16le\n", esat);
//    printf("SAST WET dot_esat = %.16le\n\n", dot_esat);
}

} // namespace ariadna