// sast_dry.cpp
#include "functions.h"

namespace ariadna {

void sast_dry(double pres, double dot_pres, double lat_geod, double height, double dpdh, double& z_d, double& dot_z_d, double& dz_ddh) {
    double f = 1.0 - cnst::GRAVITY_LAT_COEFF * std::cos(2.0 * lat_geod) - cnst::GRAVITY_HEIGHT_COEFF * height / cnst::NHMF_HEIGHT_SCALE;
    z_d = cnst::SAST_DRY_COEFF * pres / f;
    dot_z_d = z_d * dot_pres / pres;
    dz_ddh = cnst::SAST_DRY_COEFF * (dpdh - pres * cnst::GRAVITY_HEIGHT_DERIV / f) / f;

//    printf("SAST DRY pres = %.16le\n", pres);
//    printf("SAST DRY dot_pres = %.16le\n", dot_pres);
//    printf("SAST DRY lat_geod = %.16le\n", lat_geod);
//    printf("SAST DRY height = %.16le\n", height);
//    printf("SAST DRY dpdh = %.16le\n", dpdh);
//    printf("SAST DRY z_d = %.16le\n", z_d);
//    printf("SAST DRY dot_z_d = %.16le\n", dot_z_d);
//    printf("SAST DRY dz_ddh = %.16le\n\n", dz_ddh);
}

} // namespace ariadna