#include "functions.h"

namespace ariadna {
double delta(double ad1, double ad2, double bd1, double bd2, double zd2) {
    return (ad2 - ad1) * std::exp(bd1 * (zd2 - bd2));
}

double sbend(double el_rad, double temp_k, double humid_f, double press_hg) {
    // Convert zenith angle to degrees
    double z2 = 90.0 - el_rad / cnst::CDEGRAD;
    double t2 = temp_k;
    double r = humid_f; // Fractional humidity (0.0 -> 1.0)
    double p2 = press_hg;

    // Calculate corrections
    double d3 = 1.0 + delta(cnst::SBEND_Z1, z2, cnst::SBEND_C1, cnst::SBEND_C2, z2); // Исправлено
    double fp = (p2 / cnst::SBEND_P1) * (1.0 - delta(cnst::SBEND_P1, p2, cnst::SBEND_A1, cnst::SBEND_A2, z2) / d3);
    double ft = (cnst::SBEND_T1 / t2) * (1.0 - delta(cnst::SBEND_T1, t2, cnst::SBEND_B1, cnst::SBEND_B2, z2) / d3);
    double fw = 1.0 + (cnst::SBEND_WP[0] * r * std::exp((cnst::SBEND_WP[1] * t2 - cnst::SBEND_WP[2]) / (t2 - cnst::SBEND_WP[3])) / (t2 * p2));

    // Debug output
//    std::cout << "z2: " << z2 << ", d3: " << d3 << ", fp: " << fp << ", ft: " << ft << ", fw: " << fw << std::endl;

    // Calculate optical refraction
    double u = (z2 - cnst::SBEND_E[0]) / cnst::SBEND_E[1];
    double x = cnst::SBEND_E[10];
    for (int i = 0; i <= 8; ++i) {
        x = cnst::SBEND_E[10 - i] + u * x;
    }

    // Debug output
//    std::cout << "u: " << u << ", x: " << x << std::endl;

    // Combine factors and convert from arcseconds to radians
    double sbend_value = ft * fp * fw * (std::exp(x / d3) - cnst::SBEND_E[11]);
//    std::cout << "sbend_value: " << sbend_value << std::endl;

    return sbend_value * cnst::CARCRAD;
}
} // namespace ariadna