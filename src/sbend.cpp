#include "functions.h"

namespace ariadna {
double delta(double ad1, double ad2, double bd1, double bd2, double zd2) {
    return (ad2 - ad1) * std::exp(bd1 * (zd2 - bd2));
}

double sbend(double el_rad, double temp_k, double humid_f, double press_hg) {
    double z2 = 90.0 - el_rad / cnst::CDEGRAD;
    double t2 = temp_k;
    double r = humid_f; 
    double p2 = press_hg;

    double d3 = 1.0 + delta(cnst::SBEND_Z1, z2, cnst::SBEND_C1, cnst::SBEND_C2, z2); 
    double fp = (p2 / cnst::SBEND_P1) * (1.0 - delta(cnst::SBEND_P1, p2, cnst::SBEND_A1, cnst::SBEND_A2, z2) / d3);
    double ft = (cnst::SBEND_T1 / t2) * (1.0 - delta(cnst::SBEND_T1, t2, cnst::SBEND_B1, cnst::SBEND_B2, z2) / d3);
    double fw = 1.0 + (cnst::SBEND_WP[0] * r * std::exp((cnst::SBEND_WP[1] * t2 - cnst::SBEND_WP[2]) / (t2 - cnst::SBEND_WP[3])) / (t2 * p2));

    double u = (z2 - cnst::SBEND_E[0]) / cnst::SBEND_E[1];
    double x = cnst::SBEND_E[10];
    
    // ИСПРАВЛЕНИЕ: цикл с 1 до 8
    for (int i = 1; i <= 8; ++i) {
        x = cnst::SBEND_E[10 - i] + u * x;
    }

    double sbend_value = ft * fp * fw * (std::exp(x / d3) - cnst::SBEND_E[11]);
    return sbend_value * cnst::CARCRAD;
}
} // namespace ariadna