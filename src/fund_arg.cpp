#include "functions.h"

namespace ariadna {

void fund_arg(double jd, double ct, double& cent, Eigen::VectorXd& f, Eigen::VectorXd& fd) {
    // Стандартное вычисление столетий J2000
    cent = (jd - cnst::JD2000 + ct) / cnst::JUL_CENT;

    double t = cent;
    double t2 = t * t;
    double t3 = t2 * t;
    double t4 = t3 * t;

    const double* coeffs[5] = { cnst::F1, cnst::F2, cnst::F3, cnst::F4, cnst::F5 };

    for (int i = 0; i < 5; ++i) {
        const double* c = coeffs[i];
        // f = c[0]*t^4 + c[1]*t^3 + c[2]*t^2 + c[3]*t + c[4]
        f(i) = std::fmod(c[0]*t4 + c[1]*t3 + c[2]*t2 + c[3]*t + c[4], cnst::SEC360);
        // fd = 4*c[0]*t^3 + 3*c[1]*t^2 + 2*c[2]*t + c[3]
        fd(i) = 4.0*c[0]*t3 + 3.0*c[1]*t2 + 2.0*c[2]*t + c[3];
    }
}
}  // namespace ariadna