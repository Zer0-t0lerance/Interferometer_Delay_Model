#include "functions.h"

namespace ariadna {
/**
 * Computes the fundamental arguments for nutation theory (IERS Conventions 2003)
 * and the number of Julian centuries since J2000.
 *
 * jd   Julian date at zero hours UTC (days).
 * ct   Coordinate time fraction of the coordinate time day (days).
 * cent Number of Julian centuries elapsed since J2000 (centuries).
 * f    Fundamental arguments (5 elements, arcseconds):
 *             f(0): l = mean anomaly of the Moon (L - perigee longitude).
 *             f(1): l' = mean anomaly of the Sun (L' - perigee longitude).
 *             f(2): F = mean longitude of the Moon minus Omega.
 *             f(3): D = mean elongation of the Moon from the Sun.
 *             f(4): Omega = longitude of the ascending node of the Moon's orbit.
 * fd   Time derivatives of the fundamental arguments (arcseconds/century).
 */
void fund_arg(double jd, double ct, double& cent, Eigen::VectorXd& f, Eigen::VectorXd& fd) {
    // Compute Julian centuries since J2000
    double days = (jd - cnst::JD2000) + ct;
    cent = days / cnst::JUL_CENT;
    double t = cent;
//    long double t = 7.995893223819302E-02;
//    long double t2 = t * t;
//    long double t3 = t * t2;
//    long double t4 = t * t3;

    // Initialize output vectors
    f.resize(5);
    fd.resize(5);

    // Coefficients for fundamental arguments (IERS Conventions 2003)
    const double* coeffs[5] = {cnst::F1, cnst::F2, cnst::F3, cnst::F4, cnst::F5};

    double g[5] = {0.00000000000000000000L, 0.00000000000000000000L, 0.00000000000000000000L, 0.00000000000000000000L, 0.00000000000000000000L};

// Compute arguments and their derivatives in arcseconds
    for (int i = 0; i < 5; ++i) {
        // Argument: g = c0*t^4 + c1*t^3 + c2*t^2 + c3*t + c4
//        g[i] = coeffs[i][0] * t4 + coeffs[i][1] * t3 + coeffs[i][2] * t2 + coeffs[i][3] * t + coeffs[i][4];
        g[i] = coeffs[i][4] + t * (coeffs[i][3] + t * (coeffs[i][2] + t * (coeffs[i][1] + t * coeffs[i][0])));
        f(i) = fmod(g[i], cnst::SEC360); // Modulo 360 degrees in arcseconds
        // Derivative: dg/dt = 4*c0*t^3 + 3*c1*t^2 + 2*c2*t + c3
        //fd(i) = 4.0 * coeffs[i][0] * t3 + 3.0 * coeffs[i][1] * t2 + 2.0 * coeffs[i][2] * t + coeffs[i][3];
        fd(i) = (((4.0 * coeffs[i][0] * t) + 3.0 * coeffs[i][1]) * t + 2.0 * coeffs[i][2]) * t + coeffs[i][3];
    }
}
} // namespace ariadna