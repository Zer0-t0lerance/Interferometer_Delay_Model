#include "functions.h"

namespace ariadna {
void terms_lib(double cent, const Eigen::VectorXd& f, const Eigen::VectorXd& fd, Eigen::MatrixXd& dEOP_lib) {
    dEOP_lib.setZero(3, 2);
    if (f.size() != 5 || fd.size() != 5) return;

    double t = cent;
    double t2 = t * t;
    double t3 = t * t2;
    double dcent = 1.0 / (cnst::JUL_CENT * cnst::SECDAY);

    // Theta (arcsec): th5*t3 + th4*t2 + th3*t + th2*t + th1
    double theta = (cnst::THETA_COEFFS[4] * t3 + cnst::THETA_COEFFS[3] * t2 + cnst::THETA_COEFFS[2] * t + cnst::THETA_COEFFS[1] * t + cnst::THETA_COEFFS[0]) * 15.0 + cnst::SEC360 / 2.0;
    // dTheta (arcsec/sec): (3*th5*t2 + 2*th4*t + th3 + th2) * dcent * 15
    double dtheta = ((3.0 * cnst::THETA_COEFFS[4] * t2 + 2.0 * cnst::THETA_COEFFS[3] * t + cnst::THETA_COEFFS[2] + cnst::THETA_COEFFS[1]) * dcent) * 15.0;

    double dut = 0.0, dx = 0.0, dy = 0.0;
    double dutdt = 0.0, dxdt = 0.0, dydt = 0.0;

    for (int j = 0; j < 11; ++j) {
        double arg = 0.0, argdot = 0.0;
        for (int i = 0; i < 5; ++i) {
            arg += (cnst::amp_xy_data[j][i] * f(i));
            argdot += (cnst::amp_xy_data[j][i] * fd(i) * dcent);
        }
        arg += (cnst::amp_xy_data[j][5] * theta);
        argdot += (cnst::amp_xy_data[j][5] * dtheta);

        arg = std::fmod(arg, cnst::SEC360) * cnst::CARCRAD;
        if (arg < 0.0) arg += cnst::TWOPI;
        argdot *= cnst::CARCRAD;

        double s = std::sin(arg), c = std::cos(arg);

        // XY (microarcsec)
        dx += (cnst::amp_xy_data[j][6] * s + cnst::amp_xy_data[j][7] * c);
        dy += (cnst::amp_xy_data[j][8] * s + cnst::amp_xy_data[j][9] * c);

        // UT1 (microsec; LOD ignored)
        dut += (cnst::amp_ut_data[j][6] * s + cnst::amp_ut_data[j][7] * c);

        // Derivs
        double s_dot = c * argdot, c_dot = -s * argdot;
        dutdt += (cnst::amp_ut_data[j][6] * s_dot + cnst::amp_ut_data[j][7] * c_dot);
        dxdt += (cnst::amp_xy_data[j][6] * s_dot + cnst::amp_xy_data[j][7] * c_dot);
        dydt += (cnst::amp_xy_data[j][8] * s_dot + cnst::amp_xy_data[j][9] * c_dot);
    }

    dEOP_lib(0, 0) = dut * 1e-6; dEOP_lib(0, 1) = dutdt * 1e-6;
    dEOP_lib(1, 0) = dx * 1e-6;  dEOP_lib(1, 1) = dxdt * 1e-6;
    dEOP_lib(2, 0) = dy * 1e-6;  dEOP_lib(2, 1) = dydt * 1e-6;
}
} // namespace ariadna